# `--analytical-gradients`: design notes

Reverse-mode gradients of the tree log-likelihood for every substitution-model
parameter, and the optimiser built on them. Source comments of the form
`NOTE(design x.y)` point at the numbered sections below. File map:

| File | Contents |
|---|---|
| `model/phylogradient.h/.cpp` | the gradient engine (`PhyloGradient`) |
| `model/modelparammap.h/.cpp` | unconstrained parameter vector and chain rule (`ModelParamMap`) |
| `model/gradientoptimizer.h/.cpp` | capability check, gradient-check tool, self-test, later the optimiser |
| `test_scripts/ag/` | identity gate, gradient test driver, NumPy oracle |

Without `--analytical-gradients` none of this code runs; the default
optimiser is byte-for-byte the previous one (checked by
`test_scripts/ag/identity.sh` in CI).

## 1. Notation

* `S` states, `K` mixture components `m`, `C` rate categories `r`; a *class*
  `c = (m, r)` has weight `omega_c = prop_r * w_m` and rate `rate_r`.
* Component `m` has exchangeabilities `R_m` (symmetric, zero diagonal) and
  frequencies `pi_m`. Its rate matrix as the likelihood kernel sees it is
  `Q_m = nu_m R_m Pi_m` off the diagonal, rows summing to zero, with
  `nu_m = T_m / Z_m`, `Z_m = pi_m^T R_m pi_m` and `T_m` the requested mean
  rate (`total_num_subst`, 1 unless the mixture rescales it). States with
  `pi <= ZERO_FREQ = 1e-10` are dropped from `Q` and `pi` is renormalised
  over the remaining ones, exactly as `ModelMarkov::decomposeRateMatrix`.
* Eigendecomposition `Q = U Lambda U^-1`. For a reversible model
  `U^-1 = U^T Pi` (the code stores `evec = Pi^-1/2 H`, `inv_evec = H^T Pi^1/2`).
* On an edge `e` of length `t_e` and class `c`, `tau = rate_r * t_e`.
* For pattern `ptn` with count `f`: `L = sum_c omega_c L_c + ptn_invar`,
  `ptn_invar = p_inv * sum_{x in const states} pibar_x`, where `pibar` is
  the weight-averaged frequency vector (`ModelSubst::getStateFrequency(-1)`).

## 2. Why the kernel's partials are all we need

After `computeLikelihood()` every internal node holds one *inside* partial
oriented toward the root branch `(u, v)`. In the kernel's eigen coordinates
the inside partial of the subtree below edge `e` is `p~ = U^-1 p`, where `p`
is the ordinary conditional likelihood vector (no `pi`).

For any edge `e = (n -> c)` the tree likelihood can be written with the root
at `n`:

    L_c = o^T Pi e^{Q tau} p                                          (2.1)

where `o` is the *outside* partial (conditional likelihood of everything on
`n`'s side, again without `pi`). With `U^-1 = U^T Pi`:

    o^T Pi U e^{Lambda tau} U^-1 p = (U^-1 o)^T e^{Lambda tau} (U^-1 p) = o~^T e^{Lambda tau} p~   (2.2)

so the kernel's per-edge scalar product of two eigen-space partials *is*
(2.1). The engine therefore never leaves the kernel's representation: it
computes `o~` for every edge by attaching a private buffer to the reverse
neighbour `c->findNeighbor(n)` and calling `computeLikelihoodBranch`, which
fills it with the kernel's own SIMD code and scaling, and which returns the
tree log-likelihood evaluated on that edge as a free self-check
(`max_edge_logl_diff`).

The attach must happen *before* the traversal is built: `reorientPartialLh`
returns early when a neighbour already owns a partial, so the kernel never
steals an inside slot. The RAII guard restores the neighbour on every exit
path, and `theta_computed` is cleared afterwards because the kernel sets it
per edge.

## 3. Branch lengths and rates

From (2.2), with `s_i = e^{lambda_i tau}` and `a_i = o~_i p~_i`:

    dL_c/dt_e = rate_r * sum_i lambda_i s_i a_i,     dL_c/drate_r = t_e * sum_i lambda_i s_i a_i

and `dlogL/dt_e = sum_ptn (f/L) sum_c omega_c f_c dL_c/dt_e`, where `f_c` is
the kernel's per-class rescaling (section 5.3). `dlogL/drate_r` is summed
over all edges and components. Branch lengths are optimised by the existing
Newton scheme; their gradient is computed for the gradient check and the
identities only.

## 4. The rate matrix: divided differences

For a perturbation `dQ` of a diagonalisable `Q` (Daleckii-Krein):

    d(e^{Q tau}) = U [ (U^-1 dQ U) o X ] U^-1,    X_jk = (e^{lambda_j tau} - e^{lambda_k tau}) / (lambda_j - lambda_k)   (4.1)

with `X_jj = tau e^{lambda_j tau}` (the limit). Substituting in (2.1):

    dL_c = o~^T [ (U^-1 dQ U) o X ] p~ = sum_jk G_jk (U^-1 dQ U)_jk,     G = X o (o~ p~^T)    (4.2)

and, collecting `dQ_ab`:

    dlogL/dQ_ab = sum_j (U^-1)_ja sum_k G_jk U_bk,   i.e.   D = U^-T G U^T                     (4.3)

The engine accumulates, per edge and class, the rank-one matrix
`A_c = sum_ptn (f/L) omega_c f_c o~ p~^T` in the pattern loop (this is the
only `S^2` work per pattern, and it is skipped for components that need no
`Q` gradient, such as fixed C-series profiles), then once per edge and class
`G_m += X(Lambda_m, tau) o A_c`, and once per gradient (4.3).

### 4.1 Rooting and what `D` means

Every edge is accumulated with the outside partial on the side that
contains `u`, the parent endpoint of the root branch, so `D` is the
derivative of the likelihood *rooted at `u`* with all `S^2` entries of `Q`
treated as independent. For perturbations that keep `Pi Q` symmetric (every
exchangeability or frequency change, the diagonal) the likelihood is
root-independent and `D` is simply `dlogL/dQ`; a single off-diagonal entry
perturbed alone breaks reversibility and then the root matters. The oracle
roots its explicit-`Q` likelihood at `u` (`#root_side` in the dump) for
exactly this reason. Expected counts follow as `E[N_ab] = Q_ab D_ab`,
`E[T_a] = D_aa`.

### 4.2 The three regimes of `X` (`PhyloGradient::xKernel`)

`X_jk` as written suffers catastrophic cancellation when the eigenvalue gap
`d = lambda_j - lambda_k` is small compared with the eigenvalues. Writing
`e^{lambda_j tau} = e^{lambda_k tau} e^{d tau}`:

    X_jk = e^{lambda_k tau} * (e^{d tau} - 1) / d = e^{lambda_k tau} * expm1(tau d) / d       (4.4)

which is algebraically identical and, evaluated with `expm1`, accurate to
working precision for *any* non-zero `d`, however small. The regimes:

1. `|d| < 1e-200` (exactly coincident eigenvalues, e.g. JC/F81, or a gap that
   would underflow): `X_jj = tau e^{lambda_j tau}`. A wider window would
   cost relative accuracy `tau |d| / 2`, which is why the window is not
   `1e-10` as one might first choose.
2. `|d tau| < 0.1`: (4.4). Its own rounding error is about `1e-16` relative;
   the naive formula at this gap size would lose up to `log10(1/(d tau))`
   digits.
3. otherwise the direct formula; at `|d tau| >= 0.1` it loses at most one
   digit.

`--ag-selftest` evaluates the kernel against a long-double reference for
gaps `{0, 1e-12, 1e-6, 0.05/tau, 3}` and also shows the naive formula
failing at `d = 1e-12`.

## 5. Chain rule to the natural parameters (`ModelParamMap`)

Let `<D, Q> = sum_ab D_ab Q_ab` and write `pi` for the renormalised
frequencies. Because `Q = nu R Pi` with `nu = T / (pi^T R pi)`:

    dlogL/dR_ab = nu (D_ab pi_b + D_ba pi_a - D_aa pi_b - D_bb pi_a) - (2 pi_a pi_b / Z) <D, Q>       (5.1)

    dlogL/dpi_k  = nu sum_{a != k} R_ak (D_ak - D_aa) - (2 (R pi)_k / Z) <D, Q>
                 + root_term_k + w_m * inv_state_sum_k                                                 (5.2)

The first line of (5.2) is the dependence through `Q`; `root_term` is the
explicit `Pi` at the root in (2.1), accumulated at the root branch as
`sum_ptn (f/L) sum_r omega f (U o~)_k (U e^{Lambda tau} p~)_k` (note
`U o~ = o` and `U e^{Lambda tau} p~ = e^{Q tau} p`), and `inv_state_sum_k
= sum over constant patterns containing k of (f/L) p_inv` is the
invariant-site term. A raw (unnormalised) frequency entry then has
`dlogL/dpi_k^raw = (1/sum pi) (g_k - sum_j pi_j g_j)` with `g` from (5.2).

Class sums, taken once at the root branch: `dlogL/domega_c = sum_ptn (f/L)
f_c L_c`. Hence

    dlogL/dw_m   = sum_r prop_r dlogL/domega_(m,r) + sum_k pi_mk inv_state_sum_k
    dlogL/dprop_r = sum_m w_m dlogL/domega_(m,r)
    dlogL/dp_inv  = sum_ptn (f/L) ptn_invar / p_inv + sum_c [ dlogL/drate_c * rate_c - dlogL/dprop_c * prop_c ] / (1 - p_inv)

the last because `+I` scales every category rate by `1/(1-p)` and every
proportion by `(1-p)` (`RateGammaInvar::setPInvar`, `RateFreeInvar`).

`dlogL/dalpha = sum_c dlogL/drate_c * drate_c/dalpha`, with the Jacobian by
Richardson-extrapolated central differences of the gamma quantiles.
`RateGamma::computeRates` rescales relative to the rates it currently
stores, so the stored rates are restored before *each* evaluation.

The closed form (5.1)-(5.2) is cross-checked at every gradient check against
a chain rule that differentiates the `Q` builder itself by central
differences (`gradientFDQ`, reported as `qchain`), and against two zero-cost
identities:

* `sum_m <D_m, Q_m> = sum_e t_e dlogL/dt_e` (scaling every `Q` equals scaling every branch),
* `sum_m sum_k pi_mk root_term_mk = sum_c omega_c dlogL/domega_c`.

## 6. The unconstrained vector `theta`

| Block | Variables | Map |
|---|---|---|
| S | one per free rate variable of a component, in the component's own layout (DNA rate groups via `ModelDNA::param_spec`, GTR20's 189 entries with the last pinned at 1); one shared block when `--link-exchange-rates` / `+Fk` links a GTR-type matrix | `log(variable)` |
| F | per `+FO` component, one per state except the pinned largest one and states at or below `ZERO_FREQ` | `log(pi_k / pi_pinned)`, renormalised on unpack |
| W | `K-1` mixture weights | `log(w_m / w_{K-1})` |
| R | free-rate proportions `z` and rates `rho`, `K-1` each | `s = softmax(z)`, `prop = (1-p) s`, `rate_c = e^{rho_c} / ((1-p) sum_j s_j e^{rho_j})`, `rho_{K-1} = 0`, so the mean rate is 1 identically and the rate/branch-length scale redundancy is gone |
| alpha | gamma shape | `log(alpha)` |
| pinv | invariant proportion | `logit(p)` |

`unpack` writes through public setters, then `decomposeRateMatrix`,
`computePtnInvar` and `clearAllPartialLH`. Which rate entries a variable
drives is discovered numerically at construction (perturb one variable,
see which entries move; the mapping is a plain copy, so this is exact).

Jacobians used by `ModelParamMap::naturalToTheta`, with `g_r = dlogL/drate`,
`g_p = dlogL/dprop`:

    dlogL/drho_k = g_r[k] r_k - (1-p) s_k r_k sum_c g_r[c] r_c
    dlogL/dz_k   = (1-p) s_k (g_p[k] - sum_c g_p[c] s_c) - s_k ((1-p) r_k - 1) sum_c g_r[c] r_c
    dlogL/dtheta_W,k = w_k (dlogL/dw_k - sum_j w_j dlogL/dw_j)          (softmax)
    dlogL/dlog alpha = alpha dlogL/dalpha,   dlogL/dlogit p = p (1-p) dlogL/dp

## 7. Worked example: F81 on one edge

Two taxa, states `A` and `C` observed, `pi = (0.4, 0.3, 0.2, 0.1)`, F81 so
`R_ab = 1`, `Z = 1 - sum pi^2 = 0.7`, `nu = 1/0.7`. Eigenvalues are
`{0, -nu, -nu, -nu}`: three coincide, so `X` uses regime 1 on those
diagonal blocks and regime 2/3 between `0` and `-nu`; the naive formula
would divide `0/0`. With `P(t) = e^{Q t}`, `L = pi_A P_AC(t)` and
`P_AC(t) = pi_C (1 - e^{-nu t})`. Differentiating `log L` by `t` gives
`nu e^{-nu t} / (1 - e^{-nu t})`, which is what the engine's
`dlogL/dt` returns, and by `pi_C` (holding the normalisation, so through
`nu` as well) gives the combination of the `Q` term in (5.2) and the
root term with `root_term_C = 0` (the root state is `A`). The
`--ag-gradient-check-only` rows for `F81+F` on `example/example.phy` show
these agree with central differences to about `1e-9`.

## 8. Gradient check output

`--ag-gradient-check-only [--ag-gradient-check-tol 1e-4]` writes
`<prefix>.gradcheck.tsv` with one row per branch and per `theta` entry:

    iter level idx name x analytic numeric fd_h fd_h2 newton abs_err rel_err status

`numeric` is the Richardson extrapolation of the two central differences
`fd_h` (step `h = 1e-4 max(1,|x|)`) and `fd_h2` (step `h/2`); branch rows
also carry the tree's Newton derivative. `rel_err =
|a - n| / max(|a|, |n|, 1e-6 ||g||_inf)`. One summary line is printed:

    GRADCHECK iter=0 n=41 max_rel=5.1e-06 worst=pinv n_fail=0 edge_lnl_check=PASS qchain=1.9e-10 q_identity=3.3e-15 root_identity=0.0e+00 identities=PASS

`--ag-dump-gradient` additionally writes `<prefix>.aggrad.tsv` with the
model, tree, `dlogL/dt`, `dlogL/dQ` per component and the natural
gradients, which `test_scripts/ag/oracle.py` checks against an independent
NumPy implementation.

## 9. The optimiser (`GradientOptimizer`), version 1

`ModelFactory::optimizeParameters` keeps its prologue (initial likelihood,
`mlInitial`) and epilogue (rate rescaling, root position, `writeInfo`,
"took N rounds") for both paths; only the alternating loop in the middle is
replaced when `--analytical-gradients` is set and `supports()` accepts the
model (section 3.2 of the plan). `supports()` refuses, with a one-line
reason printed once, everything the engine does not cover: partition or
tree-mixture containers, heterotachy, fused `*G`/`*R`, codon, PoMo, DNA
error models, non-reversible kernels, site-specific models or rates,
`+ASC`, `-mem`, MPI with several processes, models without free parameters,
and the special DNA frequency parametrisations. Edge-linked partition
models (`-p`/`-q`) never reach the hook; the option parser clears the flag
for them with a warning. Per-partition models (`-Q`/`-S`) enter the hook
once per partition, possibly concurrently, so the optimiser keeps no global
state and prints inside a critical section.

One `optimize()` call runs, for round `k = 1, 2, ...` up to the
`num_param_iterations` cap:

1. the branch-length step of the default loop (Newton on all branches, tree
   length scaling, or nothing, according to `fixed_len`);
2. one joint BFGS minimisation of `-logL(theta)` over the whole vector of
   section 6, with the engine's gradient (`derivativeFunk`) and the
   existing `dfpmin` driver (`--ag-optalg LBFGSB` uses L-BFGS-B with 20
   retained updates instead). The bounds are numerical fences only, and
   `restartParameters` never randomises. Inside one minimisation the
   `stopEarly` hook ends the search once a step gains less than 1% of the
   largest step so far *and* less than `logl_epsilon`;

and stops when a round improves the log-likelihood by less than
`logl_epsilon`, exactly the acceptance rule of the default loop. It ends
with the same terminal branch optimisation. Before round 1, identical `+FO`
profiles (the `+Fk` default start) are made distinct, because they are a
fixed point of every gradient method: protein models take the first `k`
profiles of the smallest C-series with at least `k` classes, other data a
light log-normal jitter from a private seeded generator.

Safety nets: a non-finite likelihood inside a line search returns a huge
value (a rejected step, never a crash); a gradient with a non-finite entry
or an eigendecomposition failing `||U U^-1 - I|| < 1e-8` falls back to the
base class's finite differences for that step, warning once; and the whole
call ends by restoring the best state seen (parameters and branch lengths)
if the final log-likelihood is below both the entry score and the best
round, because `IQTree::optimizeModelParameters` aborts on a regression
larger than 1. `--ag-stats` prints the counts of likelihood and gradient
evaluations; `--ag-gradient-check [tol]` re-checks every parameter against
central differences at every `k`-th gradient evaluation
(`--ag-gradient-check-every k`) and appends the rows to
`<prefix>.gradcheck.tsv`, and `--ag-gradient-check-strict` turns a
disagreement into an error.

Not in version 1 (later stages): per-axis EM warm start, cascading
precision, multi-start, checkpointing of the optimiser's own state.

## 10. EM axes and cascading precision (version 1.1)

For vectors of at least 50 parameters (or with `--ag-force`) each round of
section 9 also runs, before the polish, one accept-or-revert M-step per
axis listed in `--ag-em-axes` (default `W,R,F`):

* **W** mixture weights: `ModelMixture::optimizeWeights`, the EM of Wang
  et al. (2008) on the class posteriors;
* **R** site rates: `RateFree::optimizeWithEM` for `+R` (rates and
  proportions, with per-category tree scaling), otherwise the rate model's
  own optimiser (Brent on alpha and/or p_inv). Afterwards the mean rate is
  moved onto the branch lengths as the default epilogue does;
* **F** profiles of `+FO` classes: posterior class memberships (summed
  over rate categories) times the site compositions, with a pseudo-count
  of 0.5; a full step is tried, then a geometric half step.

Every axis snapshots the state (theta and branch lengths), runs its
M-step, renormalises through `pack`/`unpack` (floors, mean rate 1,
`ptn_invar`) and keeps the result only if the production log-likelihood
did not decrease. The steps use IQ-TREE's own EM code where it exists; the
F step is the plan's composition heuristic, which the guard makes safe.

With `--ag-cascade on` (default off) the rounds run at decreasing precision
levels `{100, 10, 1, 0.1}` (those above `logl_epsilon`) and then at
`logl_epsilon`. At a coarse level a round is the branch step plus the EM
axes plus, with `--ag-polish per-level` (the default), the polish, and the
level ends when a round gains less than `eps` or, from the second round on,
less than 1% of the best round at that level. The target level runs the
full round of section 9 with the default loop's rule. All levels share the
`num_param_iterations` cap and the final best-state guard. The cascade is
off by default because on the LG+F4+R4 benchmark it reached the same
optimum as the single level (-4982.65 vs -4982.66) in 40 s instead of 30 s,
and without the per-level polish the EM-only coarse levels steered the
search to a worse basin (-4988.7); dropping the F axis costs 13-21
log-likelihood units on that benchmark, so all three axes stay on.

Floors (section 6): state frequencies at `min_state_freq` and mixture
weights at 1e-3 as in the default optimiser's bounds, `p_inv` at most the
fraction of constant sites; an entry held at its floor reports a zero
gradient in the flat direction.
