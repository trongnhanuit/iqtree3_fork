# `--analytical-gradients`: design notes

Reverse-mode gradients of the tree log-likelihood for every substitution-model
parameter, and the optimiser built on them. Source comments of the form
`NOTE(design x.y)` point at the numbered sections below. File map:

| File | Contents |
|---|---|
| `model/phylogradient.h/.cpp` | the gradient engine (`PhyloGradient`) |
| `model/modelparammap.h/.cpp` | unconstrained parameter vector and chain rule (`ModelParamMap`) |
| `model/gradientoptimizer.h/.cpp` | capability check, the optimiser (sections 9-12), gradient-check tool, self-test |
| `test_scripts/ag/` | identity gate, gradient test driver (suites gradcheck, oracle, threads, quality, robust), NumPy oracle, benchmark |
| `.github/workflows/analytical-gradients*.yaml` | per-push gates and the nightly long checks, called from `ci.yaml` |

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
   largest step so far *and* less than `logl_epsilon` (`--ag-no-one-percent-stop`
   drops the first condition, stopping on `logl_epsilon` alone -- faster, but a
   benchmark found it can silently lose several logL on mixture models, so it
   stays off by default);

and stops when a round improves the log-likelihood by less than
`logl_epsilon`, exactly the acceptance rule of the default loop. It ends
with the same terminal branch optimisation. Before round 1, identical `+FO`
profiles (the `+Fk` default start) are made distinct, because they are a
fixed point of every gradient method: protein models take the first `k`
profiles of the smallest C-series with at least `k` classes (up to C60);
above that, a light log-normal jitter, unless `--ag-udm-name` names an
alternative (section 11).

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

## 10. EM warm-up and cascading precision (version 1.1)

The optimiser is two phases: a cheap EM warm-up (this section), then
section 9's branch+polish loop unchanged, verbatim, as the real
optimisation. Phase 1 runs once, before phase 2, and never repeats and
never polishes: for each level in `--ag-cascade`'s schedule (`{100, 10, 1,
0.1}` (those above `logl_epsilon`) then `logl_epsilon`, if on; just
`logl_epsilon` if off, the default), a branch step (Newton, fixed at 2
iterations -- phase 1 is deliberately cheap throughout, not a ramp, since
its whole point is a rough starting point, not convergence) and then one
accept-or-revert M-step per axis listed in `--ag-em-axes` (default
`W,R,F`, subject only to that option):

* **W** mixture weights: `ModelMixture::optimizeWeights`, the EM of Wang
  et al. (2008) on the class posteriors;
* **R** site rates: `RateFree::optimizeWithEM` for `+R` (rates and
  proportions, with per-category tree scaling), otherwise the rate model's
  own optimiser (Brent on alpha and/or p_inv). Afterwards the mean rate is
  moved onto the branch lengths as the default epilogue does.
* **F** profiles of `+FO` classes: posterior class memberships (summed
  over rate categories) times the site compositions, with a pseudo-count
  of 0.5; a full step is tried, then a geometric half step;

A dedicated `--emorder` benchmark (`test_scripts/ag/bench.py`) compared
`W,R,F` against `F,W,R` on 5 randomised families (profiles, weights and,
for the GTR/GTR20 families, exchangeabilities all drawn per seed). The two
orders were indistinguishable on small mixtures (`LG+F2`, `GTR+F4`, a
2-component `MIX`), `F,W,R` was 2-5 lnL units better and 2-2.6x faster on
`LG+F4+R4`, but on the largest family tested (`GTR20+F8+I+R4`, 8 profile
classes over 20 states with linked GTR20 exchangeabilities) `W,R,F` found
an optimum 157-191 lnL units better on both seeds, at roughly 2.2-2.6x the
wall time. `W,R,F` is the default because that failure mode is worse and
the family it appears on (many linked profile classes) is the harder,
more realistic case; `--ag-em-axes F,W,R` recovers the faster order where
it is known to help. This benchmark predates the phase-1/phase-2 split
below (it ran under the old design, where each axis repeated every round
at the target level); with phase 1 now a single pass by default, axis
order sensitivity has not been re-measured under the new mechanism.

Each axis's own iteration can additionally stop once its own gain falls
below a fraction of its largest gain so far, selected per axis via
`--ag-em-stop-axes <list>` (default `W,R,F`; e.g. just `F` to narrow it, or
`""` to disable) and sized via `--ag-em-stop-frac <fraction>` (default
`0.001`, no absolute floor). This is deliberately not the joint-polish
`stopEarly` rule: the
point of the EM axes here is a cheap, rough starting point for the joint
BFGS polish, not a fully-converged answer, so over-optimising one axis
while the others are still far from correct is wasted work; a purely
relative threshold (with no `logl_epsilon` floor holding it back) lets the
axes hand off to BFGS earlier. W and R already loop internally to their
own, pre-existing convergence test (`ModelMixture::optimizeWeights`'s
per-weight 1e-4 change, `RateFree::optimizeWithEM`'s per-category
equivalent); the new test is OR'd into that loop, an addition, never a
replacement, and is gated by an explicit parameter on each function rather
than a global check, since `RateFree::optimizeWithEM` is also called by
the default (non-AG) `-optalg_qmix EM` path and must stay byte-identical
there. F had no internal loop of its own (one E-step/M-step per call), so
`emProfilesLoop` now always wraps it in the same kind of default,
always-on convergence test (max per-entry profile change < 1e-4,
mirroring W/R), with `--ag-em-stop-axes F` OR'ing in the same earlier,
relative-gain exit as W/R — consistent across all three axes now. A
dedicated `--stopfrac` benchmark (`test_scripts/ag/bench.py`) compared
default (off), 0.01 and 0.001 on 3 randomised families (simplest DNA
GTR+F4+I+G4, simple LG+F4+R4, complex GTR20+F8+I+R4, 40 taxa, 2 seeds):
0.001 never meaningfully underperformed the axes being off, and was
sometimes clearly better; 0.01 produced the single best and single fastest
results in the whole benchmark but also two real regressions (up to 542
lnL units) on two different families/seeds -- a real, seed-dependent
failure mode, not noise. `W,R,F` at `0.001` is the default as the safer
choice; `0.01` remains available as a higher-variance, occasionally
faster opt-in.

Every axis snapshots the state (theta and branch lengths), runs its
M-step, renormalises through `pack`/`unpack` (floors, mean rate 1,
`ptn_invar`) and keeps the result only if the production log-likelihood
did not decrease. The steps use IQ-TREE's own EM code where it exists; the
F step is the plan's composition heuristic, which the guard makes safe.

Phase 1's levels and per-level pass are exactly as described at the top of
this section; it is never counted against `num_param_iterations` (at most
5 passes, `{100,10,1,0.1,logl_epsilon}`, whether cascade is on or off).
Phase 2 then begins at `cur_lh` however phase 1 left it (unchanged if
`--ag-em-axes ""`) and runs section 9's loop verbatim, including its own
`num_param_iterations - 2` cap and the final best-state guard, which is
shared across both phases.

`--ag-cascade` defaults off (a single phase-1 pass at the target
precision) as the historical default from before the phase-1/phase-2
split; the coarser passes it adds are now cheap (one pass per level, not
a repeat-to-convergence loop), so this default has not been re-validated
under the new mechanism.

Floors (section 6): state frequencies at `min_state_freq` and mixture
weights at 1e-3 as in the default optimiser's bounds, `p_inv` at most the
fraction of constant sites; an entry held at its floor reports a zero
gradient in the flat direction.

**The `ndim >= 50` gate (removed).** Through Stage 4 the EM axes and
cascade only engaged for at least 50 free model parameters; a hidden
`--ag-force` flag bypassed this for testing. A dedicated `--sizegate`
benchmark (`test_scripts/ag/bench.py`) tested the gate directly: 7 small
model families (plain GTR/LG, R axis only; `MIX{GTR+FO,GTR+FO}`,
`GTR+F4`, `MIX{HKY+FO,GTR+FO}`, `LG+F2`, F/W axes; `LG+C10 -mwopt`, W
axis only) across 4 rate-heterogeneity variants, 40 taxa, comparing the
gate's default (EM off) against forcing EM on with the joint polish
unchanged either way. For R-axis-only and homogeneous-substitution
mixtures the two were identical to within 0.002 lnL at a 0-14% time cost
when forced — the gate bought nothing there, consistent with W and F
being no-ops without a mixture and R already being fully covered by the
joint gradient polish. But on `LG+F2` one seed showed the gate landing
500+ log-likelihood units short of the optimum that forcing EM reached
(a genuinely different, worse fit, not noise); `MIX{HKY+FO,GTR+FO}`
showed swings of up to +-37 lnL in both directions, suggesting a rougher
landscape rather than a clean EM benefit there. Since the gate could
silently hide a large quality loss on exactly the small profile mixtures
where symmetry-breaking matters, and bought nothing measurable anywhere
else, both the gate and `--ag-force` were removed rather than kept or
retuned to a different threshold: the EM axes and cascade now engage for
every model, subject only to `--ag-em-axes`/`--ag-cascade`. `--ag-force`
is no longer a recognised option.

Phase 1 is skipped at the final, post-search model optimisation
("FINALIZING TREE SEARCH"): the topology is fixed by then and every
preceding in-search reopt already ran a full phase 1 on this same best
tree, so a further warm-up is redundant. That one call site sets
`ModelFactory::ag_skip_em` through an RAII `AgSkipEmScope`, propagated to
every partition's own factory under `-Q`/`-S` the same way `setCheckpoint`
is; it is never checkpointed, and `--ag-stats` reports `em_levels=0` there.

## 11. Start points: warm, cold, multi-start (version 1.1)

Start-point work happens once per `ModelFactory`, on the main
optimisation only (the one that prints progress); ModelFinder candidates,
NNI refits and +I+G restarts keep their current parameters.

* **warm** (default): the parameters as IQ-TREE initialised them, with
  identical `+FO` profiles made distinct (section 9). Above C60, the
  default is jitter with a warning; `--ag-udm-name NAME` uses `NAME`'s
  components instead, loaded via `-mdef` (trying `NAME_C####` and
  `NAMEpi#`, so an unmodified UDM file from Schrempf et al. 2020,
  https://github.com/dschrempf/EDCluster, works as-is). IQ-TREE bundles
  no UDM data itself, avoiding the GPLv3 vs. GPLv2-or-later question.
  An insufficient `NAME` is an error, not a silent fallback. By default
  the K highest-weighted components are borrowed (`--ag-warm-select
  weight`); `index` borrows the first K by declared order instead;
  `sample` draws K distinct components without replacement, each
  proportional to its published weight (reject-and-redraw against the
  fixed cumulative distribution, equivalent to renormalizing over the
  remaining pool at every step, but with no rebuild) — when the family
  has exactly K components, `sample` and `weight` coincide (nothing to
  sample). Weights are read from the family's own composite entry
  (`model C10 = ...FMIX{C10pi1:1:w,...}`, or `frequency NAME =
  FMIX{NAME_C0000:1:w,...}` for a UDM file); missing weights are an
  error for `weight` and `sample` alike.
* **cold** (`--ag-start cold`): every `+FO` profile is drawn as
  `pi_k proportional to pi_emp,k * exp(0.9 z_k)`, `z ~ N(0,1)`, around the
  empirical frequencies, with equal mixture weights.
* **multi-start** (`--ag-multistart N`; `-1` selects `min(100, 10 K)` for
  `K >= 2` estimated profiles; default 0, off): `N` candidates are made from the
  start by jittering the profile entries in log-ratio space, half heavily
  (0.9) and half lightly (0.01); each costs one likelihood evaluation. The
  five best and five random heavy candidates are refined with EM rounds
  (a weight step, then `r` profile steps for each `r` in `--ag-em-ratios`,
  then a rate step), sharing `--ag-multistart-budget` likelihood
  evaluations, and the best refined candidate becomes the start of the
  optimisation. All random numbers come from a private generator seeded
  from `-seed`, so the global stream and the default path are untouched.

Multi-start is off by default because it did not pay on the LG+F4+R4
benchmark (fixed tree, same binary): warm start -4982.7 in 31 s,
multi-start with 40 candidates -4986.8 in 68 s, cold start -5000.6, cold
plus multi-start -4992.2. Ranking candidates by their initial likelihood
(before any optimisation) compares noise, and the EM refinement then
commits to a basin that the polish cannot leave; the C-series warm start
already provides the diversity the profiles need. On a simulated
two-profile mixture all three starts reach the same optimum.

## 12. Integration: checkpoints, partitions, faults (version 1.1)

* **Checkpoint.** The only optimiser state that must survive a restart is
  whether the start-point work has been done; `ModelFactory` writes it as
  `ag_init_done` inside its own checkpoint struct, and only when the flag
  is set, so the default path's checkpoint stays byte-identical. A
  checkpoint written by either path is read by the other (the missing key
  simply reads as false). `--ag-abort-after init|round:N` writes the tree
  and model checkpoint at that phase without marking the model
  optimisation finished and exits; a rerun without `-redo` restores the
  parameters and continues to the same optimum. Nothing is written from
  inside an OpenMP region: the checkpoint map is shared between
  partitions.
* **Partitions.** `-Q`/`-S` models enter the hook once per partition,
  concurrently when threads allow. The optimiser keeps no static state,
  prints through a critical section and skips checkpoint writes while in
  parallel; repeated 4-thread runs give the same log-likelihood.
* **Other callers.** ModelFinder candidates and the PMSF guide-tree fit
  use the optimiser like any other call; PMSF's site-specific second pass
  and every other unsupported configuration fall back with one NOTE.
  +I+G restarts (`--opt-gamma-inv`) run through the same hook.
* **Faults.** If anything throws inside the outside pass, the RAII guards
  restore every neighbour and the engine restores the tree's root
  bookkeeping before rethrowing; the optimiser catches the exception,
  reports it once and uses finite differences for that step. The
  environment variable `AG_TEST_FAULT_EDGE=k` injects such a throw at the
  k-th edge for the test.

## 13. Benchmark (`test_scripts/ag/bench.py`)

    bench.py <iqtree3> <out_dir> [--quick | --full] [--big] [-j N] [--threads 1,8] [--repeats N] [--seeds ...]

Arms separate the starting point from the optimiser so that speed and
quality are not conflated: `old-default`, `old-init2` (`-init_nucl_freq 2`),
`old-c10warm` (the default optimiser started from C10 profiles with
`-mfopt`), `old-em` (`-optalg_qmix EM`), `new-warm` (the flag), `new-cold`
(`--ag-start cold`) and `new-multi` (`--ag-multistart -1`). Datasets are
AliSim simulations on random Yule-Harding trees with the truth recorded
(DNA GTR+FO+G4, LG and WAG four-profile mixtures with free rates, C10 with
free rates, GTR20 with C10, and with `--big` 100 taxa x 20000 sites and an
F60 stress case) plus the repository's real alignments. Simulated data are
fitted on the true tree (`-te`) and, with `--full`, also with a tree
search; real data always with a search. Every run records the final
log-likelihood, wall time, peak memory, tree length, the optimiser's
evaluation counts and, for simulated mixtures, the RMSE of the recovered
weights and profiles after matching the classes to the truth.
`results.tsv` holds one row per run; `bench_report.py` turns it into
`summary.md` (and plots when matplotlib is available). `--quick` runs
three datasets with three arms at one thread; `--full` adds the rest with
two thread counts and three timing repeats.

## 14. Command-line summary

| Option | Meaning |
|---|---|
| `--analytical-gradients` | use the analytic optimiser for supported models (others fall back with a NOTE) |
| `--ag-start warm|cold` | C-series warm start (default) or random cold start of the profiles |
| `--ag-udm-name NAME` | above C60, warm-start from NAME's profiles (via `-mdef`) instead of jitter |
| `--ag-warm-select index|weight|sample` | borrow the first K, (default) the K highest-weighted, or K weighted-sampled reference profiles |
| `--ag-multistart N` | N jittered start candidates, refined by EM; 0 (default) off, -1 automatic |
| `--ag-em-axes W,R,F` | EM axes per round (weights, rates, profiles); empty disables EM |
| `--ag-em-stop-axes <list>` | per-axis (subset of W,R,F) relative-gain early exit; default W,R,F |
| `--ag-em-stop-frac <fraction>` | relative threshold for --ag-em-stop-axes; default 0.001 |
| `--ag-cascade on|off` | phase-1 EM warm-up at coarser precisions first, not just the target (default off) |
| `--ag-optalg BFGS|LBFGSB` | driver of the joint polish |
| `--ag-no-one-percent-stop` | joint polish stops on `logl_epsilon` alone, dropping the 1%-of-largest-step condition (default off: the 1% condition stays on) |
| `--ag-stats` | print evaluation counts, EM steps and reverts |
| `--ag-gradient-check [tol]`, `--ag-gradient-check-every k`, `--ag-gradient-check-strict` | check analytic against numerical gradients during the search |
| `--ag-gradient-check-only`, `--ag-dump-gradient`, `--ag-selftest` | one-shot checks at the initial point (tests) |
| `--ag-abort-after init|round:N` | write the checkpoint and exit at that phase (resume tests) |
