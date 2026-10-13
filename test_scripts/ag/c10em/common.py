"""
Shared helpers for the C10 re-estimation pipeline (map-reduce EM over many
alignments; design doc section 16). Standard library + NumPy only.

Model files are IQ-TREE -mdef NEXUS files written by write_mdef():
  frequency <NAME>pi<k> = <20 values>;
  model <NAME> = POISSON+G4+FMIX{<NAME>pi1:<rate>:<weight>,...};
with %.17g numbers, no whitespace inside FMIX{} and weights/rates normalised
so that sum(w) = 1 and sum(w * rate) = 1 exactly as written.
"""
import os
import re
import math
import hashlib
import numpy as np

AA = "ARNDCQEGHILKMFPSTWYV"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))


# ---------------- alignments ----------------

def read_alignment(path):
    """FASTA or PHYLIP (sequential or interleaved) -> (names, seqs dict)."""
    with open(path) as f:
        text = f.read()
    lines = [l.rstrip("\r\n") for l in text.splitlines() if l.strip()]
    if lines and lines[0].startswith(">"):
        names, seqs, cur = [], {}, None
        for l in lines:
            if l.startswith(">"):
                cur = l[1:].strip().split()[0]
                names.append(cur)
                seqs[cur] = []
            else:
                seqs[cur].append(l.strip().replace(" ", ""))
        return names, {n: "".join(s).upper() for n, s in seqs.items()}
    ntax, nsite = [int(x) for x in lines[0].split()[:2]]
    names, seqs = [], {}
    for l in lines[1:]:
        parts = l.split(None, 1)
        if len(names) < ntax and len(parts) == 2:
            names.append(parts[0])
            seqs[parts[0]] = parts[1].replace(" ", "").upper()
        else:   # interleaved continuation: append to the first sequence still short
            for n in names:
                if len(seqs[n]) < nsite:
                    seqs[n] += l.replace(" ", "").upper()
                    break
    for n in names:
        if len(seqs[n]) != nsite:
            raise ValueError("%s: sequence %s has length %d, expected %d" % (path, n, len(seqs[n]), nsite))
    return names, seqs


def write_phylip(path, names, seqs):
    with open(path, "w") as f:
        f.write("%d %d\n" % (len(names), len(seqs[names[0]])))
        for n in names:
            f.write("%s %s\n" % (n, seqs[n]))


def site_counts(seqs):
    """per site: counts of the 20 amino acids (gaps/ambiguities ignored), shape (nsite, 20)"""
    idx = {a: i for i, a in enumerate(AA)}
    nsite = len(next(iter(seqs.values())))
    out = np.zeros((nsite, 20))
    for s in seqs.values():
        for j, ch in enumerate(s):
            i = idx.get(ch)
            if i is not None:
                out[j, i] += 1
    return out


# ---------------- trees ----------------

class Node:
    __slots__ = ("name", "children", "length")

    def __init__(self, name=None, length=None):
        self.name, self.children, self.length = name, [], length


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                n.children.append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                pos += 1   # ')'
                break
        start = pos
        while pos < len(s) and s[pos] not in ":,)(":
            pos += 1
        n.name = s[start:pos].strip().strip("'") or None
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",)(":
                pos += 1
            n.length = float(s[start:pos])
        return n
    return node()


def write_newick(n, top=True):
    if n.children:
        out = "(" + ",".join(write_newick(c, False) for c in n.children) + ")"
        if n.name:
            out += n.name
    else:
        out = n.name
    if not top and n.length is not None:
        out += ":%.10g" % n.length
    return out + (";" if top else "")


def leaves(n):
    return [n.name] if not n.children else [x for c in n.children for x in leaves(c)]


def prune(n, keep):
    """drop leaves not in `keep`, collapse unary nodes (lengths add); returns the new root or None"""
    if not n.children:
        return n if n.name in keep else None
    kids = [k for k in (prune(c, keep) for c in n.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        k = kids[0]
        if k.length is not None and n.length is not None:
            k.length += n.length
        elif n.length is not None:
            k.length = n.length
        return k
    n.children = kids
    return n


def scale_tree(s, c):
    """multiply every branch length in a newick string by c"""
    return re.sub(r":([0-9.eE+-]+)", lambda m: ":%.10g" % (float(m.group(1)) * c), s)


# ---------------- model files ----------------

def phi(pi):
    """POISSON event-rate factor: expected substitutions per redraw, 1 - sum(pi^2)"""
    return 1.0 - float(np.sum(np.asarray(pi) ** 2))


def fmt(x):
    s = "%.17g" % x
    assert "e+" not in s and "E+" not in s, s
    return s


def write_mdef(path, name, profiles, weights, rates):
    """write an -mdef file; weights and rates are renormalised (sum w = 1, sum w*r = 1)"""
    profiles = [np.asarray(p, float) / np.sum(p) for p in profiles]
    w = np.asarray(weights, float)
    w = w / w.sum()
    r = np.asarray(rates, float)
    r = r / float(np.dot(w, r))
    lines = ["#nexus", "begin models;"]
    for k, p in enumerate(profiles):
        lines.append("frequency %spi%d = %s;" % (name, k + 1, " ".join(fmt(x) for x in p)))
    terms = ",".join("%spi%d:%s:%s" % (name, k + 1, fmt(r[k]), fmt(w[k])) for k in range(len(profiles)))
    lines.append("model %s = POISSON+G4+FMIX{%s};" % (name, terms))
    lines.append("end;")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.replace(tmp, path)


def read_mdef(path, name=None):
    """-> (name, profiles list, weights array, rates array) for a POISSON+G4+FMIX model file"""
    freqs, model = {}, None
    with open(path) as f:
        for l in f:
            l = l.strip()
            m = re.match(r"frequency\s+(\S+)\s*=\s*([0-9.eE+\-\s]+);", l)
            if m:
                freqs[m.group(1)] = np.array([float(x) for x in m.group(2).split()])
                continue
            m = re.match(r"model\s+(\S+)\s*=\s*(.+);", l)
            if m and (name is None or m.group(1) == name):
                model = (m.group(1), m.group(2))
    if model is None:
        raise ValueError("no model %s in %s" % (name or "", path))
    body = re.search(r"FMIX\{([^}]*)\}", model[1]).group(1)
    profiles, w, r = [], [], []
    for term in body.split(","):
        fname, rate, weight = term.split(":")
        profiles.append(freqs[fname])
        r.append(float(rate))
        w.append(float(weight))
    return model[0], profiles, np.array(w), np.array(r)


def published_c10():
    """the published C10 profiles and weights from example/models.nex"""
    _, profiles, w, _ = read_mdef(os.path.join(REPO, "example", "models.nex"), "C10")
    return profiles, w


def paper_rates(profiles):
    return np.array([phi(p) for p in profiles])


# ---------------- --ag-dump-em output ----------------

def read_agem(path):
    """parse <prefix>.agem.tsv; returns None if the file is missing or has no #END"""
    if not os.path.exists(path):
        return None
    d = {"classes": []}
    end = False
    with open(path) as f:
        for l in f:
            l = l.rstrip("\n")
            if not l.startswith("#"):
                continue
            key, _, rest = l[1:].partition("\t")
            v = rest.split("\t") if rest else []
            if key == "END":
                end = True
            elif key == "class":
                d["classes"].append({"weight": float(v[2]), "tns": float(v[4])})
            elif key in ("freq", "mass", "root", "subin", "dwell"):
                d["classes"][int(v[0])][key] = np.array([float(x) for x in v[1:]])
            elif key in ("cat_rates", "cat_props"):
                d[key] = np.array([float(x) for x in v])
            elif key == "tree":
                d["tree"] = rest.strip()
            elif key in ("logl", "logl_init", "alpha"):
                d[key] = float(v[0])
            elif key in ("nsite", "nstates", "nmix", "ncat"):
                d[key] = int(v[0])
            else:
                d[key] = rest
    return d if end else None


def echo_matches(d, profiles, weights, rates, tol=1e-10):
    """the dump's echoes must equal the model IQ-TREE was given (catches stale files)"""
    if d is None or len(d["classes"]) != len(profiles):
        return False
    w = np.asarray(weights) / np.sum(weights)
    r = np.asarray(rates) / float(np.dot(w, rates))
    for k, cl in enumerate(d["classes"]):
        p = np.asarray(profiles[k]) / np.sum(profiles[k])
        if abs(cl["weight"] - w[k]) > tol * max(1.0, w[k]) or abs(cl["tns"] - r[k]) > tol * max(1.0, r[k]):
            return False
        if np.max(np.abs(cl["freq"] - p)) > tol:
            return False
    return True


# ---------------- profile matching ----------------

def sym_kl(p, q):
    p = np.asarray(p) + 1e-12
    q = np.asarray(q) + 1e-12
    return float(np.sum((p - q) * np.log(p / q)))


def match_profiles(A, B):
    """exact minimum-cost assignment of profiles A[i] -> B[j] by symmetric KL (bitmask DP).
    Returns (pairs [(i, j, cost)], mutual-nearest-neighbour flags)."""
    n, m = len(A), len(B)
    assert n <= m <= 16, "matching supports up to 16 classes"
    C = np.array([[sym_kl(a, b) for b in B] for a in A])
    INF = float("inf")
    best = {0: (0.0, [])}
    for i in range(n):
        nxt = {}
        for mask, (cost, assign) in best.items():
            for j in range(m):
                if mask & (1 << j):
                    continue
                nm, nc = mask | (1 << j), cost + C[i, j]
                if nc < nxt.get(nm, (INF, None))[0]:
                    nxt[nm] = (nc, assign + [j])
        best = nxt
    cost, assign = min(best.values(), key=lambda x: x[0])
    pairs = [(i, assign[i], C[i, assign[i]]) for i in range(n)]
    mnn = [int(np.argmin(C[i])) == j and int(np.argmin(C[:, j])) == i for i, j, _ in pairs]
    return pairs, mnn


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_list(path):
    with open(path) as f:
        return [l.split()[0] for l in f if l.strip() and not l.startswith("#")]
