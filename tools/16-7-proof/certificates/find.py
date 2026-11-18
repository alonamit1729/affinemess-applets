"""Search for certificates of non-realisability of distance patterns (pilot).

A pattern on n points with k colours assigns to each pair {i, j} a colour c(i, j) in 0..k-1.  If points
p_0..p_{n-1} of the plane have this pattern, with squared distances x_0..x_{k-1} (distinct, positive),
then for every base point b the matrix

        G^b_{ij} = x_{c(b,i)} + x_{c(b,j)} - x_{c(i,j)}        (i, j != b;  G^b_{ii} = 2 x_{c(b,i)})

is twice the Gram matrix of the vectors p_i - p_b, so it has rank <= 2 and all its 3x3 minors vanish.
The minors are homogeneous cubics f_1, f_2, ... in x.  A certificate is an identity

        d * prod_c x_c^{e_c} * prod_{a<b} (x_a - x_b)^{e_ab}  =  sum_i a_i f_i        (d a non-zero integer)

with polynomials a_i: the right side vanishes at x, the left side does not.

This module finds such identities: linear algebra modulo a prime to find the product, exact linear
algebra over Q for the a_i.  Nothing here has to be trusted: check.py verifies each identity.

Needs numpy and python-flint.
"""
import itertools
import math
from functools import lru_cache

import numpy as np
from flint import fmpz_mat, nmod_mat

P = 1073741827                      # search prime (products of two residues fit in int64)


# ------------------------------------------------------------------ monomials
@lru_cache(None)
def monomials(k, d):
    """Exponent tuples of the monomials of degree d in k variables, in a fixed order."""
    if k == 1:
        return [(d,)]
    return [(e,) + rest for e in range(d, -1, -1) for rest in monomials(k - 1, d - e)]


@lru_cache(None)
def mon_index(k, d):
    return {m: i for i, m in enumerate(monomials(k, d))}


@lru_cache(None)
def mul_index(k, d, a):
    """For each monomial m of degree d, the index of m * x_a among the monomials of degree d + 1."""
    idx = mon_index(k, d + 1)
    return np.array([idx[m[:a] + (m[a] + 1,) + m[a + 1:]] for m in monomials(k, d)], dtype=np.int64)


def factors_of(k):
    """The forbidden linear forms: ('x', c) is x_c, ('d', a, b) is x_a - x_b."""
    return [("x", c) for c in range(k)] + [("d", a, b) for a, b in itertools.combinations(range(k), 2)]


# ------------------------------------------------------------------ polynomials as dicts
def pmul(f, g):
    out = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            v = out.get(m, 0) + c1 * c2
            if v:
                out[m] = v
            else:
                out.pop(m, None)
    return out


def padd(f, g, s=1):
    out = dict(f)
    for m, c in g.items():
        v = out.get(m, 0) + s * c
        if v:
            out[m] = v
        else:
            out.pop(m, None)
    return out


def linear(k, coeffs):
    """{variable: coefficient} -> polynomial dict."""
    out = {}
    for v, c in coeffs.items():
        if c:
            out[tuple(1 if i == v else 0 for i in range(k))] = c
    return out


def factor_poly(k, f):
    return linear(k, {f[1]: 1}) if f[0] == "x" else linear(k, {f[1]: 1, f[2]: -1})


def product_poly(k, exps):
    """exps: list of (factor, exponent)."""
    out = {tuple([0] * k): 1}
    for f, e in exps:
        for _ in range(e):
            out = pmul(out, factor_poly(k, f))
    return out


# ------------------------------------------------------------------ the minors
def colour_table(n, col):
    c = [[None] * n for _ in range(n)]
    for e, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        c[i][j] = c[j][i] = col[e]
    return c


def minor(n, col, b, R, C, k=None, c=None):
    """The 3x3 minor of G^b on rows R and columns C, as a polynomial dict."""
    k = k or max(col) + 1
    c = c or colour_table(n, col)

    def entry(i, j):
        if i == j:
            return linear(k, {c[b][i]: 2})
        co = {}
        for v, s in ((c[b][i], 1), (c[b][j], 1), (c[i][j], -1)):
            co[v] = co.get(v, 0) + s
        return linear(k, co)

    M = [[entry(i, j) for j in C] for i in R]
    out = {}
    for perm, sign in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
        out = padd(out, pmul(pmul(M[0][perm[0]], M[1][perm[1]]), M[2][perm[2]]), sign)
    return out


def minors(n, col):
    """All distinct non-zero minors, principal ones (R = C) first: list of ((b, R, C), poly)."""
    k = max(col) + 1
    c = colour_table(n, col)
    out, seen = [], set()
    labels = []
    for b in range(n):
        others = [i for i in range(n) if i != b]
        triples = list(itertools.combinations(others, 3))
        labels += [(0, b, R, R) for R in triples]
        labels += [(1, b, R, C) for R in triples for C in triples if R < C]
    for _, b, R, C in sorted(labels):
        f = minor(n, col, b, R, C, k, c)
        if not f:
            continue
        key = frozenset(f.items())
        neg = frozenset((m, -v) for m, v in f.items())
        if key in seen or neg in seen:
            continue
        seen.add(key)
        out.append(((b, R, C), f))
    return out


# ------------------------------------------------------------------ linear algebra modulo P
def to_flint(A):
    A = np.asarray(A, dtype=np.int64) % P
    return nmod_mat(A.shape[0], A.shape[1], A.ravel().tolist(), P)


def from_flint(M, rows=None):
    r = M.nrows() if rows is None else rows
    c = M.ncols()
    flat = M.entries()
    return np.array([int(v) for v in flat[: r * c]], dtype=np.int64).reshape(r, c)


def rref(A):
    """Reduced row echelon basis of the row space (numpy int64), and its rank."""
    if A.shape[0] == 0 or A.shape[1] == 0:
        return np.zeros((0, A.shape[1]), dtype=np.int64), 0
    R, rank = to_flint(A).rref()
    return from_flint(R, rank), rank


def rank(A):
    if A.shape[0] == 0 or A.shape[1] == 0:
        return 0
    return to_flint(A).rank()


def nullspace(A, ncols):
    """Rows spanning {v : A v = 0}."""
    if A.shape[0] == 0:
        return np.eye(ncols, dtype=np.int64)
    X, nullity = to_flint(A).nullspace()
    full = from_flint(X)
    return full[:, :nullity].T.copy()


def vec(k, d, f):
    idx = mon_index(k, d)
    v = np.zeros(len(idx), dtype=np.int64)
    for m, c in f.items():
        v[idx[m]] = c % P
    return v


def shift(k, d, B, a):
    """Rows of B are forms of degree d; returns the rows multiplied by x_a (degree d + 1)."""
    out = np.zeros((B.shape[0], len(monomials(k, d + 1))), dtype=np.int64)
    out[:, mul_index(k, d, a)] = B
    return out


def ideal_basis(k, gens, D):
    """Echelon bases (mod P) of the degree-d parts of the ideal for d = 3..D: {d: array}."""
    B = {3: rref(np.array([vec(k, 3, f) for f in gens]))[0]}
    for d in range(4, D + 1):
        B[d] = rref(np.vstack([shift(k, d - 1, B[d - 1], a) for a in range(k)]))[0]
    return B


def pull(k, d, N, f):
    """N: functionals on forms of degree d (rows).  Returns the functionals h -> N(f * h) on degree d - 1."""
    if f[0] == "x":
        return N[:, mul_index(k, d - 1, f[1])]
    return (N[:, mul_index(k, d - 1, f[1])] - N[:, mul_index(k, d - 1, f[2])]) % P


def find_product(k, D, N, budget=300):
    """A product of D forbidden factors killed by all functionals in N (the functionals on forms of
    degree D that vanish on the ideal): such a product lies in the ideal.  Depth-first search, trying
    the factor that leaves the fewest functionals first.  Returns a list of factors, or None."""
    fac = factors_of(k)
    seen, nodes = set(), [0]

    def search(d, N, chosen):
        if N.shape[0] == 0:
            return chosen + [0] * d               # every form of degree d is in the ideal: pad with x_0
        if d == 0:
            return None
        key = tuple(sorted(chosen))
        if key in seen:
            return None
        seen.add(key)
        nodes[0] += 1
        if nodes[0] > budget:
            return None
        options = []
        for j in range(len(fac)):
            M = pull(k, d, N, fac[j])
            options.append((rank(M), j, M))
        options.sort(key=lambda t: (t[0], t[1]))
        for r, j, M in options:
            if d == 1 and r > 0:
                break
            got = search(d - 1, rref(M)[0], chosen + [j])
            if got is not None:
                return got
        return None

    got = search(D, N, [])
    return None if got is None else [fac[j] for j in got]


# ------------------------------------------------------------------ exact cofactors
def pivots(R):
    """Pivot column of each row of a reduced echelon matrix (numpy)."""
    return [int(np.flatnonzero(row)[0]) for row in R]


def solve_identity(n, col, D, prod, gens, tries=1, rng=None):
    """Given a product of D forbidden factors that lies in the ideal modulo P, find integers d != 0
    and polynomials a_i with  d * product = sum a_i * gens[i]  exactly.  Returns a certificate dict
    (see check.py for the format) or None if the exact identity does not exist on the chosen support."""
    from flint import fmpq_mat
    k = max(col) + 1
    g = product_poly(k, [(f, 1) for f in prod])
    idx = mon_index(k, D)
    gvec = vec(k, D, g)
    mult = monomials(k, D - 3)
    use = range(len(gens))
    if D >= 5:                                   # many multiples: an independent set of minors is enough
        R3, _ = rref(np.array([vec(k, 3, f) for _, f in gens]).T)
        use = pivots(R3)
    cols = [(i, m) for i in use for m in mult]                         # the products m * f_i
    polys = {}

    def column(i, m):
        if (i, m) not in polys:
            polys[i, m] = {tuple(a + b for a, b in zip(mm, m)): c for mm, c in gens[i][1].items()}
        return polys[i, m]

    C = np.zeros((len(cols), len(idx)), dtype=np.int64)
    for j, (i, m) in enumerate(cols):
        for mm, c in column(i, m).items():
            C[j, idx[mm]] = c % P
    best = None
    order0 = list(range(len(cols)))
    for attempt in range(tries):
        order = order0 if attempt == 0 else list(rng.permutation(len(cols)))
        # earliest independent columns in this order: pivots of the echelon form of the transpose
        R, r = rref(np.hstack([C[order].T, gvec.reshape(-1, 1)]))
        piv = pivots(R)
        if piv and piv[-1] == len(order):
            return None                                                   # g is not in the span modulo P
        support = [order[pc] for row, pc in zip(R, piv) if row[-1] != 0]
        if best is None or len(support) < len(best):
            best = support
    S = sorted(best)
    # exact solution on the support
    A = fmpq_mat(len(idx), len(S) + 1, [0] * (len(idx) * (len(S) + 1)))
    for j, s in enumerate(S):
        for mm, c in column(*cols[s]).items():
            A[idx[mm], j] = c
    for mm, c in g.items():
        A[idx[mm], len(S)] = c
    R, r = A.rref()
    if r != len(S):
        return None
    sol = [R[j, len(S)] for j in range(len(S))]
    if any(R[j, j] != 1 for j in range(len(S))):
        return None
    den = 1
    for v in sol:
        den = den * int(v.q) // math.gcd(den, int(v.q))
    ints = [int(v.p) * (den // int(v.q)) for v in sol]
    common = math.gcd(den, *ints)
    den, ints = den // common, [t // common for t in ints]
    terms = {}
    for s, t in zip(S, ints):
        if t:
            i, m = cols[s]
            terms.setdefault(i, []).append([t, list(m)])
    exps = {}
    for f in prod:
        exps[f] = exps.get(f, 0) + 1
    cert = {"n": n, "pattern": list(col), "scale": den,
            "product": [[list(f[1:]), e] for f, e in sorted(exps.items(), key=lambda t: (len(t[0]), t[0]))],
            "terms": [[gens[i][0][0], list(gens[i][0][1]), list(gens[i][0][2]), a] for i, a in sorted(terms.items())]}
    return cert


def search(n, col, dmax=6, tries=1, rng=None, budget=300):
    """Look for a certificate with a product of degree 3, 4, ..., dmax.  Returns (certificate or None,
    the number of independent functionals vanishing on the ideal in each degree tried, the degree)."""
    k = max(col) + 1
    gens = minors(n, col)
    B = {3: rref(np.array([vec(k, 3, f) for _, f in gens]))[0]}
    hf = []
    for D in range(3, dmax + 1):
        if D > 3:
            B[D] = rref(np.vstack([shift(k, D - 1, B[D - 1], a) for a in range(k)]))[0]
        N = nullspace(B[D], len(monomials(k, D)))
        hf.append(N.shape[0])
        prod = find_product(k, D, N, budget)
        if prod is not None:
            cert = solve_identity(n, col, D, prod, gens, tries, rng)
            return cert, hf, D
    return None, hf, None


def size_of(cert):
    import json
    return len(json.dumps(cert, separators=(",", ":")))
