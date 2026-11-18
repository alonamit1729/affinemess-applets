"""The search of find.py, organised for speed (the certificates are of the same kind and are checked
by the same check.py).

Differences from find.py:
  * the minors are computed with numpy for all base points, rows and columns at once;
  * for products of degree 3, 4 and 5 ALL products of forbidden factors are tested against the ideal at
    once (one matrix product), so the search in these degrees is exhaustive, not a depth-first search;
  * among the products found, a few are tried and the smallest certificate is kept.

Needs numpy and python-flint.
"""
import itertools
import json
import math
from functools import lru_cache

import numpy as np
from flint import fmpz_mat

from find import P, factors_of, from_flint, mon_index, monomials, mul_index, nullspace, pivots, rref, to_flint

PERMS = (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1))


# ------------------------------------------------------------------ minors
@lru_cache(None)
def minor_labels(n):
    """(base point, rows, columns) of all minors, principal ones first; and the same as index arrays."""
    labels = []
    for b in range(n):
        others = [i for i in range(n) if i != b]
        triples = list(itertools.combinations(others, 3))
        labels += [(0, b, R, R) for R in triples]
        labels += [(1, b, R, C) for R in triples for C in triples if R < C]
    labels.sort()
    return ([(b, R, C) for _, b, R, C in labels], np.array([l[1] for l in labels]),
            np.array([l[2] for l in labels]), np.array([l[3] for l in labels]))


@lru_cache(None)
def cubic_table(k):
    """For each (a, b, c): the index of the monomial x_a x_b x_c among the cubics."""
    idx = mon_index(k, 3)
    out = np.zeros((k, k, k), dtype=np.int64)
    for a, b, c in itertools.product(range(k), repeat=3):
        e = [0] * k
        e[a] += 1
        e[b] += 1
        e[c] += 1
        out[a, b, c] = idx[tuple(e)]
    return out.ravel()


def minors(n, col):
    """All distinct non-zero minors of the pattern: (labels, integer matrix with one cubic per row).
    Same list, in the same order, as find.minors."""
    k = max(col) + 1
    labels, base, R, C = minor_labels(n)
    c = np.zeros((n, n), dtype=np.int64)
    c[np.triu_indices(n, 1)] = col
    c = c + c.T
    A = np.eye(k, dtype=np.int64)[c]
    A[np.arange(n), np.arange(n)] = 0
    G = A[:, :, None, :] + A[:, None, :, :] - A[None, :, :, :]          # G[b, i, j] as a linear form
    U = G[base[:, None, None], R[:, :, None], C[:, None, :]]            # the 3x3 matrices of linear forms
    cube = np.zeros((len(labels), k ** 3), dtype=np.int64)
    for perm, sign in PERMS:
        cube += sign * np.einsum("ma,mb,mc->mabc", U[:, 0, perm[0]], U[:, 1, perm[1]], U[:, 2, perm[2]]).reshape(len(labels), -1)
    F = np.zeros((len(monomials(k, 3)), len(labels)), dtype=np.int64)
    np.add.at(F, cubic_table(k), cube.T)
    F = F.T
    first = (F != 0).argmax(axis=1)
    lead = F[np.arange(len(F)), first]
    _, index = np.unique(F * np.sign(lead)[:, None], axis=0, return_index=True)    # equal up to sign
    keep = sorted(i for i in index.tolist() if lead[i] != 0)
    return [labels[i] for i in keep], F[keep]


def independent_minors(F):
    """(indices of the first maximal independent set of rows of F modulo P, functionals vanishing on
    all rows).  Usually the first rows already span everything, so only a small matrix is reduced."""
    block = min(len(F), 96)
    S = pivots(rref(F[:block].T)[0])
    N = nullspace(F[S], F.shape[1])
    if block < len(F) and N.shape[0] and (np.mod(F[block:].astype(np.float64) @ N.T.astype(np.float64), P) != 0).any():
        S = pivots(rref(F.T)[0])
        N = nullspace(F[S], F.shape[1])
    return np.array(S), N


# ------------------------------------------------------------------ all products of d forbidden factors
@lru_cache(None)
def forms(k):
    fac = factors_of(k)
    out = np.zeros((len(fac), k), dtype=np.int64)
    for j, f in enumerate(fac):
        if f[0] == "x":
            out[j, f[1]] = 1
        else:
            out[j, f[1]], out[j, f[2]] = 1, -1
    return out


@lru_cache(None)
def product_table(k, d):
    """(list of multisets of d factor indices, matrix of their products as forms of degree d)."""
    L = forms(k)
    if d == 1:
        return [(j,) for j in range(len(L))], L.copy()
    prev, Q = product_table(k, d - 1)
    where = {c: i for i, c in enumerate(prev)}
    combos = list(itertools.combinations_with_replacement(range(len(L)), d))
    out = np.zeros((len(combos), len(monomials(k, d))), dtype=np.int64)
    rows = np.array([where[c[:-1]] for c in combos])
    last = np.array([c[-1] for c in combos])
    for a in range(k):
        out[:, mul_index(k, d - 1, a)] += Q[rows] * L[last, a][:, None]
    return combos, out


@lru_cache(None)
def distinct_factors(k, d):
    return np.array([len(set(c)) for c in product_table(k, d)[0]])


@lru_cache(None)
def product_table_float(k, d):
    return product_table(k, d)[1].astype(np.float64)


def products_in_ideal(k, d, N):
    """Indices (in product_table) of the products killed by all functionals in N."""
    Q = product_table_float(k, d)
    if N.shape[0] == 0:
        return np.arange(len(Q))
    alive = np.arange(len(Q))
    for start in range(0, N.shape[0], 4):                                # a few functionals at a time
        Z = Q[alive] @ N[start:start + 4].T.astype(np.float64)           # exact: all values below 2^53
        alive = alive[(np.mod(Z, P) == 0).all(axis=1)]
        if len(alive) == 0:
            break
    return alive


# ------------------------------------------------------------------ exact multipliers
@lru_cache(None)
def shift_table(k, d):
    """For each monomial m of degree d - 3: the positions of m * (cubic monomials) in degree d."""
    idx = mon_index(k, d)
    return [np.array([idx[tuple(a + b for a, b in zip(mm, m))] for mm in monomials(k, 3)]) for m in monomials(k, d - 3)]


def solve(k, d, g, F, use):
    """Integers `scale` and multipliers with  scale * g = sum multiplier * minor  (g: integer vector of
    a form of degree d; `use`: an independent set of rows of F).
    Returns (scale, {minor index: [[coefficient, exponents], ...]}) or None."""
    mult = monomials(k, d - 3)
    if d == 3:
        A = F[use]                                                      # the solution is unique
        S = list(range(len(use)))
    else:
        C = np.zeros((len(use) * len(mult), len(monomials(k, d))), dtype=np.int64)
        for j, positions in enumerate(shift_table(k, d)):
            C[j::len(mult), positions] = F[use]
        R, r = rref(np.hstack([C.T, (g % P).reshape(-1, 1)]))
        piv = pivots(R)
        if piv and piv[-1] == C.shape[0]:
            return None
        S = [pc for row, pc in zip(R, piv) if row[-1] != 0]
        if not S:
            return None
        A = C[S]                                                        # len(S) x monomials, exact integers
    rows = pivots(rref(A)[0])                                           # independent rows of the transpose
    if len(rows) != len(S):
        return None
    sq = fmpz_mat(len(S), len(S), A[:, rows].T.ravel().tolist())
    try:
        X = sq.solve(fmpz_mat(len(S), 1, [int(v) for v in g[rows]]))
    except ZeroDivisionError:
        return None
    sol = [X[i, 0] for i in range(len(S))]
    den = 1
    for v in sol:
        den = den * int(v.q) // math.gcd(den, int(v.q))
    ints = [int(v.p) * (den // int(v.q)) for v in sol]
    common = math.gcd(den, *ints)
    den, ints = den // common, [t // common for t in ints]
    # exact check on all coefficients (Python integers)
    total = [0] * A.shape[1]
    for t, row in zip(ints, A.tolist()):
        if t:
            for j, v in enumerate(row):
                if v:
                    total[j] += t * v
    if total != [den * int(v) for v in g]:
        return None
    terms = {}
    for s, t in zip(S, ints):
        if t:
            terms.setdefault(int(use[s // len(mult)]), []).append([int(t), list(mult[s % len(mult)])])
    return den, terms


def certificate(n, col, d, combo, scale, terms, labels):
    fac = factors_of(max(col) + 1)
    exps = {}
    for j in combo:
        exps[j] = exps.get(j, 0) + 1
    return {"n": n, "pattern": [int(c) for c in col], "scale": scale,
            "product": [[list(fac[j][1:]), e] for j, e in sorted(exps.items())],
            "terms": [[labels[i][0], list(labels[i][1]), list(labels[i][2]), a] for i, a in sorted(terms.items())]}


def size_of(cert):
    return len(json.dumps(cert, separators=(",", ":")))


def search(n, col, dmax=5, tries=1):
    """Certificate with a product of the smallest possible degree d <= dmax (d <= 5), or None.
    Returns (certificate or None, degree or None)."""
    col = tuple(int(c) for c in col)
    k = max(col) + 1
    labels, F = minors(n, col)
    use, N = independent_minors(F)
    B = None
    for d in range(3, min(dmax, 5) + 1):
        if d > 3:
            if B is None:
                B = rref(F[use])[0]
            stack = np.zeros((k * B.shape[0], len(monomials(k, d))), dtype=np.int64)
            for a in range(k):
                stack[a * B.shape[0]:(a + 1) * B.shape[0], mul_index(k, d - 1, a)] = B
            B = rref(stack)[0]
            N = nullspace(B, len(monomials(k, d)))
        hits = products_in_ideal(k, d, N)
        if len(hits) == 0:
            continue
        combos, Q = product_table(k, d)
        # prefer products with few distinct factors; with tries > 1, keep the smallest certificate
        order = hits[np.lexsort((hits, distinct_factors(k, d)[hits]))][:tries].tolist()
        best = None
        for i in order:
            got = solve(k, d, Q[i], F, use)
            if got is None:
                continue
            cert = certificate(n, col, d, combos[i], got[0], got[1], labels)
            if best is None or size_of(cert) < size_of(best):
                best = cert
        if best is not None:
            return best, d
    return None, None
