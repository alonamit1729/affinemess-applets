"""Alternative realization route for heavy patterns: msolve real roots of the Gram-minor system (variables s1..s_{k-1}, t),
then PSD / rank-2 check and coordinate reconstruction at 100 digits."""
import os, sys, json, time
from fractions import Fraction
from msolve_util import *
from realize import system_gram, colour_names
import mpmath as mp
from extend import mp, Classes, pattern, normalise, TOL
import itertools

def gram_realize(n, col, precision=400, timeout=14000):
    V, P = system_gram(n, col)
    t0 = time.time()
    out = msolve_real_raw(V, P, precision=precision, timeout=timeout)
    sols = parse_real(out)
    return V, sols, time.time() - t0

def sets_from_gram(n, col, V, sols):
    """for each real s-solution build the squared-distance matrix, check PSD rank 2, return point sets + diagnostics."""
    names = colour_names(col)
    E = list(itertools.combinations(range(n), 2))
    res = []
    for sol in sols:
        vv = {}
        for name, (a, b) in zip(V, sol):
            m = (Fraction(a) + Fraction(b)) / 2
            vv[name] = mp.mpf(m.numerator) / mp.mpf(m.denominator)
        val = {c: (mp.mpf(1) if names[c] is None else vv[names[c]]) for c in names}
        D = mp.zeros(n, n)
        for e, (i, j) in enumerate(E):
            D[i, j] = D[j, i] = val[col[e]]
        # Gram matrix w.r.t. point 0: G_ij = (D0i + D0j - Dij)/2, i,j = 1..n-1
        G = mp.zeros(n - 1, n - 1)
        for i in range(1, n):
            for j in range(1, n):
                G[i - 1, j - 1] = (D[0, i] + D[0, j] - D[i, j]) / 2
        ev = mp.eigsy(G)[0]
        ev = sorted([ev[i] for i in range(n - 1)])
        res.append({'vals': val, 'eig': ev, 'psd': ev[0] > -TOL and ev[-3] < TOL and ev[-2] > TOL, 'G': G})
    return res

def coords_from_gram(G, n):
    """points p0=0, p1..p_{n-1} from a PSD rank-2 Gram matrix (eigen-decomposition)."""
    ev, U = mp.eigsy(G)
    idx = sorted(range(n - 1), key=lambda i: -ev[i])[:2]
    pts = [(mp.mpf(0), mp.mpf(0))]
    for i in range(n - 1):
        pts.append(tuple(mp.sqrt(max(ev[k], mp.mpf(0))) * U[i, k] for k in idx))
    return pts
