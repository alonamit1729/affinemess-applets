"""A check on the 'kept' side.  It is not part of any proof; it tests that the search misses nothing.

For a candidate without a certificate, look (modulo the prime P of find.py) for a point at which all
minors vanish and no forbidden factor does.  If there is one, no certificate can exist, unless P divides
its scale.

How: the functionals that vanish on the ideal in degrees 5 and 6 describe the solution set.  If it is
finite, multiplication by a linear form acts on a finite-dimensional space, and the product q of all
forbidden factors vanishes at every solution exactly when its multiplication operator is nilpotent.
If the solution set contains curves inside forbidden hyperplanes, divide the ideal by one or two
forbidden factors first (pull the functionals back from degree 7 or 8).

    python3 kept_check.py          (about five minutes for the 50 kept candidates of the pilot)
"""
import itertools, json, random, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from find import *


def not_nilpotent(k, N5, N6, rng):
    """N5, N6: functionals vanishing on an ideal in degrees 5, 6, of the same rank h (finitely many
    points).  True if the product of all forbidden factors does not vanish at all of them."""
    h = N5.shape[0]
    R5, _ = rref(N5)
    piv = pivots(R5)

    def coords(form):
        Pm = np.zeros((N6.shape[0], R5.shape[1]), dtype=np.int64)
        for a, c in enumerate(form):
            if c:
                Pm = (Pm + c * N6[:, mul_index(k, 5, a)]) % P
        X = Pm[:, piv]
        if not ((X.astype(object) @ R5.astype(object)) % P == Pm).all():
            return None
        return to_flint(X)

    XL = coords([rng.randrange(1, P) for _ in range(k)])
    if XL is None or XL.rank() < h:
        return None
    XLinv = XL.inv()
    T = None
    for f in factors_of(k):
        form = [0] * k
        if f[0] == "x":
            form[f[1]] = 1
        else:
            form[f[1]], form[f[2]] = 1, P - 1
        X = coords(form)
        if X is None:
            return None
        T = XLinv * X if T is None else T * (XLinv * X)
    for _ in range(max(1, h.bit_length())):
        T = T * T
    return T.rank() > 0


def witness2(n, col, rng):
    k = max(col) + 1
    gens = [f for _, f in minors(n, col)]
    B = ideal_basis(k, gens, 8)
    N = {d: nullspace(B[d], len(monomials(k, d))) for d in (5, 6, 7, 8)}
    fac = factors_of(k)
    tried = 0
    for depth in (0, 1, 2):
        for combo in itertools.combinations_with_replacement(range(len(fac)), depth):
            M5, M6 = N[5 + depth], N[6 + depth]
            for step, j in enumerate(combo):
                M5 = rref(pull(k, 5 + depth - step, M5, fac[j]))[0]
                M6 = rref(pull(k, 6 + depth - step, M6, fac[j]))[0]
            if M5.shape[0] == 0 or M5.shape[0] != M6.shape[0]:
                continue
            tried += 1
            if not_nilpotent(k, M5, M6, rng):
                return ("witness", depth, [fac[j] for j in combo], M5.shape[0])
    return ("none", tried)


if __name__ == "__main__":
    import os
    from collections import Counter
    rng = random.Random(11)
    here = os.path.dirname(os.path.abspath(__file__))
    kept = json.load(open(os.path.join(here, "pilot6_kept.json")))["patterns"]
    results = [witness2(6, tuple(p), rng) for p in kept]
    print(len(kept), "kept candidates;  (outcome, number of factors divided out) -> candidates:")
    print(dict(Counter((w[0], w[1]) if w[0] == "witness" else (w[0],) for w in results)))
