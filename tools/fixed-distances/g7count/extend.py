"""Extension of a rigid core A (a finite planar set with |Delta(A)| <= 6 distances, given to 100 digits) to all
supersets X = A ∪ Q with |Q| = k_extra points and |Delta(X)| <= |Delta(A)| + budget, budget = 7 - |Delta(A)|.

All arithmetic in mpmath at 100 decimal digits on SQUARED distances; two values are "equal" iff they differ by
less than TOL = 1e-40 (sets are scaled so that the minimum squared distance is 1).  Soundness: the coordinates of
every core are algebraic numbers of small degree given to ~1e-95; a genuine equality is never missed, and a spurious
one (two distinct algebraic numbers within 1e-40) would only add a candidate set, which is then checked exactly.

Candidates for a point q of Q (every distance from q to A lies in Delta(A) ∪ N, |N| <= budget <= 2, |A| >= 10):
  * q has >= 2 distances to A in Delta(A): q is an intersection of two circles centred at points of A with radii
    in Delta(A);
  * otherwise >= |A| - 1 >= 9 of its distances lie in N, so one value occurs >= 5 times and q is the circumcentre of
    a non-collinear triple of A.
"""
import itertools, math
import mpmath as mp

mp.mp.dps = 100
TOL = mp.mpf('1e-40')


def d2(p, q):
    dx = p[0] - q[0]; dy = p[1] - q[1]
    return dx * dx + dy * dy


def normalise(P):
    """scale so that the minimum squared distance is 1 (translation/rotation irrelevant)."""
    m = min(d2(p, q) for p, q in itertools.combinations(P, 2))
    s = 1 / mp.sqrt(m)
    return [(p[0] * s, p[1] * s) for p in P]


class Values:
    """a growing list of distinct values (clustered with TOL), hashed on floor(x*10^35)."""
    SC = mp.mpf(10) ** 35

    def __init__(self, vals=()):
        self.v = []; self.grid = {}
        for x in vals: self.add(x)

    def find(self, x):
        k = int(mp.floor(x * self.SC))
        for kk in (k - 1, k, k + 1):
            for i in self.grid.get(kk, ()):
                if abs(x - self.v[i]) < TOL: return i
        return None

    def add(self, x):
        i = self.find(x)
        if i is None:
            self.v.append(x); i = len(self.v) - 1
            self.grid.setdefault(int(mp.floor(x * self.SC)), []).append(i)
        return i

    def __len__(self): return len(self.v)


def dist_values(P):
    V = Values()
    for p, q in itertools.combinations(P, 2): V.add(d2(p, q))
    return V


def circle_pts(p, r2p, q, r2q):
    dd = d2(p, q)
    if dd < TOL: return []
    d = mp.sqrt(dd)
    a = (r2p - r2q + dd) / (2 * d)
    h2 = r2p - a * a
    ux = (q[0] - p[0]) / d; uy = (q[1] - p[1]) / d
    mx = p[0] + a * ux; my = p[1] + a * uy
    if h2 < -TOL: return []
    if h2 < TOL: return [(mx, my)]
    h = mp.sqrt(h2)
    return [(mx - h * uy, my + h * ux), (mx + h * uy, my - h * ux)]


def circumcentre(a, b, c):
    D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    if abs(D) < TOL: return None
    A2 = a[0] ** 2 + a[1] ** 2; B2 = b[0] ** 2 + b[1] ** 2; C2 = c[0] ** 2 + c[1] ** 2
    ux = (A2 * (b[1] - c[1]) + B2 * (c[1] - a[1]) + C2 * (a[1] - b[1])) / D
    uy = (A2 * (c[0] - b[0]) + B2 * (a[0] - c[0]) + C2 * (b[0] - a[0])) / D
    return (ux, uy)


class PointSet:
    """dedupe of points with tolerance (grid hashing)."""
    def __init__(self):
        self.grid = {}; self.pts = []

    def _key(self, p):
        return (int(mp.floor(p[0] * mp.mpf(10) ** 30)), int(mp.floor(p[1] * mp.mpf(10) ** 30)))

    def add(self, p):
        kx, ky = self._key(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in self.grid.get((kx + dx, ky + dy), ()):
                    q = self.pts[i]
                    if abs(q[0] - p[0]) < TOL and abs(q[1] - p[1]) < TOL: return i
        self.pts.append(p); self.grid.setdefault((kx, ky), []).append(len(self.pts) - 1)
        return len(self.pts) - 1


def candidates(A, budget):
    """all candidate points q (not in A) with at most `budget` distance values to A outside Delta(A).
    Returns (vals, survivors) where vals is a Values object (first len(Delta(A)) entries = Delta(A)) and
    survivors = list of (point, newids) with newids = frozenset of value ids outside Delta(A)."""
    vals = dist_values(A); kA = len(vals)
    S = list(vals.v)
    ps = PointSet()
    for a in A: ps.add(a)
    nA = len(A)
    for i, j in itertools.combinations(range(nA), 2):
        for r1 in S:
            for r2 in S:
                for q in circle_pts(A[i], r1, A[j], r2): ps.add(q)
    for t in itertools.combinations(range(nA), 3):
        c = circumcentre(A[t[0]], A[t[1]], A[t[2]])
        if c is not None: ps.add(c)
    surv = []
    for q in ps.pts[nA:]:
        loc = Values(); ok = True
        for a in A:
            x = d2(q, a)
            if x < TOL: ok = False; break
            i = vals.find(x)
            if i is None or i >= kA:
                loc.add(x)
                if len(loc) > budget: ok = False; break
        if not ok: continue
        newids = frozenset(vals.add(x) for x in loc.v)   # register the new values globally
        surv.append((q, newids))
    return vals, kA, surv


def extensions(A, k_extra, budget):
    """all k_extra-subsets Q of candidates with |Delta(A ∪ Q)| <= |Delta(A)| + budget.
    Returns (list of point lists Q, stats dict)."""
    vals, kA, surv = candidates(A, budget)
    n = len(surv)
    # pairwise values
    pairv = [[None] * n for _ in range(n)]   # None = in Delta(A); else value id
    for i in range(n):
        for j in range(i + 1, n):
            x = d2(surv[i][0], surv[j][0])
            if x < TOL:
                pairv[i][j] = pairv[j][i] = -1      # coincident points (should not happen after dedupe)
                continue
            k = vals.add(x)
            v = None if k < kA else k
            pairv[i][j] = pairv[j][i] = v
    newids = [s[1] for s in surv]
    results = []
    cnt = [0]

    def dfs(chosen, N, start):
        if len(chosen) == k_extra:
            results.append(list(chosen)); return
        need = k_extra - len(chosen)
        for i in range(start, n - need + 1):
            cnt[0] += 1
            N2 = N | newids[i]
            if len(N2) > budget: continue
            ok = True
            for c in chosen:
                v = pairv[c][i]
                if v == -1: ok = False; break
                if v is not None:
                    N2 = N2 | {v}
                    if len(N2) > budget: ok = False; break
            if not ok: continue
            chosen.append(i); dfs(chosen, N2, i + 1); chosen.pop()
    dfs([], frozenset(), 0)
    Qs = [[surv[i][0] for i in r] for r in results]
    return Qs, {'ncand': n, 'nvals_total': len(vals), 'dfs_nodes': cnt[0], 'nsol': len(results)}


# ---------------------------------------------------------------- congruence / invariants
def inv_key(P, digits=20):
    """similarity invariant: sorted multiset of squared distances scaled by min, rounded."""
    P = normalise(P)
    ds = sorted(d2(p, q) for p, q in itertools.combinations(P, 2))
    return tuple(mp.nstr(x, digits) for x in ds)


def congruent(P, Q, tol=mp.mpf('1e-30')):
    """similarity test by vertex bijection (backtracking)."""
    n = len(P)
    if len(Q) != n: return False
    P = normalise(P); Q = normalise(Q)
    D1 = [[d2(a, b) for b in P] for a in P]; D2 = [[d2(a, b) for b in Q] for a in Q]
    r1 = [sorted(row) for row in D1]; r2 = [sorted(row) for row in D2]
    def rowmatch(x, y): return all(abs(u - v) < tol for u, v in zip(x, y))
    cand = [[j for j in range(n) if rowmatch(r1[i], r2[j])] for i in range(n)]
    if any(not c for c in cand): return False
    order = sorted(range(n), key=lambda i: len(cand[i]))
    perm = [-1] * n; used = [False] * n
    def bt(k):
        if k == n: return True
        i = order[k]
        for j in cand[i]:
            if used[j]: continue
            if all(abs(D1[i][order[t]] - D2[j][perm[order[t]]]) < tol for t in range(k)):
                perm[i] = j; used[j] = True
                if bt(k + 1): return True
                used[j] = False; perm[i] = -1
        return False
    return bt(0)


class Classes:
    """collection of point sets up to similarity."""
    def __init__(self):
        self.buckets = {}; self.reps = []

    def add(self, P, info=None):
        k = inv_key(P)
        for idx in self.buckets.get(k, []):
            if congruent(self.reps[idx][0], P):
                self.reps[idx][2] += 1
                return idx, False
        self.reps.append([P, info, 1]); idx = len(self.reps) - 1
        self.buckets.setdefault(k, []).append(idx)
        return idx, True


def pattern(P):
    """rgf colour tuple of the point set (edges in itertools.combinations order)."""
    V = Values(); col = []
    for p, q in itertools.combinations(P, 2): col.append(V.add(d2(p, q)))
    return tuple(col)


def to_floats(P):
    return [[float(p[0]), float(p[1])] for p in P]


def to_str(P, digits=60):
    return [[mp.nstr(p[0], digits), mp.nstr(p[1], digits)] for p in P]
