"""Symmetry group of a finite planar point set: all permutations p of the points with
|p(i)p(j)| = |ij| for all i,j (relative tolerance), found by backtracking; then classified as
trivial / mirror / C_k / D_k by counting orientation-reversing elements."""
import math, itertools
import numpy as np

def dmat(P):
    P = np.asarray(P, float)
    return np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)

def automorphisms(P, rtol=1e-6):
    P = np.asarray(P, float); n = len(P); D = dmat(P); tol = rtol * max(D.max(), 1e-300)
    perms = []
    img = [-1] * n
    def rec(i):
        if i == n:
            perms.append(tuple(img)); return
        used = set(img[:i])
        for c in range(n):
            if c in used: continue
            ok = True
            for j in range(i):
                if abs(D[img[j], c] - D[j, i]) > tol: ok = False; break
            if ok:
                img[i] = c; rec(i + 1); img[i] = -1
    rec(0)
    return perms

def orientation(P, perm, rtol=1e-6):
    """+1 if the isometry realising perm is a rotation/translation, -1 if a reflection, 0 if undecidable (collinear)."""
    P = np.asarray(P, float); n = len(P); D = dmat(P); sc = max(D.max(), 1e-300)
    for i, j, k in itertools.combinations(range(n), 3):
        a = (P[j] - P[i]); b = (P[k] - P[i]); cr = a[0] * b[1] - a[1] * b[0]
        if abs(cr) > 1e-6 * sc * sc:
            a2 = P[perm[j]] - P[perm[i]]; b2 = P[perm[k]] - P[perm[i]]; cr2 = a2[0] * b2[1] - a2[1] * b2[0]
            return 1 if cr * cr2 > 0 else -1
    return 0

def classify(P, rtol=1e-6):
    """returns dict(order, type, perms). type in 'trivial','mirror','C_k','D_k','collinear(order m)'."""
    perms = automorphisms(P, rtol)
    m = len(perms)
    if m == 1: return {'order': 1, 'type': 'trivial', 'perms': perms}
    ors = [orientation(P, p, rtol) for p in perms]
    if any(o == 0 for o in ors):
        return {'order': m, 'type': f'collinear(order {m})', 'perms': perms}
    r = sum(1 for o in ors if o < 0)
    if r == 0: t = f'C_{m}'
    else:
        k = m // 2
        t = 'mirror' if k == 1 else f'D_{k}'
    return {'order': m, 'type': t, 'perms': perms}

def e_dim(n, d): return 2 * n - 4 - n * (n - 1) // 2 + d

if __name__ == '__main__':
    # sanity checks
    import cmath
    sq = [[0, 0], [1, 0], [1, 1], [0, 1]]; print('square', classify(sq)['type'])
    pent = [[math.cos(2 * math.pi * k / 5), math.sin(2 * math.pi * k / 5)] for k in range(5)]; print('pentagon', classify(pent)['type'])
    tri = [[0, 0], [1, 0], [0.3, 0.7]]; print('generic tri', classify(tri)['type'])
    iso = [[0, 0], [1, 0], [0.5, 0.7]]; print('isoceles', classify(iso)['type'])
    rot = [[1, 0], [-1, 0], [0.3, 0.6], [-0.3, -0.6], [0.2, 0.1], [-0.2, -0.1]]; print('central', classify(rot)['type'])
    c3 = []
    for r, th in ((1, 0.1), (0.5, 0.9)):
        for k in range(3): c3.append([r * math.cos(th + 2 * math.pi * k / 3), r * math.sin(th + 2 * math.pi * k / 3)])
    print('C3', classify(c3)['type'])
    print('collinear', classify([[0, 0], [1, 0], [3, 0]])['type'], classify([[0, 0], [1, 0], [2, 0]])['type'])
