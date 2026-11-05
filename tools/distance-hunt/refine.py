import numpy as np
def refine(coords, cls):
    z = np.array([complex(float(x), float(y)) for x, y in coords])
    n = len(z); w = (z - z[0]) / (z[1] - z[0])
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    c01 = cls[0]; free = sorted(set(cls) - {c01})
    v = np.concatenate([np.ravel(np.column_stack([w[2:].real, w[2:].imag])),
                        [np.mean([abs(w[i]-w[j])**2 for (i, j), c in zip(pairs, cls) if c == f]) for f in free]])
    m = len(v)
    for it in range(60):
        P = np.vstack([[0, 0], [1, 0], v[:2*(n-2)].reshape(-1, 2)])
        T = {c01: 1.0}; T.update({f: v[2*(n-2)+k] for k, f in enumerate(free)})
        F = np.zeros(len(pairs)); J = np.zeros((len(pairs), m))
        for r, ((i, j), c) in enumerate(zip(pairs, cls)):
            d = P[i] - P[j]; F[r] = d @ d - T[c]
            if i >= 2: J[r, 2*(i-2):2*(i-2)+2] += 2*d
            if j >= 2: J[r, 2*(j-2):2*(j-2)+2] -= 2*d
            if c != c01: J[r, 2*(n-2)+free.index(c)] -= 1
        step = np.linalg.lstsq(J, F, rcond=None)[0]
        v = v - step
        if np.max(np.abs(step)) < 1e-15: break
    P = np.vstack([[0, 0], [1, 0], v[:2*(n-2)].reshape(-1, 2)])
    zz = P[:, 0] + 1j * P[:, 1]
    return list(zz * (z[1] - z[0]) + z[0]), float(np.max(np.abs(F)))
