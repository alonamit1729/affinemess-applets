"""Check certificates of non-realisability of distance patterns.  Pure Python, integers only.

A pattern on n points is a list of colours, one for each pair (0,1), (0,2), ..., (n-2,n-1); the colours
are 0..k-1.  Suppose points p_0..p_{n-1} of the plane have this pattern: pairs of the same colour are at
the same distance, and x_c > 0 is the squared distance of colour c, with x_a != x_b for a != b.

For a base point b the numbers G_ij = x_c(b,i) + x_c(b,j) - x_c(i,j) (with G_ii = 2 x_c(b,i)) are twice
the inner products of the vectors p_i - p_b and p_j - p_b.  These vectors lie in a plane, so the matrix
G has rank at most 2 and each of its 3x3 minors is zero.

A certificate consists of a non-zero integer `scale`, a list `product` of factors x_c or (x_a - x_b)
with exponents, and a list `terms` of minors (b, rows, columns), each with a polynomial multiplier.
check() verifies the polynomial identity

        scale * product  =  sum of multiplier * minor.

At the point x the right side is 0 and the left side is not.  So no such points exist.
"""
from itertools import combinations


def multiply(f, g):
    out = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}


def add(f, g, sign=1):
    out = dict(f)
    for m, c in g.items():
        out[m] = out.get(m, 0) + sign * c
    return {m: c for m, c in out.items() if c}


def check(cert):
    n, pattern = cert["n"], cert["pattern"]
    pairs = list(combinations(range(n), 2))
    assert len(pattern) == len(pairs) and all(isinstance(c, int) for c in pattern)
    k = max(pattern) + 1
    assert set(pattern) == set(range(k))
    colour = {}
    for (i, j), c in zip(pairs, pattern):
        colour[i, j] = colour[j, i] = c
    one = {(0,) * k: 1}
    var = lambda c: {tuple(int(i == c) for i in range(k)): 1}

    def entry(b, i, j):                       # G_ij for the base point b, a linear form
        if i == j:
            return add(var(colour[b, i]), var(colour[b, i]))
        return add(add(var(colour[b, i]), var(colour[b, j])), var(colour[i, j]), -1)

    def minor(b, rows, cols):                 # 3x3 determinant, expanded along the first row
        m = [[entry(b, i, j) for j in cols] for i in rows]
        det2 = lambda r1, r2, c1, c2: add(multiply(m[r1][c1], m[r2][c2]), multiply(m[r1][c2], m[r2][c1]), -1)
        out = multiply(m[0][0], det2(1, 2, 1, 2))
        out = add(out, multiply(m[0][1], det2(1, 2, 0, 2)), -1)
        return add(out, multiply(m[0][2], det2(1, 2, 0, 1)))

    left = {(0,) * k: cert["scale"]}
    assert isinstance(cert["scale"], int) and cert["scale"] != 0
    for factor, exponent in cert["product"]:
        assert isinstance(exponent, int) and exponent >= 1 and all(isinstance(c, int) and 0 <= c < k for c in factor)
        if len(factor) == 1:
            form = var(factor[0])
        else:
            a, b = factor
            assert a != b
            form = add(var(a), var(b), -1)
        for _ in range(exponent):
            left = multiply(left, form)

    right = {}
    for b, rows, cols, multiplier in cert["terms"]:
        for triple in (rows, cols):
            assert len(set(triple)) == 3 and b not in triple and all(isinstance(i, int) and 0 <= i < n for i in triple)
        assert 0 <= b < n
        a = {}
        for coefficient, exponents in multiplier:
            assert isinstance(coefficient, int) and len(exponents) == k and all(isinstance(e, int) and e >= 0 for e in exponents)
            a = add(a, {tuple(exponents): coefficient})
        right = add(right, multiply(a, minor(b, rows, cols)))
    assert left == right and left, "the identity does not hold"
    return True


def read(name):
    """Certificates from a file with one JSON object per line (.jsonl or .jsonl.gz)."""
    import gzip, json
    opener = gzip.open if name.endswith(".gz") else open
    with opener(name, "rt") as fh:
        return [json.loads(line) for line in fh if line.strip()]


if __name__ == "__main__":
    import sys
    total = 0
    for name in sys.argv[1:]:
        certs = read(name)
        for c in certs:
            check(c)
        total += len(certs)
        print(f"{name}: {len(certs)} certificates, every identity holds")
    print(total, "certificates checked")
