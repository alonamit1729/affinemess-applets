"""Check certificates of non-realisability of distance patterns.  Pure Python, integers only.

A pattern on n points is a list of colours, one for each pair (0,1), (0,2), ..., (n-2,n-1); the
colours are 0..k-1.  Suppose points p_0..p_{n-1} of the plane have this pattern: pairs of the
same colour are at the same distance, x_c > 0 is the squared distance of colour c, and
x_a != x_b for a != b.

For a base point b the numbers G_ij = x_c(b,i) + x_c(b,j) - x_c(i,j) (with G_ii = 2 x_c(b,i))
are twice the inner products of the vectors p_i - p_b and p_j - p_b.  These vectors lie in a
plane, so the matrix G has rank at most 2 and each of its 3x3 minors is zero.

A certificate consists of a non-zero integer `scale`, a list `product` of factors x_c or
(x_a - x_b) with exponents, and a list `terms` of minors (b, rows, columns), each with a
polynomial multiplier.  check() verifies the polynomial identity

        scale * product  =  sum of multiplier * minor.

At the point x the right side is 0 and the left side is not.  So no such points exist.

A polynomial is a dictionary {exponents: coefficient}.

    python3 check.py [--jobs K] FILE...        (files with one certificate per line)
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
    """True if the certificate is valid; an AssertionError if it is not."""
    n, pattern = cert["n"], cert["pattern"]
    pairs = list(combinations(range(n), 2))
    assert len(pattern) == len(pairs) and all(isinstance(c, int) for c in pattern)
    k = max(pattern) + 1
    assert set(pattern) == set(range(k))
    colour = {}
    for (i, j), c in zip(pairs, pattern):
        colour[i, j] = colour[j, i] = c

    def var(c):                               # the polynomial x_c
        return {tuple(int(i == c) for i in range(k)): 1}

    def entry(b, i, j):                       # G_ij for the base point b, a linear form
        if i == j:
            return add(var(colour[b, i]), var(colour[b, i]))
        return add(add(var(colour[b, i]), var(colour[b, j])), var(colour[i, j]), -1)

    def minor(b, rows, cols):                 # 3x3 determinant, expanded along the first row
        m = [[entry(b, i, j) for j in cols] for i in rows]

        def det2(c1, c2):                     # rows 1 and 2, columns c1 and c2
            return add(multiply(m[1][c1], m[2][c2]), multiply(m[1][c2], m[2][c1]), -1)

        out = multiply(m[0][0], det2(1, 2))
        out = add(out, multiply(m[0][1], det2(0, 2)), -1)
        return add(out, multiply(m[0][2], det2(0, 1)))

    # left side: scale * product of factors
    assert isinstance(cert["scale"], int) and cert["scale"] != 0
    left = {(0,) * k: cert["scale"]}
    for factor, exponent in cert["product"]:
        assert isinstance(exponent, int) and exponent >= 1 and len(factor) in (1, 2)
        assert all(isinstance(c, int) and 0 <= c < k for c in factor)
        if len(factor) == 1:
            form = var(factor[0])                              # x_c
        else:
            assert factor[0] != factor[1]
            form = add(var(factor[0]), var(factor[1]), -1)     # x_a - x_b
        for _ in range(exponent):
            left = multiply(left, form)

    # right side: sum of multiplier * minor
    right = {}
    for b, rows, cols, multiplier in cert["terms"]:
        assert isinstance(b, int) and 0 <= b < n
        for triple in (rows, cols):
            assert len(set(triple)) == 3 and b not in triple
            assert all(isinstance(i, int) and 0 <= i < n for i in triple)
        a = {}
        for coefficient, exponents in multiplier:
            assert isinstance(coefficient, int) and len(exponents) == k
            assert all(isinstance(e, int) and e >= 0 for e in exponents)
            a = add(a, {tuple(exponents): coefficient})
        right = add(right, multiply(a, minor(b, rows, cols)))

    assert left == right, "the identity does not hold"
    return True


def lines(name):
    """The lines of a file with one certificate (a JSON object) per line (.jsonl, .jsonl.gz)."""
    import gzip
    name = str(name)
    opener = gzip.open if name.endswith(".gz") else open
    with opener(name, "rt") as fh:
        for line in fh:
            if line.strip():
                yield line


def read(name):
    import json
    return [json.loads(line) for line in lines(name)]


def check_line(line):
    import json
    cert = json.loads(line)
    check(cert)
    return cert["n"]


if __name__ == "__main__":
    # one certificate at a time: a large file is never loaded whole
    import sys
    from collections import Counter
    args = sys.argv[1:]
    jobs = 1
    if args and args[0] == "--jobs":
        jobs, args = int(args[1]), args[2:]
    total = Counter()
    for name in args:
        if jobs > 1:
            from multiprocessing import Pool
            with Pool(jobs) as pool:
                count = Counter(pool.imap_unordered(check_line, lines(name), chunksize=500))
        else:
            count = Counter(map(check_line, lines(name)))
        total.update(count)
        print(f"{name}: {sum(count.values())} certificates, every identity holds", flush=True)
    print(sum(total.values()), "certificates checked; by number of points:",
          dict(sorted(total.items())))
