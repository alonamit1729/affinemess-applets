"""fewdist -- exact tools for planar point sets with few distinct distances.

Pure Python plus SymPy.  The notebook `two-shapes.ipynb` explains every idea used here;
this file is the single place where the code lives, so that the heavy runs and the small
examples in the notebook execute exactly the same functions.

Contents
    1. Patterns: edge-colourings of the complete graph K_n, canonical forms
    2. Generating patterns one point at a time
    3. Exact numbers: rationals and small number fields, with certified signs
    4. From a pattern to all of its exact realisations
    5. Extending a point set without creating too many new distances
    6. Classifying the point sets that use a prescribed list of distances
    7. Small helpers for saving results
"""
from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product

import sympy as sp

# =====================================================================================
# 1. PATTERNS
# =====================================================================================
# A *pattern* on n points is a colouring of the n(n-1)/2 pairs of points.  It is stored as
# a tuple of colour numbers, one per pair, the pairs being listed in the order
# (0,1), (0,2), ..., (0,n-1), (1,2), ..., (n-2,n-1).  Colours have no names of their own:
# only the partition of the pairs into colour classes matters.


@lru_cache(maxsize=None)
def pairs(n):
    """The pairs (i, j), i < j, of n points, in the fixed order used throughout."""
    return tuple(combinations(range(n), 2))


def n_points(pattern):
    """Number of points of a pattern (inverse of n -> n(n-1)/2)."""
    n = (1 + math.isqrt(1 + 8 * len(pattern))) // 2
    assert n * (n - 1) // 2 == len(pattern), "not a colouring of a complete graph"
    return n


def rename(colours):
    """Rename colours 0, 1, 2, ... in order of first appearance."""
    names = {}
    return tuple(names.setdefault(c, len(names)) for c in colours)


def squared_distance(p, q):
    """Squared Euclidean distance of two points given by exact coordinates."""
    return sum((a - b) ** 2 for a, b in zip(p, q))


def pattern_of(points, d2=squared_distance):
    """The distance pattern of a point set: equal colour <=> equal distance."""
    return rename(d2(points[i], points[j]) for i, j in pairs(len(points)))


def colour_matrix(pattern):
    """The pattern as a symmetric n x n table M[i][j] = colour of the pair {i, j}."""
    n = n_points(pattern)
    M = [[-1] * n for _ in range(n)]
    for (i, j), c in zip(pairs(n), pattern):
        M[i][j] = M[j][i] = c
    return M


def subpattern(pattern, points):
    """The pattern induced on a subset of the points."""
    M = colour_matrix(pattern)
    return rename(M[a][b] for a, b in combinations(points, 2))


# ---- canonical forms ---------------------------------------------------------------
# Two patterns are *equivalent* when one becomes the other after renumbering the points
# (and renaming the colours).  `canonical` returns the same tuple for all equivalent
# patterns and different tuples for inequivalent ones.
#
# Method ("individualise and refine", the idea behind nauty):
#   * refine: sort the points into classes by what they see -- how many pairs of each
#     colour, towards which classes -- and repeat until nothing changes;
#   * if some class still has several points, try each of them as "the first one" and
#     refine again;
#   * every branch ends with all points in different classes, i.e. with an ordering of
#     the points; the smallest colour string over all these orderings is the answer.
# Nothing in the procedure looks at the original numbering, so equivalent patterns
# produce the same set of strings.


def _refine(M, n, cls, named):
    k = 1 + max(max(row) for row in M)
    while True:
        if named:                       # colour numbers are meaningful: use them as they are
            crank = list(range(k))
        else:                           # colours are anonymous: describe each by the classes it joins
            seen = [[] for _ in range(k)]
            for i in range(n):
                ci, Mi = cls[i], M[i]
                for j in range(i + 1, n):
                    cj = cls[j]
                    seen[Mi[j]].append((ci, cj) if ci <= cj else (cj, ci))
            keys = [tuple(sorted(s)) for s in seen]
            rank = {key: r for r, key in enumerate(sorted(set(keys)))}
            crank = [rank[key] for key in keys]
        sig = [(cls[v], tuple(sorted((crank[M[v][w]], cls[w]) for w in range(n) if w != v)))
               for v in range(n)]
        rank = {s: r for r, s in enumerate(sorted(set(sig)))}
        if len(rank) == max(cls) + 1:   # no class was split: stable
            return cls
        cls = [rank[s] for s in sig]


def _orderings(M, n, named):
    """The orderings of the points at the leaves of the individualise-and-refine tree."""
    stack = [_refine(M, n, [0] * n, named)]
    while stack:
        cls = stack.pop()
        r = max(cls) + 1
        if r == n:
            order = [0] * n
            for v, c in enumerate(cls):
                order[c] = v
            yield order
            continue
        size = Counter(cls)
        t = min(c for c in range(r) if size[c] > 1)          # first class with several points
        for v in range(n):
            if cls[v] == t:                                    # v goes first within its class
                new = [c + 1 if (c > t or (c == t and u != v)) else c for u, c in enumerate(cls)]
                stack.append(_refine(M, n, new, named))


def _strings(pattern, named):
    n = n_points(pattern)
    M = colour_matrix(pattern)
    for order in _orderings(M, n, named):
        s = tuple(M[order[a]][order[b]] for a, b in pairs(n))
        yield s if named else rename(s)


def canonical(pattern, named=False):
    """Canonical form of a pattern.  With named=True the colour numbers are kept as they
    are (use this when the colours stand for actual lengths)."""
    return min(_strings(pattern, named))


def symmetry_count(pattern, named=False):
    """Number of renumberings of the points that map the pattern to itself."""
    strings = list(_strings(pattern, named))
    return strings.count(min(strings))


def same_pattern(p, q):
    return len(p) == len(q) and canonical(p) == canonical(q)


# =====================================================================================
# 2. GENERATING PATTERNS ONE POINT AT A TIME
# =====================================================================================

def colourings(m, max_colours, used=0):
    """All ways to colour m objects in a row with anonymous colours, at most max_colours
    in total, when `used` colours 0..used-1 are already in play.  A new colour always
    gets the next free number, so no colouring is produced twice under renaming."""
    def rec(prefix, top):
        if len(prefix) == m:
            yield tuple(prefix)
            return
        for c in range(min(top + 1, max_colours)):
            prefix.append(c)
            yield from rec(prefix, max(top, c + 1))
            prefix.pop()
    yield from rec([], used)


def all_patterns(n, max_colours=6):
    """Every pattern on n points with at most max_colours colours, one per equivalence
    class.  Brute force: only for n <= 4."""
    return sorted({canonical(s) for s in colourings(n * (n - 1) // 2, max_colours)})


def add_a_point(parents, allowed, max_colours=6):
    """All patterns on n+1 points that (a) extend one of the `parents` (patterns on n points)
    by a new point and (b) have every n-point sub-pattern in `allowed` (canonical forms).
    This is one step of the hereditary generation of Section 3 of the notebook."""
    allowed = set(allowed)
    n = n_points(parents[0])
    old_pairs, new_pairs = pairs(n), pairs(n + 1)
    index = {e: k for k, e in enumerate(old_pairs)}
    # for each sub-pattern obtained by dropping one old point: where its pairs come from
    drops = [[e for e in combinations([v for v in range(n + 1) if v != d], 2)] for d in range(n)]
    table = {}                                  # labelled sub-pattern -> allowed?  (memo)
    children = set()
    for P in parents:
        for new in colourings(n, max_colours, used=max(P) + 1):
            colour = lambda a, b: new[a] if b == n else P[index[(a, b)]]
            ok = True
            for sub_pairs in drops:
                sub = rename(colour(a, b) for a, b in sub_pairs)
                if sub not in table:
                    table[sub] = canonical(sub) in allowed
                if not table[sub]:
                    ok = False
                    break
            if ok:
                children.add(canonical(tuple(colour(a, b) for a, b in new_pairs)))
    return sorted(children)


# =====================================================================================
# 3. EXACT NUMBERS
# =====================================================================================
# Lengths are either rational numbers (Python Fractions) or elements of a number field
# Q(t), t a root of an irreducible polynomial f.  An element of Q(t) is a polynomial in t
# of degree < deg f with rational coefficients; arithmetic is polynomial arithmetic
# followed by reduction modulo f.  Such arithmetic decides *equality* exactly.  To decide
# *signs* one must say which real root of f is meant; `sign` does this with interval
# arithmetic on an isolating interval that is narrowed until the answer is certain.


class NumberField:
    def __init__(self, f):
        """f: coefficients of a monic irreducible polynomial, constant term first."""
        self.f = tuple(Fraction(c) for c in f)
        assert self.f[-1] == 1 and len(self.f) >= 3
        self.deg = d = len(self.f) - 1
        # t^d, t^(d+1), ..., t^(2d-2) written in terms of 1, t, ..., t^(d-1)
        self._high = []
        cur = [-c for c in self.f[:-1]]
        for _ in range(d - 1):
            self._high.append(tuple(cur))
            top, cur = cur[-1], [Fraction(0)] + cur[:-1]
            cur = [a - top * b for a, b in zip(cur, self.f[:-1])]
        self.one = Alg(self, (Fraction(1),) + (Fraction(0),) * (d - 1))
        self.t = Alg(self, (Fraction(0), Fraction(1)) + (Fraction(0),) * (d - 2))
        poly = sp.Poly([sp.Rational(c.numerator, c.denominator) for c in reversed(self.f)], sp.Symbol('t'))
        self.sympy_poly = poly
        # isolating intervals of the real roots, in increasing order
        self.intervals = [[Fraction(int(a.p), int(a.q)), Fraction(int(b.p), int(b.q))]
                          for (a, b), _ in poly.intervals()]

    def element(self, coefficients):
        c = tuple(Fraction(x) for x in coefficients)
        return Alg(self, c + (Fraction(0),) * (self.deg - len(c)))

    @property
    def real_roots(self):
        return len(self.intervals)

    def _eval_f(self, x):
        v = Fraction(0)
        for c in reversed(self.f):
            v = v * x + c
        return v

    def _narrow(self, i):
        lo, hi = self.intervals[i]
        mid = (lo + hi) / 2
        fm = self._eval_f(mid)
        if fm == 0:                      # cannot happen for irreducible f of degree >= 2
            raise ArithmeticError("rational root")
        if (self._eval_f(lo) > 0) == (fm > 0):
            self.intervals[i][0] = mid
        else:
            self.intervals[i][1] = mid

    def enclose(self, a, i):
        """An interval [L, H] that contains the value of a at the i-th real root of f
        (Horner's rule in interval arithmetic on the current isolating interval)."""
        lo, hi = self.intervals[i]
        L = H = Fraction(0)
        for c in reversed(a.c):
            ends = (L * lo, L * hi, H * lo, H * hi)
            L, H = min(ends) + c, max(ends) + c
        return L, H

    def sign(self, a, i):
        """Certified sign of the element a at the i-th real root of f."""
        if not isinstance(a, Alg):
            return (a > 0) - (a < 0)
        if not a:
            return 0
        while True:
            L, H = self.enclose(a, i)
            if L > 0:
                return 1
            if H < 0:
                return -1
            self._narrow(i)

    def approx(self, a, i, digits=15):
        """Floating-point value of a at the i-th real root."""
        eps = Fraction(1, 10 ** (digits + 2))
        while self.intervals[i][1] - self.intervals[i][0] > eps:
            self._narrow(i)
        x = (self.intervals[i][0] + self.intervals[i][1]) / 2
        if not isinstance(a, Alg):
            return float(a)
        v = Fraction(0)
        for c in reversed(a.c):
            v = v * x + c
        return float(v)

    def __repr__(self):
        return f"Q(t), t a root of {self.sympy_poly.as_expr()}"


class Alg:
    """An element of a NumberField: a polynomial in t of degree < deg f."""
    __slots__ = ("K", "c", "_hash")

    def __init__(self, K, c):
        self.K, self.c, self._hash = K, c, None

    def _coerce(self, o):
        if isinstance(o, Alg):
            return o.c
        return (Fraction(o),) + (Fraction(0),) * (self.K.deg - 1)

    def __add__(self, o):
        return Alg(self.K, tuple(a + b for a, b in zip(self.c, self._coerce(o))))
    __radd__ = __add__

    def __sub__(self, o):
        return Alg(self.K, tuple(a - b for a, b in zip(self.c, self._coerce(o))))

    def __rsub__(self, o):
        return Alg(self.K, tuple(b - a for a, b in zip(self.c, self._coerce(o))))

    def __neg__(self):
        return Alg(self.K, tuple(-a for a in self.c))

    def __mul__(self, o):
        if not isinstance(o, Alg):
            o = Fraction(o)
            return Alg(self.K, tuple(a * o for a in self.c))
        d = self.K.deg
        prod = [Fraction(0)] * (2 * d - 1)
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    if b:
                        prod[i + j] += a * b
        low = prod[:d]
        for k in range(d, 2 * d - 1):
            if prod[k]:
                for i, h in enumerate(self.K._high[k - d]):
                    low[i] += prod[k] * h
        return Alg(self.K, tuple(low))
    __rmul__ = __mul__

    def inverse(self):
        """Solve a * x = 1 as a linear system for the coefficients of x."""
        d = self.K.deg
        cols, power = [], self                      # columns: a, a*t, a*t^2, ...
        for _ in range(d):
            cols.append(power.c)
            power = power * self.K.t
        rows = [[cols[j][i] for j in range(d)] + [Fraction(int(i == 0))] for i in range(d)]
        for col in range(d):                        # Gauss-Jordan elimination
            piv = next(r for r in range(col, d) if rows[r][col] != 0)
            rows[col], rows[piv] = rows[piv], rows[col]
            p = rows[col][col]
            rows[col] = [x / p for x in rows[col]]
            for r in range(d):
                if r != col and rows[r][col] != 0:
                    m = rows[r][col]
                    rows[r] = [x - m * y for x, y in zip(rows[r], rows[col])]
        return Alg(self.K, tuple(rows[i][d] for i in range(d)))

    def __truediv__(self, o):
        if isinstance(o, Alg):
            return self * o.inverse()
        return Alg(self.K, tuple(a / Fraction(o) for a in self.c))

    def __rtruediv__(self, o):
        return self.inverse() * o

    def __eq__(self, o):
        return self.c == self._coerce(o)

    def __hash__(self):
        if self._hash is None:
            rational = not any(self.c[1:])
            self._hash = hash(self.c[0]) if rational else hash(self.c)
        return self._hash

    def __bool__(self):
        return any(self.c)

    def __repr__(self):
        return str(sum(sp.Rational(a.numerator, a.denominator) * sp.Symbol('t') ** k
                       for k, a in enumerate(self.c)))


def field_of(x):
    """The NumberField an exact number lives in (None for a rational)."""
    return x.K if isinstance(x, Alg) else None


# =====================================================================================
# 4. FROM A PATTERN TO ITS EXACT REALISATIONS
# =====================================================================================
# Unknowns: the squared lengths.  The colour of the pair {0,1} is the unit (squared length
# 1); the other colours get variables.  For points p_0, ..., p_{n-1} with squared
# distances D_ij the vectors v_i = p_i - p_0 have dot products
#       <v_i, v_j> = (D_0i + D_0j - D_ij) / 2.
# The points lie in a plane exactly when this (n-1) x (n-1) "Gram matrix" has rank <= 2,
# and then every 3 x 3 sub-determinant vanishes.


def unit_colour(pattern):
    return pattern[0]


def length_symbols(pattern):
    """One symbol per colour other than the unit colour; returns (symbols, value per colour)."""
    k = max(pattern) + 1
    xs = sp.symbols(f"x1:{k}") if k > 1 else ()
    it = iter(xs)
    values = [sp.Integer(1) if c == unit_colour(pattern) else next(it) for c in range(k)]
    return tuple(xs), values


def distance_table(pattern, values):
    """D[i][j] = value of the colour of the pair {i, j} (0 on the diagonal)."""
    n = n_points(pattern)
    zero = values[0] - values[0]
    D = [[zero] * n for _ in range(n)]
    for (i, j), c in zip(pairs(n), pattern):
        D[i][j] = D[j][i] = values[c]
    return D


def gram2(D, base=0):
    """Twice the Gram matrix of the vectors p_i - p_base."""
    others = [i for i in range(len(D)) if i != base]
    return [[D[base][i] + D[base][j] - D[i][j] for j in others] for i in others]


def _det3(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def planarity_equations(pattern, minors="principal"):
    """Polynomial equations in the unknown squared lengths that every planar realisation
    satisfies: 3 x 3 minors of the Gram matrix based at point 0.

    minors="principal": only the minors on rows = columns {i,j,k}; each says that the four
        points 0, i, j, k are coplanar (it is their Cayley-Menger determinant).  These are
        necessary conditions, and enough for our purposes because every candidate
        solution is checked afterwards by rebuilding the points.
    minors="all": every 3 x 3 minor, i.e. exactly "rank <= 2" (Szollosi-Ostergard's system).
    Returns (symbols, list of polynomials); a polynomial may be a non-zero constant, which
    means the pattern is impossible."""
    xs, values = length_symbols(pattern)
    if xs:
        ring = sp.ring(xs, sp.QQ)[0]
        values = [ring.from_expr(v) for v in values]
    else:
        values = [sp.Integer(1) for _ in values]
    G = gram2(distance_table(pattern, values))
    m = len(G)
    triples = list(combinations(range(m), 3))
    index_sets = ([(I, I) for I in triples] if minors == "principal"
                  else [(I, J) for I in triples for J in triples if I <= J])
    found = set()
    for I, J in index_sets:
        d = _det3([[G[i][j] for j in J] for i in I])
        if d != 0:
            found.add(d.monic() if xs else sp.Integer(1))
    polys = sorted((p.as_expr() if xs else p for p in found), key=lambda e: (sp.count_ops(e), str(e)))
    return xs, polys


def _nondegeneracy_factors(pattern):
    """Linear factors whose product must be non-zero: every length, and every difference."""
    xs, values = length_symbols(pattern)
    return list(xs) + [a - b for a, b in combinations(values, 2)]


def _saturated(pattern, minors):
    """Groebner basis of <equations, w*q - 1> where q is the product of the non-degeneracy
    factors (reduced along the way).  Returns (symbols, w, basis) or None if no solutions."""
    xs, eqs = planarity_equations(pattern, minors)
    if any(e.is_number for e in eqs):            # a non-zero constant: 1 = 0
        return None
    if not xs:
        return (), None, []
    w = sp.Symbol("w")
    if eqs:
        G = sp.groebner(eqs, *xs, order="grevlex", domain=sp.QQ)
        if G.exprs == [1]:
            return None
        q = sp.Integer(1)
        for factor in _nondegeneracy_factors(pattern):
            q = G.reduce(sp.expand(q * factor))[1]
            if q == 0:                           # the product vanishes on every solution
                return None
        base = list(G.exprs)
    else:
        q, base = sp.expand(sp.prod(_nondegeneracy_factors(pattern))), []
    S = sp.groebner(base + [w * q - 1], w, *xs, order="grevlex", domain=sp.QQ)
    if S.exprs == [1]:
        return None
    return xs, w, S


def complex_realisable(pattern):
    """Szollosi-Ostergard's test: is there a (possibly complex) assignment of pairwise
    different non-zero squared lengths making the Gram matrix of rank <= 2?

    With a single equation F = 0 (always the case for four points) the answer is read off
    the factorisation of F: the solution set lies inside the forbidden hyperplanes exactly
    when F is, up to a constant, a product of the forbidden linear factors."""
    xs, eqs = planarity_equations(pattern, "all")
    if len(eqs) == 1 and not eqs[0].is_number:
        forbidden = [sp.Poly(f, *xs) for f in _nondegeneracy_factors(pattern)]
        for g, _ in sp.factor_list(eqs[0], *xs)[1]:
            g = sp.Poly(g, *xs)
            if not any(g.monic() == f.monic() for f in forbidden):
                return True
        return False
    return _saturated(pattern, "all") is not None


def _solve_zero_dimensional(B, xs):
    """All solutions of a lex Groebner basis B with finitely many solutions.

    Linear algebra in the quotient ring Q[xs]/<B>: its dimension N is the number of
    solutions (B is radical in every case we meet; this is asserted).  Pick a linear form
    L taking N different values on them; then 1, L, ..., L^(N-1) is a basis, the minimal
    polynomial m of L has the values of L as roots, and every variable is a polynomial in L.
    Returns [(m_i, [variable as polynomial in t]), ...], one entry per irreducible factor."""
    lead = [sp.Poly(b, *xs).monoms()[0] for b in B]        # leading monomials (lex order)
    std, frontier = [tuple([0] * len(xs))], [tuple([0] * len(xs))]
    while frontier:                               # monomials not divisible by a leading monomial
        new = []
        for mono in frontier:
            for k in range(len(xs)):
                m2 = mono[:k] + (mono[k] + 1,) + mono[k + 1:]
                if m2 not in std and not any(all(a >= b for a, b in zip(m2, l)) for l in lead):
                    std.append(m2)
                    new.append(m2)
        frontier = new
        if len(std) > 200:
            raise ValueError("infinitely many solutions?")
    N = len(std)

    def coords(expr):
        r = sp.Poly(sp.reduced(sp.expand(expr), B, *xs, order="lex")[1], *xs)
        d = dict(r.terms())
        assert set(d) <= set(std)
        return [d.get(e, sp.Integer(0)) for e in std]

    t = sp.Symbol("t")
    forms = [x for x in reversed(xs)] + [sum((c ** j) * x for j, x in enumerate(xs)) for c in (2, 3, 5, 7)]
    for L in forms:
        powers = [coords(L ** k) for k in range(N + 1)]
        V = sp.Matrix(N, N, lambda r, c: powers[c][r])
        if V.det() == 0:
            continue                              # L does not separate the solutions
        top = V.LUsolve(sp.Matrix(powers[N]))
        m = sp.Poly(t ** N - sum(top[k] * t ** k for k in range(N)), t)
        if sp.gcd(m, m.diff(t)).degree() > 0:
            continue
        as_poly = [list(V.LUsolve(sp.Matrix(coords(x)))) for x in xs]
        out = []
        for f, mult in m.factor_list()[1]:
            f = f.monic()
            out.append((f, [sp.Poly(sum(c * t ** k for k, c in enumerate(col)), t).rem(f) for col in as_poly]))
        return out
    raise ValueError("no separating linear form found")


class LengthSolution:
    """One solution of the length equations of a pattern, up to the choice of a real root:
    `lengths[c]` is the squared length of colour c, a Fraction or an element of `field`."""

    def __init__(self, pattern, lengths):
        self.pattern = tuple(pattern)
        self.lengths = tuple(lengths)
        self.field = next((x.K for x in self.lengths if isinstance(x, Alg)), None)

    @property
    def embeddings(self):
        """Indices of the real roots at which this is a genuine planar point set."""
        if self.field is None:
            return [0] if self.realisation(0) is not None else []
        return [i for i in range(self.field.real_roots) if self.realisation(i) is not None]

    def realisation(self, root=0):
        """The point set at the given real root, or None if the lengths are not all positive
        there or the plane is not Euclidean there."""
        K = self.field
        sign = (lambda a: K.sign(a, root)) if K else (lambda a: (a > 0) - (a < 0))
        if any(sign(v) <= 0 for v in self.lengths):
            return None
        conf = place(self.pattern, self.lengths)
        if conf is None or sign(conf.A * conf.B - conf.C * conf.C) <= 0:
            return None
        conf.root = root
        return conf

    def to_json(self):
        s = lambda q: str(Fraction(q))
        if self.field is None:
            return {"f": None, "lengths": [s(v) for v in self.lengths]}
        K = self.field
        return {"f": [s(c) for c in K.f],
                "lengths": [[s(c) for c in (v.c if isinstance(v, Alg) else K.element([v]).c)] for v in self.lengths]}

    @staticmethod
    def from_json(pattern, d, fields=None):
        if d["f"] is None:
            return LengthSolution(pattern, [Fraction(v) for v in d["lengths"]])
        key = tuple(d["f"])
        fields = fields if fields is not None else {}
        K = fields.get(key) or fields.setdefault(key, NumberField([Fraction(c) for c in d["f"]]))
        return LengthSolution(pattern, [K.element([Fraction(c) for c in v]) for v in d["lengths"]])


_MINIMAL = {}


def real_number_id(x, root=0):
    """A name for the real number "x at the given root" that does not depend on how x is
    written: ("Q", value) for a rational number, otherwise (minimal polynomial of x over Q,
    index of x among the real roots of that polynomial).  Two real algebraic numbers are
    equal exactly when their names are equal."""
    if not isinstance(x, Alg):
        return ("Q", Fraction(x))
    if not any(x.c[1:]):
        return ("Q", x.c[0])
    K = x.K
    if (K.f, x.c) not in _MINIMAL:
        cols, power = [], K.one                      # matrix of "multiply by x" on 1, t, t^2, ...
        for _ in range(K.deg):
            cols.append((x * power).c)
            power = power * K.t
        M = sp.Matrix(K.deg, K.deg, lambda i, j: sp.Rational(cols[j][i].numerator, cols[j][i].denominator))
        t = sp.Symbol("t")
        (mu, _), = sp.Poly(M.charpoly(t).as_expr(), t).factor_list()[1]   # characteristic polynomial = mu^k
        mu = mu.monic()
        frac = lambda q: Fraction(int(q.p), int(q.q))
        _MINIMAL[K.f, x.c] = (tuple(frac(c) for c in mu.all_coeffs()),
                              [(frac(a), frac(b)) for (a, b), _ in mu.intervals()])
    coeffs, isolating = _MINIMAL[K.f, x.c]
    while True:                                      # which real root of mu is x?  narrow until certain
        L, H = K.enclose(x, root)
        hits = [k for k, (a, b) in enumerate(isolating) if not (H < a or L > b)]
        if len(hits) == 1:
            return coeffs, hits[0]
        K._narrow(root)


def shape_key(solution, root=0):
    """A key that is equal for two realisations exactly when the point sets are similar:
    the lengths divided by the shortest one, in increasing order, together with the canonical
    form of the pattern in which every colour is named by the rank of its length."""
    from functools import cmp_to_key
    K, lengths = solution.field, solution.lengths
    sign = (lambda a: K.sign(a, root)) if K else (lambda a: (a > 0) - (a < 0))
    order = sorted(range(len(lengths)), key=cmp_to_key(lambda a, b: sign(lengths[a] - lengths[b])))
    rank = {c: r for r, c in enumerate(order)}
    ratios = tuple(real_number_id(lengths[c] / lengths[order[0]], root) for c in order)
    return ratios, canonical(tuple(rank[c] for c in solution.pattern), named=True)


def solve_lengths(pattern):
    """Every assignment of pairwise different, non-zero squared lengths (unit colour = 1)
    that satisfies the planarity equations of the pattern, as a list of LengthSolution.
    One LengthSolution stands for deg(f) numerical solutions, one per root of its field.
    Raises ValueError if there are infinitely many (a flexible pattern)."""
    sat = _saturated(pattern, "principal")
    if sat is None:
        return []
    xs, w, S = sat
    k = max(pattern) + 1
    if not xs:
        return [LengthSolution(pattern, [Fraction(1)] * k)]
    if not S.is_zero_dimensional:
        raise ValueError("infinitely many solutions: the pattern is flexible")
    B = [g for g in S.fglm("lex").exprs if not g.has(w)]      # eliminate the auxiliary variable w
    solutions = []
    for f, polys in _solve_zero_dimensional(B, xs):
        coeffs = lambda p: [Fraction(int(c.p), int(c.q)) for c in reversed(p.all_coeffs())]
        if f.degree() == 1:
            root = -f.all_coeffs()[1] / f.all_coeffs()[0]
            vals = iter(Fraction(int(v.p), int(v.q)) for v in (p.eval(root) for p in polys))
            unit = Fraction(1)
        else:
            K = NumberField(coeffs(f))
            vals = iter(K.element(coeffs(p)) for p in polys)
            unit = K.one
        lengths = [unit if c == unit_colour(pattern) else next(vals) for c in range(k)]
        solutions.append(LengthSolution(pattern, lengths))
    return solutions


# ---- from lengths back to points ------------------------------------------------------

class Configuration:
    """A point set in exact coordinates with respect to a basis e1, e2 of the plane that need
    not be orthonormal: <e1,e1> = A, <e2,e2> = B, <e1,e2> = C.  The squared length of the
    vector (x, y) is A x^2 + 2 C x y + B y^2.  Everything is exact."""

    def __init__(self, points, A, B, C, root=0):
        self.points, self.A, self.B, self.C, self.root = list(points), A, B, C, root

    def dot(self, u, v):
        return self.A * u[0] * v[0] + self.B * u[1] * v[1] + self.C * (u[0] * v[1] + u[1] * v[0])

    def d2(self, p, q):
        w = (p[0] - q[0], p[1] - q[1])
        return self.dot(w, w)

    def distances(self):
        return {self.d2(p, q) for p, q in combinations(self.points, 2)}

    def collinear(self, p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) == (q[1] - p[1]) * (r[0] - p[0])

    def with_points(self, extra):
        return Configuration(self.points + list(extra), self.A, self.B, self.C, self.root)

    def xy(self, points=None):
        """Ordinary (floating-point) coordinates, for drawing."""
        K = field_of(self.A) or field_of(self.C) or next(
            (field_of(c) for p in self.points for c in p if field_of(c)), None)
        val = (lambda a: K.approx(a, self.root)) if K else float
        A, B, C = val(self.A), val(self.B), val(self.C)
        e1, e2 = (math.sqrt(A), 0.0), (C / math.sqrt(A), math.sqrt(max(B - C * C / A, 0.0)))
        pts = self.points if points is None else points
        return [(val(x) * e1[0] + val(y) * e2[0], val(y) * e2[1]) for x, y in pts]


def place(pattern, lengths):
    """Rebuild the points from the squared lengths, or return None if they do not fit in a
    plane.  Point 0 is the origin; two other points a, b with independent position vectors
    serve as basis; every point is then found from its dot products with those two; finally
    ALL squared distances are compared with the prescribed ones."""
    n = n_points(pattern)
    D = distance_table(pattern, list(lengths))
    G = [[(D[0][i] + D[0][j] - D[i][j]) / 2 for j in range(n)] for i in range(n)]
    basis = next(((a, b) for a, b in combinations(range(1, n), 2)
                  if G[a][a] * G[b][b] != G[a][b] * G[a][b]), None)
    if basis is None:
        return None                                 # all points on a line
    a, b = basis
    A, B, C = G[a][a], G[b][b], G[a][b]
    H = A * B - C * C
    pts = [((B * G[a][j] - C * G[b][j]) / H, (A * G[b][j] - C * G[a][j]) / H) for j in range(n)]
    conf = Configuration(pts, A, B, C)
    if any(conf.d2(pts[i], pts[j]) != D[i][j] for i, j in pairs(n)):
        return None
    return conf


# =====================================================================================
# 5. EXTENDING A POINT SET
# =====================================================================================

def anchor_triples(conf, k):
    """A short list of non-collinear triples of points of `conf` such that every k of its
    points contain one of the triples (chosen greedily)."""
    pts, n = conf.points, len(conf.points)
    triples = [T for T in combinations(range(n), 3) if not conf.collinear(*(pts[i] for i in T))]
    mask = lambda T: sum(1 << i for i in T)
    masks = [mask(T) for T in triples]
    todo = {mask(S) for S in combinations(range(n), k)}
    chosen = []
    while todo:
        best = max(range(len(triples)), key=lambda j: sum((masks[j] & W) == masks[j] for W in todo))
        covered = {W for W in todo if (masks[best] & W) == masks[best]}
        if not covered:
            raise ValueError(f"some {k} points are collinear: the completeness argument fails")
        todo -= covered
        chosen.append(triples[best])
    return chosen


def _locate(conf, T, r, s, t):
    """The unique point whose squared distances to the three points T = (u, v, w) differ
    like r, s, t:  2<y, v-u> = |v|^2 - |u|^2 + r - s  and the same with w, t."""
    u, v, w = (conf.points[i] for i in T)
    rows = []
    for p in (v, w):
        e = (p[0] - u[0], p[1] - u[1])
        rows.append((conf.A * e[0] + conf.C * e[1], conf.C * e[0] + conf.B * e[1]))
    (a, b), (c, d) = rows
    det = a * d - b * c
    zero = (0, 0)
    quv, quw = conf.d2(v, zero) - conf.d2(u, zero), conf.d2(w, zero) - conf.d2(u, zero)
    r1, r2 = (quv + r - s) / 2, (quw + r - t) / 2
    return ((d * r1 - b * r2) / det, (a * r2 - c * r1) / det)


def candidate_points(conf, budget, exhaustive=False):
    """Every point of the plane that can be added to `conf` creating at most `budget` new
    squared distances *to the points of conf*.  Returns {point: set of new distances}.

    Why the list is complete (n = number of points of conf, at most 5 of them collinear is
    checked by anchor_triples):  look at the n distances from a new point y.
      * If some value occurs three times, y is the circumcentre of three points of conf.
      * Otherwise each new value occurs at most twice, so at least n - 2*budget distances
        are old; any n - 2*budget points of conf contain an anchor triple T, and y is
        determined by its (old) distances to T.
    exhaustive=True uses all non-collinear triples instead of the anchor triples (slower;
    the result must be the same, which is a test of the argument)."""
    old = conf.distances()
    lengths = sorted(old, key=repr)
    pts, n = conf.points, len(conf.points)
    seen, good = set(pts), {}
    all_triples = [T for T in combinations(range(n), 3) if not conf.collinear(*(pts[i] for i in T))]
    anchors = all_triples if exhaustive else anchor_triples(conf, n - 2 * budget)

    def consider(y):
        if y in seen:
            return
        seen.add(y)
        new = set()
        for p in pts:
            d = conf.d2(y, p)
            if d not in old:
                new.add(d)
                if len(new) > budget:
                    return
        good[y] = frozenset(new)

    zero = lengths[0] - lengths[0]
    for T in all_triples:                          # circumcentres
        consider(_locate(conf, T, zero, zero, zero))
    for T in anchors:                              # points at old distances from an anchor triple
        for r, s, t in product(lengths, repeat=3):
            consider(_locate(conf, T, r, s, t))
    return good


def grow(conf, budget, exhaustive=False):
    """All ways to add points to `conf` so that at most `budget` new squared distances appear
    in total.  Returns (candidates, subsets): `candidates` is the list of addable points and
    `subsets` the list of all index-tuples of candidates that can be added together."""
    old = conf.distances()
    cand = list(candidate_points(conf, budget, exhaustive).items())
    subsets = []

    def search(start, chosen, new):
        subsets.append(tuple(chosen))
        for j in range(start, len(cand)):
            y, new_y = cand[j]
            total = new | new_y
            if len(total) > budget:
                continue
            for c in chosen:
                d = conf.d2(y, cand[c][0])
                if d not in old:
                    total = total | {d}
            if len(total) <= budget:
                search(j + 1, chosen + [j], total)

    search(0, [], frozenset())
    return [y for y, _ in cand], subsets


def proportional(values, target):
    """Is the set `values` equal to c * target for some constant c?"""
    values = set(values)
    return len(values) == len(target) and any(values == {u * m for m in target} for u in values)


# =====================================================================================
# 6. POINT SETS WITH A PRESCRIBED LIST OF SQUARED DISTANCES
# =====================================================================================

def sets_with_distances(palette, size):
    """Every planar set of `size` points all of whose squared distances lie in `palette`
    (a list of positive rationals), up to congruence.  Returns {canonical form: Configuration}.

    Such a set is not collinear (size > number of distances + 1), so it contains a genuine
    triangle; its squared sides (a <= b <= c) are in the palette.  Put the triangle at
    (0,0), (1,0), (0,1) in the basis of two of its sides.  Any further point is fixed by its
    three squared distances (r, s, t) to the corners -- two linear equations -- so there are at
    most |palette|^3 candidates per triangle.  Then look for size-3 candidates that are
    pairwise at allowed distances."""
    palette = [Fraction(x) for x in palette]
    allowed = set(palette)
    rank = {v: i for i, v in enumerate(sorted(allowed))}
    found = {}
    for a, b, c in combinations_with_replacement_sorted(palette):
        C = (a + b - c) / 2
        if a * b - C * C <= 0:
            continue                                   # not a triangle
        frame = Configuration([(Fraction(0), Fraction(0)), (Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))], a, b, C)
        cand = []
        for r, s, t in product(palette, repeat=3):
            y = _locate(frame, (0, 1, 2), r, s, t)
            if frame.d2(y, frame.points[0]) == r:      # then the distances to the other corners are s, t
                cand.append(y)
        later = [0] * len(cand)                        # bit j set: candidates i < j are compatible
        for i, j in combinations(range(len(cand)), 2):
            if frame.d2(cand[i], cand[j]) in allowed:
                later[i] |= 1 << j
        need = size - 3

        def cliques(pool, chosen):
            if len(chosen) == need:
                yield list(chosen)
                return
            while pool.bit_count() + len(chosen) >= need:
                low = pool & -pool
                j = low.bit_length() - 1
                pool ^= low
                yield from cliques(pool & later[j], chosen + [j])

        for clique in cliques((1 << len(cand)) - 1, []):
            conf = frame.with_points(cand[j] for j in clique)
            key = canonical(tuple(rank[conf.d2(p, q)] for p, q in combinations(conf.points, 2)), named=True)
            found.setdefault(key, conf)
    return found


def combinations_with_replacement_sorted(values):
    values = sorted(values)
    return [(a, b, c) for i, a in enumerate(values) for j, b in enumerate(values[i:], i) for c in values[j:]]


def congruence_key(conf):
    """Canonical form of the full table of squared distances: equal for two point sets with
    rational squared distances exactly when they are congruent."""
    values = sorted(conf.distances())
    rank = {v: i for i, v in enumerate(values)}
    return tuple(values), canonical(tuple(rank[conf.d2(p, q)] for p, q in combinations(conf.points, 2)), named=True)


# =====================================================================================
# 7. SAVING AND LOADING
# =====================================================================================

def load_patterns(path):
    with open(path) as fh:
        return [tuple(p) for p in json.load(fh)["patterns"]]


def save_solutions(path, solutions_by_pattern):
    with open(path, "w") as fh:
        json.dump([[s.to_json() for s in sols] for sols in solutions_by_pattern], fh, separators=(",", ":"))


def load_solutions(path, patterns):
    fields = {}
    with open(path) as fh:
        data = json.load(fh)
    return [[LengthSolution.from_json(p, d, fields) for d in sols] for p, sols in zip(patterns, data)]


def solve_many(patterns, processes=None):
    """solve_lengths for a list of patterns, in parallel."""
    import multiprocessing as mp
    with mp.Pool(processes) as pool:
        raw = pool.map(_solve_json, patterns, chunksize=1)
    fields = {}
    return [[LengthSolution.from_json(p, d, fields) for d in sols] for p, sols in zip(patterns, raw)]


def _solve_json(pattern):
    return [s.to_json() for s in solve_lengths(pattern)]
