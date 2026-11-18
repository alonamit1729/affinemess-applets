"""Fifteen points in the plane with at most seven distances: all of them, up to similarity.

Every such set other than the regular 15-gon contains nine points with at most six distances
(paper, section 8).  This program takes every nine-point set with at most six distances and
finds every way to add six points without exceeding seven distances.

  * The 635 nine-point patterns of certificates/L9.json are solved (fewdist.solve_lengths);
    the solutions are saved in solutions9.json.  631 patterns have finitely many solutions,
    two (124 and 520) have a one-parameter family of them.
  * Each rigid realisation K is extended exactly (candidate points, then a depth-first search).
  * Each family K(s) is analysed with s symbolic: no candidate point works for all s, so a
    fifteen-point set can contain K(s0) only for s0 among the roots of finitely many
    polynomials; K(s0) is then extended exactly for each of them.

All fifteen-point sets found have squared distances proportional to 1, 3, 4, 7, 9, 12, 13, and
fewdist.sets_with_distances lists those.  Result: 15 sets, plus the regular 15-gon.

    python3 fifteen.py [--jobs K] [--resolve]      (about 25 minutes on two cores)
"""
import json
import multiprocessing as mp
import os
import sys
import time
from collections import Counter
from fractions import Fraction
from functools import cmp_to_key
from itertools import combinations, product

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import fewdist as fd                                            # noqa: E402

NEED = 6                                                         # 9 + 6 = 15 points
PALETTE = [1, 3, 4, 7, 9, 12, 13]


# ---- shapes ---------------------------------------------------------------------------------

def shape_key(conf, root, field):
    """Equal for two point sets exactly when they are similar (as fewdist.shape_key)."""
    sign = (lambda a: field.sign(a, root)) if field else (lambda a: (a > 0) - (a < 0))
    vals = sorted(conf.distances(), key=cmp_to_key(lambda a, b: sign(a - b)))
    rank = {v: r for r, v in enumerate(vals)}
    ratios = tuple(fd.real_number_id(v / vals[0], root) for v in vals)
    pattern = tuple(rank[conf.d2(p, q)] for p, q in combinations(conf.points, 2))
    return ratios, fd.canonical(pattern, named=True)


# ---- candidate points -----------------------------------------------------------------------

def _bisector(conf, p, q):
    """The perpendicular bisector of p, q as a row (a, b, c): a*y0 + b*y1 = c."""
    e = (q[0] - p[0], q[1] - p[1])
    z = (0, 0)
    return (2 * (conf.A * e[0] + conf.C * e[1]), 2 * (conf.C * e[0] + conf.B * e[1]),
            conf.d2(q, z) - conf.d2(p, z))


def _meet(row1, row2):
    (a, b, c), (d, e, f) = row1, row2
    det = a * e - b * d
    return None if det == 0 else ((e * c - b * f) / det, (a * f - d * c) / det)


def _matchings(items):
    if not items:
        yield []
        return
    for k in range(1, len(items)):
        for m in _matchings(items[1:k] + items[k + 1:]):
            yield [(items[0], items[k])] + m


def candidate_points(conf, budget):
    """Every point that can be added to conf with at most `budget` new squared distances to it.

    Let y be such a point.  If a value occurs three times among its distances to conf, y is
    the circumcentre of three points (not collinear).  Otherwise at most 2*budget distances
    are new, so the set S of points of conf at old distances has at least n - 2*budget points.
      * If S contains three non-collinear points, y is determined by its old distances to them
        (two linear equations).  When no n - 2*budget points of conf are collinear, a short
        list of anchor triples suffices (fewdist.anchor_triples); otherwise all triples are
        used.
      * If S is collinear, then |S| = n - 2*budget (checked: conf has no more collinear points)
        and the other 2*budget points pair up at the `budget` new values, each taken twice;
        y lies on the perpendicular bisectors of these pairs.  Two bisectors that meet give y.
        If all of them are parallel and distinct there is no such y; if they coincide, y is on
        that line and on circles about the points of S, which is solved exactly below.
    """
    P = conf.points
    n = len(P)
    old = conf.distances()
    vals = sorted(old, key=repr)
    zero = vals[0] - vals[0]
    good, seen = {}, set(P)

    def consider(y):
        if y is None or y in seen:
            return
        seen.add(y)
        new = set()
        for p in P:
            d = conf.d2(y, p)
            if d not in old:
                new.add(d)
                if len(new) > budget:
                    return
        good[y] = frozenset(new)

    noncol = [T for T in combinations(range(n), 3) if not conf.collinear(*(P[t] for t in T))]
    lines = {}
    for T in combinations(range(n), 3):
        if conf.collinear(*(P[t] for t in T)):
            line = frozenset(k for k in range(n) if conf.collinear(P[T[0]], P[T[1]], P[k]))
            lines[line] = True
    most = max((len(L) for L in lines), default=2)
    k = n - 2 * budget
    for T in noncol:
        consider(fd._locate(conf, T, zero, zero, zero))
    anchors = fd.anchor_triples(conf, k) if most < k else noncol
    for T in anchors:
        for r, s, t in product(vals, repeat=3):
            y = fd._locate(conf, T, r, s, t)
            if conf.d2(y, P[T[0]]) == r:
                consider(y)
    if most >= k:
        assert most == k, "more collinear points than the argument allows"
        for S in combinations(range(n), k):
            if any(T in noncol for T in combinations(S, 3)):
                continue
            rest = [i for i in range(n) if i not in S]
            for m in _matchings(rest):
                rows = [_bisector(conf, P[a], P[b]) for a, b in m]
                meets = [_meet(r1, r2) for r1, r2 in combinations(rows, 2)]
                for y in meets:
                    consider(y)
                if all(y is None for y in meets):
                    for y in _on_common_bisector(conf, rows, [P[i] for i in S], vals):
                        consider(y)
    return good


def _on_common_bisector(conf, rows, S, vals):
    """Points y on the common line of the (parallel) rows, at old squared distances from all
    points of S.  Returns points in the field of conf; raises if irrational points occur."""
    a0, b0, c0 = rows[0]
    if not all(a0 * r[1] == b0 * r[0] and a0 * r[2] == c0 * r[0] and b0 * r[2] == c0 * r[1]
               for r in rows[1:]):
        return []                                       # parallel and distinct: no point
    if a0 != 0:
        y0, e = (c0 / a0, 0 * c0), (-b0 / a0, 1 + 0 * b0)
    else:
        y0, e = (0 * c0, c0 / b0), (1 + 0 * a0, 0 * a0)
    quads = []                                          # |y0 + t e - u|^2 = al t^2 + be t + ga
    for u in S:
        w = (y0[0] - u[0], y0[1] - u[1])
        quads.append((conf.dot(e, e), 2 * conf.dot(w, e), conf.dot(w, w)))
    out = []
    for targets in product(vals, repeat=len(S)):
        (al, b1, g1) = quads[0]
        g1 = g1 - targets[0]
        lin = [(q[1] - b1, q[2] - targets[j] - g1) for j, q in enumerate(quads) if j > 0]
        if all(p == 0 and q == 0 for p, q in lin):
            raise ValueError("a whole quadratic of solutions: irrational points possible")
        if any(p == 0 and q != 0 for p, q in lin):
            continue
        ts = {-q / p for p, q in lin if p != 0}
        if len(ts) == 1:
            t = ts.pop()
            if al * t * t + b1 * t + g1 == 0:
                out.append((y0[0] + t * e[0], y0[1] + t * e[1]))
    return out


def completions(conf, budget, need=NEED):
    """All sets of `need` candidate points that can be added together within the budget."""
    old = conf.distances()
    cand = list(candidate_points(conf, budget).items())
    found = []

    def search(start, chosen, new):
        if len(chosen) == need:
            found.append(tuple(chosen))
            return
        if len(cand) - start < need - len(chosen):
            return
        for j in range(start, len(cand)):
            y, new_y = cand[j]
            total = new | new_y
            if len(total) > budget:
                continue
            for c in chosen:
                d = conf.d2(y, cand[c][0])
                if d not in old:
                    total = total | {d}
                    if len(total) > budget:
                        break
            if len(total) <= budget:
                search(j + 1, chosen + [j], total)

    search(0, [], frozenset())
    return [c[0] for c in cand], found


def extend(sol):
    """The shapes of all fifteen-point sets containing the nine-point solution `sol`."""
    roots = sol.embeddings
    if not roots:
        return {}, 0
    K = fd.place(sol.pattern, sol.lengths)
    cands, comps = completions(K, 7 - len(K.distances()))
    shapes = {}
    for root in roots:
        for c in comps:
            shapes.setdefault(shape_key(K.with_points(cands[j] for j in c), root, sol.field), 0)
    return shapes, len(roots)


# ---- the nine-point patterns ----------------------------------------------------------------

def load_L9():
    return fd.load_patterns(os.path.join(HERE, "..", "certificates", "L9.json"))


def solve_all(patterns, path, resolve=False):
    if os.path.exists(path) and not resolve:
        data = json.load(open(path))["by_pattern"]
    else:
        data = []
        for p in patterns:
            try:
                data.append({"solutions": [s.to_json() for s in fd.solve_lengths(p)]})
            except ValueError:
                data.append({"flexible": True})
        json.dump({"patterns": len(data), "by_pattern": data}, open(path, "w"),
                  separators=(",", ":"))
    fields, out = {}, []
    for p, d in zip(patterns, data):
        out.append(None if d.get("flexible") else
                   [fd.LengthSolution.from_json(p, x, fields) for x in d["solutions"]])
    return out


_WORK = {}


def _extend_index(args):
    i, k = args
    shapes, cores = extend(_WORK["solutions"][i][k])
    return list(shapes), cores


# ---- the two flexible families --------------------------------------------------------------

QS, s = sp.field("s", sp.QQ)


def family(pattern, index):
    """The squared lengths of a flexible pattern as elements of Q(s), after checking that the
    parametrisation covers every solution of its planarity equations."""
    xs, values = fd.length_symbols(pattern)
    sat = fd._saturated(pattern, "principal")
    xsym, w, S = sat
    basis = {sp.Poly(g, *xs).monic() for g in sp.groebner(list(S.exprs), w, *xs,
                                                           order="lex").exprs if not g.has(w)}
    X = dict(zip(map(str, xs), xs))
    if index == 124:    # x1 = (x4+1)/2, x2 = x1 - x3, 12x3^2 - 6x3x4 - 6x3 + x4^2 - 2x4 + 1 = 0
        x1, x2, x3, x4 = (X[f"x{j}"] for j in range(1, 5))
        expect = [2 * x1 - x4 - 1, 2 * x2 + 2 * x3 - x4 - 1,
                  12 * x3**2 - 6 * x3 * x4 - 6 * x3 + x4**2 - 2 * x4 + 1]
        # the conic 3 x3 (4 - x3) = w^2, w = x4 - 3 x3 - 1, through (x3, w) = (0, 0): w = s x3
        v3 = 12 / (s**2 + 3)
        v4 = 3 * v3 + 1 + s * v3
        value = {x1: (v4 + 1) / 2, x2: (v4 + 1) / 2 - v3, x3: v3, x4: v4}
        note = "every solution except x3 = 0 (a zero length), the limit s -> infinity"
    elif index == 520:  # x1 + x2 = (x5+1)/2, (4x2 - x5 - 1)^2 = 4x5, x3, x4 linear in x5
        x1, x2, x3, x4, x5 = (X[f"x{j}"] for j in range(1, 6))
        expect = [2 * x1 + 2 * x2 - x5 - 1,
                  16 * x2**2 - 8 * x2 * x5 - 8 * x2 + x5**2 - 2 * x5 + 1,
                  4 * x3 - x5 - 3, 4 * x4 - 3 * x5 - 1]
        value = {x5: s**2, x2: (s + 1)**2 / 4, x3: (s**2 + 3) / 4, x4: (3 * s**2 + 1) / 4}
        value[x1] = (s**2 + 1) / 2 - value[x2]
        note = "every solution (s = +-sqrt(x5) covers both branches)"
    else:
        raise ValueError(index)
    assert basis == {sp.Poly(e, *xs).monic() for e in expect}, basis
    lengths = [QS(1) if v == 1 else value[v] for v in values]
    as_expr = lambda v: sp.sympify(str(v), locals={"s": sp.Symbol("s")})
    for e in expect:                                    # the parametrisation solves the basis
        assert sp.simplify(e.subs({x: as_expr(value[x]) for x in xs})) == 0
    return lengths, note


def _numerator(f):
    e = f.numer.as_expr()
    t = sp.Symbol("s")
    poly = sp.Poly(e, t) if e.free_symbols else sp.Poly(e, t, domain="QQ")
    return poly.clear_denoms()[1].primitive()[1] if not poly.is_zero else None


def special_values(pattern, lengths):
    """Irreducible polynomials (sympy Polys in s) whose roots are the only parameters s0 at
    which K(s0) can be part of a fifteen-point set with at most seven distances."""
    K = fd.place(pattern, lengths)
    P = K.points
    n = len(P)
    old = K.distances()
    vals = sorted(old, key=str)
    b = 7 - len(old)
    polys = [_numerator(K.A * K.B - K.C * K.C)]              # the basis of `place` degenerates
    noncol = [T for T in combinations(range(n), 3) if not K.collinear(*(P[t] for t in T))]
    for T in noncol:                                           # a triple becomes collinear
        u, v, w = (P[t] for t in T)
        polys.append(_numerator((v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])))
    zero = vals[0] - vals[0]
    ys = [fd._locate(K, T, zero, zero, zero) for T in noncol]
    for T in fd.anchor_triples(K, n - 2 * b):
        for r, t, u in product(vals, repeat=3):
            y = fd._locate(K, T, r, t, u)
            e = K.d2(y, P[T[0]]) - r
            if e == 0:
                ys.append(y)
            else:
                polys.append(_numerator(e))                    # this point exists only there
    generic = 0
    for y in set(ys):
        if y in P:
            continue
        ds = [K.d2(y, q) for q in P]
        newv = sorted({d for d in ds if d not in old}, key=str)
        if len(newv) <= b:
            generic += 1
            continue
        # valid at s0 only if, among any b+1 of these values, two agree or one becomes old
        g = None
        for group in combinations(newv, b + 1):
            f = QS(1)
            for v1, v2 in combinations(group, 2):
                f *= v1 - v2
            for v1 in group:
                for x in vals:
                    f *= v1 - x
            q = _numerator(f)
            g = q if g is None else g.gcd(q)
            if g.degree() == 0:
                break
        polys.append(g)
    assert generic == 0, "a candidate point that works for every s"
    factors = set()
    for q in polys:
        if q is None or q.degree() <= 0:
            continue
        for f, _ in q.factor_list()[1]:
            if f.degree() > 0 and f.count_roots() > 0:
                factors.add(f.monic())
    return factors


def at_value(pattern, lengths, f):
    """The nine-point set K(s0), s0 a root of f, as a LengthSolution (or None)."""
    c = [Fraction(int(x.p), int(x.q)) for x in reversed(f.all_coeffs())]
    if len(c) == 2:
        s0, one = -c[0], Fraction(1)
    else:
        F = fd.NumberField(c)
        s0, one = F.t, F.one

    def ev(g):
        e = sp.Poly(g.as_expr(), sp.Symbol("s")) if g.as_expr().free_symbols else None
        cs = e.all_coeffs() if e else [g.as_expr()]
        v = one * 0
        for x in cs:
            x = sp.Rational(x)
            v = v * s0 + one * Fraction(int(x.p), int(x.q))
        return v

    try:
        L = [ev(x.numer) / ev(x.denom) for x in lengths]
    except ZeroDivisionError:
        return None
    if any(x == 0 for x in L) or len(set(L)) < len(L):
        return None
    return fd.LengthSolution(pattern, L)


# ---- main -----------------------------------------------------------------------------------

def main():
    args = sys.argv[1:]
    jobs = int(args[args.index("--jobs") + 1]) if "--jobs" in args else 1
    t0 = time.time()
    patterns = load_L9()
    solutions = solve_all(patterns, os.path.join(HERE, "solutions9.json"), "--resolve" in args)
    flexible = [i for i, x in enumerate(solutions) if x is None]
    print(f"{len(patterns)} nine-point patterns: {len(patterns) - len(flexible)} with finitely "
          f"many solutions, flexible: {flexible}  ({time.time() - t0:.0f} s)", flush=True)

    shapes, cores = set(), 0
    _WORK["solutions"] = solutions
    tasks = [(i, k) for i, x in enumerate(solutions) if x for k in range(len(x))]
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
    with ctx.Pool(jobs) as pool:
        for keys, c in pool.imap_unordered(_extend_index, tasks, chunksize=8):
            shapes.update(keys)
            cores += c
    print(f"rigid nine-point sets: {cores} realisations extended; {len(shapes)} fifteen-point "
          f"shapes  ({time.time() - t0:.0f} s)", flush=True)

    for i in flexible:
        lengths, note = family(patterns[i], i)
        factors = special_values(patterns[i], lengths)
        members = 0
        for f in factors:
            sol = at_value(patterns[i], lengths, f)
            if sol is None or not sol.embeddings:
                continue
            members += 1
            keys, _ = extend(sol)
            shapes.update(keys)
        print(f"family {i} (parametrisation covers {note}): {len(factors)} special values, "
              f"{members} of them genuine; {len(shapes)} shapes so far  "
              f"({time.time() - t0:.0f} s)", flush=True)

    target = tuple(("Q", Fraction(v)) for v in PALETTE)
    assert all(k[0] == target for k in shapes), "a fifteen-point set with another palette"
    found = fd.sets_with_distances(PALETTE, 15)
    assert {k[1] for k in shapes} == {fd.canonical(key, named=True) for key in found}
    result = []
    for key, conf in found.items():
        result.append({"A": str(conf.A), "B": str(conf.B), "C": str(conf.C),
                       "points": [[str(x), str(y)] for x, y in conf.points],
                       "census": {str(k): v for k, v in sorted(Counter(
                           conf.d2(p, q) for p, q in combinations(conf.points, 2)).items())}})
    json.dump({"sets": result, "plus": "the regular 15-gon"},
              open(os.path.join(HERE, "result.json"), "w"), indent=1)
    print(f"all {len(shapes)} fifteen-point sets have squared distances proportional to "
          f"{PALETTE}, and they are exactly the {len(found)} sets with these squared distances "
          f"(result.json); with the regular 15-gon, {len(found) + 1} sets  "
          f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
