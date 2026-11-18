# Fourteen points, seven distances

Up to similarity there are exactly **84** sets of 14 points in the plane with at most seven
distinct distances:

* the regular 14-gon;
* the regular 15-gon minus a vertex;
* the regular 13-gon with its centre;
* 81 pieces of the triangular lattice. These are all the 14-point sets whose squared distances
  are proportional to 1, 3, 4, 7, 9, 12, 13 (listed in `result.json`). Four of them are the exceptional
  fifteen-point set F of `../fifteen` minus a point.

No 14-point set with seven distances belongs to a continuous family. All 84 are rigid.

## The proof

**Reduction.** Let X have 14 points and at most seven distances. Let h be the number of vertices
of its convex hull, D its diameter, and G_D the graph of pairs at distance D, with m vertices.
The lemma numbers are those of `../paper/seven-distances.pdf`.

* If h ≥ 14, Lemma 10 gives X = R₁₄ or R₁₅ minus a vertex.
* Otherwise m ≤ h ≤ 13. A set A made of an independent set of G_D together with the points
  outside G_D avoids the distance D, so it has at most six distances, and |A| = 14 − m + (size of
  the independent set).
  * If G_D is not a single cycle, Lemma 8(3) gives an independent set with at least ⌈m/2⌉
    vertices, so |A| ≥ 14 − ⌊m/2⌋ ≥ 8.
  * If G_D is a single (odd) cycle with m ≤ 11, every other vertex gives |A| = 14 − (m+1)/2 ≥ 8.
  * If G_D is a 13-cycle, its vertices form a regular 13-gon (Lemma 11), which has six distances.

So every 14-point set with at most seven distances, other than R₁₄ and R₁₅ minus a vertex,
contains eight points with at most six distances. These have one of the 1565 patterns of
`../certificates/L8.json` (the certified generation of the paper, §2–3).

**Eight points.** Of the 1565 patterns, 1551 have finitely many solutions, which give 2171 real
point sets. The other 14 have a curve of solutions. In each case the lexicographic Gröbner basis
consists of linear equations and one conic in two of the squared lengths. The conic is
parametrised by the lines through a rational point on it. The basis is checked against the
parametrisation, and the points the parametrisation misses are treated separately.

**Adding six points.** Let K be an eight-point set with |Δ(K)| distances. The budget for new
distances is b = 7 − |Δ(K)| ∈ {1, 2, 3}. Let y be a point that can be added, S the set of points
of K at old distances from y, and r = 8 − |S|.

1. If one value occurs three times among the eight distances from y, then y is the circumcentre of
   three points.
2. Otherwise every value occurs at most twice. Then r ≤ 2b, and at least r − b of the new values
   occur exactly twice.
3. If S contains three non-collinear points, y is fixed by its old distances to an anchor triple,
   one of a short list contained in every set of 8 − 2b non-collinear points.
4. Otherwise S is collinear and has at least two points. For each new value taken twice, y lies on
   the perpendicular bisector of the two points.
   * Two such bisectors that meet fix y.
   * If all of them coincide in one line L, intersect L with the line |y−u|² − |y−v|² = ρ − σ
     for u, v ∈ S. If that line is parallel to L, then y = f ± h·ν with h² = ρ − |f−u|²
     (ν normal to uv).
   * If no new value is taken twice, S has at least 8 − b collinear points, and y = f ± h·ν in the
     same way.

The points f ± h·ν can be irrational over the field of K. If h is in the field (decided
exactly: by SymPy factoring x² − h² over the number field, with the square root then checked),
they join the ordinary candidates. Otherwise their distances a + h·b are compared exactly as pairs
(a, b). The few that survive are combined with all other candidates in a search that uses
certified interval arithmetic (mpmath.iv). A combination is discarded only when the intervals
force more than b distinct new values. No combination survives.

The exact candidates are combined by a depth-first search, as in `../fifteen`.

**The families.** As in `../fifteen`, everything is done over Q(s). Some candidate points of these
families work for every s. For every six of them, the new values they produce are compared with s
symbolic. Every other candidate, every degeneration of the construction, and every six-point choice
among the candidates valid for all s gives a polynomial in s. A fourteen-point set can contain K(s₀)
only if s₀ is a root of one of these polynomials. At each real root, K(s₀) is extended exactly in
the number field Q(s₀).

**Result.** The extensions give 82 sets: 81 with squared distances proportional to
1, 3, 4, 7, 9, 12, 13, which are exactly the sets found by `fewdist.sets_with_distances` for that
palette, and the regular 13-gon with its centre.

## Running it

    python3 fourteen.py --jobs 2          # 68 minutes on two cores (run.log)
    python3 fourteen.py --resolve         # also re-solves the 1565 patterns (+ 20 minutes)

| file | what it is |
|---|---|
| `fourteen.py` | the program |
| `solutions8.json` | the solutions of the 1565 eight-point patterns (the basis for the 14 flexible ones) |
| `result.json` | the 81 lattice sets (Gram basis A, B, C, exact points, distance census) |
| `run.log` | the output of a full run |
