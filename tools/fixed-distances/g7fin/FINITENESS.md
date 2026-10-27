# Finiteness of maximum d-distance sets: the plan, and the state of the d = 7 case

Conjecture F(d): the set M_d of d-distance sets with g(d) points, modulo similarity, is finite.
It follows from the Erdős–Fishburn lattice conjecture (lattice pieces are countable; a positive-dimensional
semi-algebraic set is not), and it is known for d ≤ 6 from the classifications. This note sets up the argument
that proves F(7) modulo one computational input, and records what is done.

## 1. Preliminaries

**Semi-algebraicity.** M_d is a finite union of colouring strata (one for each realizable colouring of K_n with
d colours, n = g(d)); each stratum is a real algebraic set minus a closed set (the inequations). So M_d is
semi-algebraic and "finite" means "real dimension 0".

**Constant-colouring lemma.** If M_d has positive real dimension, some stratum has positive dimension and
therefore contains a smooth arc X(t), t ∈ (0,1), on which the colouring, hence the diameter class and the
diameter graph G_D, are constant. Along such an arc at least one ratio of distance values varies (K_n with all
distances fixed is rigid).

**Pinning lemma.** Let X(t) be such an arc and A ⊂ X a labelled subset whose shape is constant along the arc
(A is *rigid*). Call a colour class *fixed* if it occurs inside A (its value is then a constant multiple of A's
scale) and *free* otherwise. A point x ∈ X \ A is determined up to finitely many positions by A unless all of
its distances to A lie in free classes, each of which contains at most two points of A, all pairs sharing one
perpendicular bisector. (Two fixed-class distances: two circles. One fixed class plus a free pair: circle and
line. A free class on three non-collinear points: x is their circumcentre. Two free pairs with different
bisectors: a point.) Hence x is pinned whenever |A| > 2·(number of free classes). With d = 7 classes and |A| ≥ 6
(so A has k ≥ 3 classes, free ≤ 4) every point outside A is pinned, and the whole of X(t) is constant: **a
family of d-distance sets cannot contain a rigid subset of ≥ 6 points when d ≤ 7.**

**Core lemma (Wei).** Let m = |X_D| be the number of diameter endpoints. G_D has no even cycle and at most one
cycle; if G_D ≠ C_m its independence number is ≥ ⌈m/2⌉, and for C_m (m odd) it is (m−1)/2. An independent set
I together with X \ X_D is a set A with no diameter pair, so A has at most d − 1 distances and
|A| ≥ n − m + α(G_D).

So: **F(d) holds as soon as, for every m that the diameter-graph analysis allows, every (n − m + α)-point set
with ≤ d − 1 distances is rigid** (zero-dimensional in its colouring stratum). This is the induction: finiteness
at the top of column d reduces to rigidity of near-maximum sets in column d − 1.

**Open lemma (not needed here).** g(d−1) < g(d). It is what Strategy B (degeneration ends) would need; it is
true for d ≤ 7 and I know no general proof.

## 2. The case d = 7, n = 16

Let X(t) be a family of 16-point 7-distance sets, m = |X_D|.

| m | α(G_D) | core size | what is needed | status |
|---|---|---|---|---|
| ≤ 10 | ≥ ⌈m/2⌉ | ≥ 11 | 11-point ≤6-distance sets rigid | **done**: 62 patterns, finitely many realizations (Amit–Astra Lemma 6; re-checked here, `pat_rigid 11`) |
| 11, G_D ≠ C₁₁ | ≥ 6 | ≥ 11 | same | done |
| 11, G_D = C₁₁ | 5 | 10 | 10-point ≤6-distance sets rigid | 225 patterns known and certified rigid; SÖ count 228 |
| 12 | ≥ 6 | ≥ 10 | same | same |
| 13, G_D ≠ C₁₃ | ≥ 7 | ≥ 10 | same | same |
| 13, G_D = C₁₃ | 6 | 9 | not needed: X_D is a 13-point set with diameter graph C₁₃ and ≤ 7 distances, so by Amit–Astra Lemma 4 it is the regular 13-gon, a rigid 13-point subset | **done** |
| 14 | — | — | X_D is a convex 14-gon with ≤ 7 distances, hence exactly 7, hence R₁₄ or R₁₅ minus a vertex (Erdős–Fishburn): rigid with all 7 classes fixed | **done** |
| 15 | — | — | convex 15-gon with ≤ 7 distances: R₁₅ | **done** |
| 16 | — | — | convex 16-gon needs ≥ 8 distances | impossible |

Therefore: **F(7) ⇐ every 10-point planar set with at most six distances is rigid.**

## 3. The ten-point input

Szöllősi–Östergård, Table 9 (arXiv:1804.06040): the orderly generation of candidate Gram patterns in the plane with
at most six distinct distances gives 228 candidates at order 10 (and 62, 15, 3 at orders 11, 12, 13; the
order-11 count is the one used in the g(7) proof). "Candidate" means: passes the complex rank condition (5) with
distinct nonzero values, hereditarily. Every pattern of a real 10-point set with ≤ 6 distances is a candidate.

Known realizable patterns at order 10 (`g7fin/pat10.py`): 10-point subsets of R₁₃, R₁₂ ∪ O, R₁₁ ∪ O, R₁₀ ∪ O and of
the eleven 12-point lattice cores, the applet's (10,5) and (10,6) catalogs, an independent enumeration of all
10-point triangular-lattice subsets with ≤ 6 distances (199 similarity classes; stable when the scope is widened
to minimum norm ≤ 7 and diameter ≤ 6 × minimum distance), and the 24 six-or-fewer-distance endpoints of the
two-concentric-pentagon families (`g7fin/pent_ends.py`). Up to equivalence (vertex and colour permutations,
confirmed with VF2 inside hash buckets): **225 patterns** (19 with five distances, 206 with six); at order 11 the
same procedure gives exactly the 62 of the paper, and `pat_rigid 11` certifies all 62 zero-dimensional over ℚ
(an independent check of the rigidity part of Amit–Astra's Lemma 6).

Rigidity of each known pattern (`g7fin/pat_rigid.py`): the squared-distance equations with the scale fixed and a
Rabinowitsch variable for "values distinct and nonzero"; msolve's grevlex Gröbner basis modulo a prime screens
for zero-dimensionality (a pure power of every variable among the leading monomials), and Singular's `dim(std(I))`
over ℚ certifies it. Result: **all 225 known patterns are zero-dimensional** (`pat_rigid_10.json`, 5–70 s each).

**Gap.** 228 − 225 = 3 candidates are unaccounted for: either complex-realizable only, or real 10-point sets not
in any known family of examples (not lattice, not circular, not in the catalogs). Closing it needs either the
candidate list itself (not published; the authors could supply it) or a re-run of the orderly generation
(levels 4–10; the expensive level is n = 6 with 537 230 raw candidates). Until then F(7) is proved only for
families whose cores avoid those three patterns, i.e. the theorem is: *a one-parameter family of 16-point
seven-distance sets would have to contain a 10-point core realizing one of three specific unknown Gram patterns
(six distances, not a subset of any known 11-point set, not a lattice piece, not circular, not a pentagon
configuration).*

**A sharper reformulation of what a family must look like.** By the pinning lemma, any 11-point subset of X(t)
missing one of the seven colour classes is a ≤6-distance set, hence rigid, hence X(t) is rigid. So in a family
every colour class must meet every 11-subset: every class graph has vertex-cover number ≥ 6 (matching number ≥ 3).
Likewise, once the 10-point classification is complete, every class pair must meet every 10-subset.

## 4. Strategy B is dead in its naive form

Test (`g7fin/extend_ends.py 4`): along the 25 one-parameter families of (5,4), near their collision ends
(min distance / diameter = 0.02 … 0.2) 28 % of sampled configurations admit a one-point extension keeping the
palette, the same fraction as at interior points (24 %) — extendability is a property of the whole family
(algebraic), not of the end. So nearing a collision gives maximality nothing to bite on. The rigid sets of the
cell are extendable 40 % of the time. The degeneration structure might still force algebraicity, but the
"unfold near the end" idea does not work.
