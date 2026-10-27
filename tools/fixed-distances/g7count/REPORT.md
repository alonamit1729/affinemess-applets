# How many 16-point planar sets have exactly seven distinct distances?

**Answer: exactly two, up to similarity** — the two triangular-lattice pieces of Erdős–Fishburn (§5).
This is the complete classification of the extremal sets for g(7) = 16.  It rests on the finiteness argument of
`g7fin/FINITENESS.md` (re-verified in §1: every 16-point 7-distance set contains a rigid 10-point core with ≤ 6
distances, or is a regular 14/15-gon plus 1–2 points), on msolve's exact real-root isolation of all real
realizations of the 228 ten-point patterns, and on an exhaustive extension search in 100-digit arithmetic.
Everything is in `scratchpad/g7count/`; nothing outside it was modified.

## 1. Where every 16-point 7-distance set lives (verified)

Let X be a 16-point set with exactly 7 distances, X_D its m diameter endpoints, G_D the diameter graph on X_D.
Elementary facts (proofs re-derived here, not just cited):

* every diameter endpoint is an extreme point of X (a convex combination of other points is at distance < D from
  the partner), so no diameter segment contains a third point of X, and two diameter segments are never collinear;
* any two diameter segments meet (share an endpoint or cross);
* **no even cycle**: if c₁…c₂ᵣ were a cycle, the path c₃…c₂ᵣ (2r−3 edges, each crossing the segment c₁c₂ exactly
  once, no vertex on the line c₁c₂) puts c₃ and c₂ᵣ on opposite sides of the line c₁c₂, so the edges c₂c₃ and c₁c₂ᵣ
  cannot meet;
* **cycle lemma**: an edge pq vertex-disjoint from a cycle C must cross each of the |C| (odd) edges of C once, but
  p, q are extreme points of X, hence outside conv(C), so pq crosses the closed curve C an even number of times —
  contradiction.  Hence every edge meets the cycle; two distinct cycles would give an edge disjoint from one of them
  or an even cycle, so G_D has at most one cycle;
* therefore G_D is a forest (α ≥ ⌈m/2⌉) or one odd cycle with pendant leaves (α ≥ ⌈m/2⌉ if there is at least one
  leaf; α = (m−1)/2 for G_D = C_m).  This is Wei's bound, confirmed.

A core A = I ∪ (X ∖ X_D), I independent in G_D, contains no diameter pair, so |Δ(A)| ≤ 6 and |A| = 16 − m + α:

| m | α(G_D) | min |A| | consequence |
|---|---|---|---|
| ≤ 10 | ≥ (m−1)/2 | ≥ 11 | X contains an 11-point ≤6-distance set (fast pass), hence also a 10-point one |
| 11 | 5 (C₁₁) or ≥ 6 | 10 | X contains a 10-point ≤6-distance set |
| 12 | ≥ 6 | 10 | same |
| 13, G_D ≠ C₁₃ | ≥ 7 | 10 | same |
| 13, G_D = C₁₃ | 6 | 9 | X_D is a 13-point set with diameter graph C₁₃ and ≤ 7 distances, hence R₁₃ (Amit–Astra Lemma 4); R₁₃ has 6 distances, so any 10 of its points form a 10-point ≤6-distance set |
| 14 | — | — | X_D is a convex 14-gon with exactly 7 distances: R₁₄ or R₁₅ − v (Erdős–Fishburn); X = X_D ∪ {q₁, q₂} and every distance from qᵢ to X_D is a chord length of X_D (§4) |
| 15 | — | — | X_D = R₁₅, X = R₁₅ ∪ {q}, same condition (§4) |
| 16 | — | — | a convex 16-gon has ≥ 8 distances: impossible |

Every 10-point planar set with ≤ 6 distances has its colour pattern among the 228 of
`g7fin/gen/levels/L10.json` (regenerated Szöllősi–Östergård list), all rigid.  So apart from the polygon cases,
X = A ∪ Q with A one of the finitely many real realizations of one of the 228 patterns, and the six points of Q
have all their distances (to A and among themselves) in Δ(A) ∪ N, |N| = 7 − |Δ(A)| ∈ {1, 2}.

**Candidates for a point q ∈ Q** (|A| = 10, |N| ≤ 2): if q has ≥ 2 distances to A in Δ(A) it is an intersection of
two circles centred at points of A with radii in Δ(A) (≤ 2 points per choice); otherwise ≥ 9 of its 10 distances
lie in N, one value occurs ≥ 5 times, and q is the circumcentre of a non-collinear triple of A.  Both families are
finite and enumerated completely (`extend.py: candidates`, ~3100 circle points + 120 circumcentres per core);
a depth-first search over the surviving candidates (those with ≤ budget new values to A) enumerates every
6-subset whose pairwise distances stay in Δ(A) ∪ N with the running set N never exceeding the budget
(`extend.py: extensions`).  The same routine with 5 extra points serves the 11-point cores.

## 2. Fast pass: 16-point supersets of the 68 eleven-point ≤6-distance sets (`fast11.py`, 18 s)

The 62 eleven-point patterns (`g7fin/patterns_11.json`, Amit–Astra Lemma 6) are realized by 68 point sets:
52 one-point deletions of the eleven 12-point lattice cores (`g7/cores.json`, dedup by congruence) and 16 circular
sets (R₁₃ − {0,k}, k = 1..6; R₁₂ − v; R₁₂ ∪ O − {0,k}, k = 1..6; R₁₀ ∪ O; R₁₁; R₁₁ − v ∪ O); `cores11.py` checks that
these 68 sets realize all 62 patterns (the six R₁₃ − {0,k} share one pattern because multiplication by a unit of ℤ₁₃
permutes the chord classes).  For each core every 16-point superset with ≤ 7 distances was enumerated, with no
assumption on which value is the diameter (unlike the paper's Lemma 7, which only bounds |X|).

Result (`fast11_result.json`): ≤ 14 surviving candidates per core, 0–9 labelled 5-extensions per core, and
**exactly 2 similarity classes of 16-point sets, both with exactly 7 distances, congruent to the two known sets**
(`r12/r12_subsets_7.json["16,7"]`).  By the table this already classifies all X with m ≤ 10 and m = 11, G_D ≠ C₁₁,
given Lemma 6; the complete pass makes Lemma 6 unnecessary.

## 3. Complete pass: all real realizations of the 228 ten-point patterns

**Realizations** (`realize_all.py`, `msolve_util.py`): msolve real root isolation (`-p 400`, isolating intervals of
width 2⁻⁴⁰⁰) on the coordinate system of `g7fin/gen/realize.py` (p₀ = (0,0), p₁ = (x₁,0), value 1 on the colour of
edge 01, free values s_c, Rabinowitsch variable t).  218 of the 228 patterns finished this way (total msolve time
5554 s, median ≈ 2 s, max 658 s; 1077 labelled real solutions, 4 per pattern typically = reflection × sign of x₁).
Ten patterns (213, 215, 217, 219, 221, 223, 224, 225, 226, 227; the ones msolve handles with an added generic
linear form) were too heavy for this machine (pattern 224 was OOM-killed at 3.8 GB, pattern 213 ran > 40 min); for
them, and as an independent cross-check for all 228, the **Gram route** (`gram_util.py`, `gram_all.py`) was used:
msolve real roots of Szöllősi–Östergård's own system (all 3×3 minors of the (n−1)×(n−1) Gram matrix in the ≤ 5 free
values plus t; 0.3–1.0 s per pattern, 139 s in all), then for each real solution the Gram matrix's eigenvalues at
100 digits (seven ≈ 10⁻¹⁰⁰, two clearly positive: PSD of rank 2, i.e. a genuine Euclidean realization — none was
indefinite) and coordinates by eigendecomposition.  On the 218 patterns done both ways the two routes give
congruent realizations, class by class (`crosscheck` in the session log); every reconstructed set reproduces its
pattern exactly (asserted).

**Count of 10-point sets**: the 228 patterns have **257 real realizations up to similarity** (214 patterns with one,
9 with two, 2 with four, 1 with five, 2 with six — the multi-realization patterns are those with colour
automorphisms, e.g. R₁₃-type ones).  Cross-checks (`crosscheck.py`): all 210 ten-point sets of the R₁₂+R₁₂ catalog
with ≤ 6 distances (`r12/r12_subsets_7.json["10,5"]`, `["10,6"]`) and all 236 similarity classes of 10-subsets of
the 68 eleven-point cores are among the 257; no two patterns share a realization.

**Extension** (`complete10.py` → `ext10/pat_<i>.json`, `collect.py` → `complete_result.json`; 65 s in all): every
one of the 257 cores was extended by 6 points with budget 7 − |Δ(A)|.  ≤ 16 surviving candidates per core, ≤ 1084
DFS nodes, 826 labelled 6-extensions from 185 cores, every extension has exactly 7 distances, and up to similarity
they are **exactly 2 sets**: E₁ (694 labelled hits from 184 cores) and E₂ (132 hits from 89 cores), congruent to
the two known sets.  Consistency: E₁ has exactly 184 and E₂ exactly 89 similarity classes of 10-subsets with
≤ 6 distances, i.e. every admissible core of each set produced it and nothing else did.

## 4. The polygon cases (`polygons.py`, 1.4 s)

For X_D ∈ {R₁₄, R₁₅ − v, R₁₅} all seven distances of X occur inside X_D, so an extra point q has all 14 or 15 of its
distances among the chord lengths of X_D; such a q has ≥ 2 distances in Δ(X_D) and is a circle–circle intersection,
so `candidates(X_D, 0)` lists all of them: **R₁₄: none; R₁₅: none; R₁₅ − v: exactly one, namely v itself** (which
gives R₁₅, 15 points).  (The axis argument predicts the same: a circle about q meets the circumcircle in ≤ 2 points,
so 14 vertices need all 7 values twice and q on a symmetry axis; 15 vertices force q = centre, and the circumradius
is not a chord length of R₁₄ or R₁₅.)  Hence no 16-point 7-distance set has m ≥ 14.

## 5. The two sets and their certification (`certify16.py`, `certify16.json`)

Triangular-lattice coordinates (a, b) ↦ (a + b/2, b√3/2), Q(a,b) = a² + ab + b²; both sets have squared-distance
spectrum {1, 3, 4, 7, 9, 12, 13} (minimum distance 1, diameter √13), both are subsets of R₁₂ + R₁₂ (the two entries
of `r12_subsets_7.json["16,7"]`), and neither contains a 13-point 6-distance set:

* **E₁** (symmetry group of order 2): b = −2: a ∈ {0, 1}; b = −1: {−1, 0, 1, 2}; b = 0: {−2, −1, 0, 1};
  b = 1: {−2, −1, 0, 1}; b = 2: {−2, −1}.
* **E₂** (symmetry group of order 6): b = −2: a ∈ {1, 2}; b = −1: {−1, 0, 1, 2}; b = 0: {−2, −1, 0, 1};
  b = 1: {−2, −1, 0, 1}; b = 2: {−1, 0}.

Both are the 12-point block of three rows of four (b = −1, 0, 1) with a 2-point row above and below; the two short
rows sit on the same side (E₁) or on opposite sides (E₂).

Certification: msolve real root isolation (precision 120) on each 16-point colour pattern (36 variables, 121
equations) is zero-dimensional with **exactly 4 real labelled solutions, all congruent to the set found** (4 =
reflection × sign of x₁), so each pattern is rigid with a unique realization (62 s and 51 s).

## 6. Soundness assumptions

1. Completeness of the 228-pattern list and rigidity of all 228 (`g7fin`, regenerated orderly search agreeing with
   Szöllősi–Östergård at every order); msolve's real root isolation (RUR + exact isolation) then returns every real
   solution.  As in the rest of the project, msolve's rational parametrization is multi-modular with probabilistic
   verification.  The Gram route adds the (exact) equivalence rank-2 symmetric ⇔ V^T V over ℂ and the PSD test at
   100 digits (eigenvalues separated from 0 by ≥ 0.7 in all cases).
2. Extension arithmetic: mpmath at 100 decimal digits (gmpy backend); two squared distances of a set scaled to
   minimum distance 1 are equal iff they differ by < 10⁻⁴⁰; coordinates are algebraic numbers of small degree
   known to ~10⁻⁹⁵, so no true coincidence is missed, and a spurious one could only add candidate sets, each then
   certified exactly (§5).  Candidate points are deduplicated at 10⁻⁴⁰; tangent circles (|h²| < 10⁻⁴⁰) give their
   contact point.
3. Cited theorems: Erdős–Fishburn (convex n-gons with ⌊n/2⌋ distances are regular n-gons or regular (n+1)-gons minus
   a vertex; m = 14, 15) and Amit–Astra 2026 Lemma 4 (13 points with diameter graph C₁₃ and ≤ 7 distances form R₁₃).
   Lemma 6 (11-point classification) is used only by the fast pass.

## 7. Files and timings

`extend.py` (candidates, DFS, congruence, 100-digit arithmetic), `cores11.py` (68 eleven-point sets),
`fast11.py` → `fast11_result.json` (18 s), `msolve_util.py`, `realize_all.py` → `real10/pat_<i>.json` (coordinate
route, 5554 s; the ten heavy patterns hold coordinates reconstructed from the Gram route, field `source`),
`gram_util.py`, `gram_all.py` → `real10g/pat_<i>.json` (Gram route, all 228, 139 s), `complete10.py` →
`ext10/pat_<i>.json` (65 s), `collect.py` → `complete_result.json`, `crosscheck.py`, `polygons.py` (1.4 s),
`certify16.py` → `certify16.json` (113 s); logs `real10/w*.log`, `gram_all.log`, `ext10.log`, `certify16.log`.
Wall-clock for the whole computation ≈ 1 h 50 min on the two cores.
