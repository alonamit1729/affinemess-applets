# fixed-distances: build pipeline

Source of `../../fixed-distances.html` (Planar sets with a fixed number of distances).
Everything is reproducible from this folder; the page is `fixed-distances.template.html` with
`data/applet_data.json` injected by `build_html.py`.

## Rebuild
    cd tools/fixed-distances
    python3 build_data2.py      # assembles data/applet_data.json from the results below (~10 min)
    python3 build_html.py       # writes fixed-distances.html next to this README; copy it to the repo root
Requires Python 3 with numpy, scipy, sympy, networkx. Scripts read/write files in the current directory, so run them
inside this folder after copying `data/*` here (or symlink).

## Where the numbers come from
- `strata.py n d seed starts` — exhaustive colouring types of K_n with d colours (up to S_n), numerical realizability
  (random-start least squares), local dimension, rigid solutions clustered by congruence. Results `data/strata_*.json`
  (`repair.py` recomputes dimensions at singular points).
- `certify54.py`, `isolated54.py` — exact certification of (5,4) with msolve + Singular: `data/cert54.json`,
  `data/cert54_rigid_exact.json` (468 rigid, exact similarity classes via `dedup54.py`), `data/isolated54.json` (empty: no real-isolated points on the curves).
- `certify5d.py d` — same pipeline for (5,d), d = 5,6,7: `data/cert55.json`, `data/cert56.json`, `data/cert57.json` (complete).
- `stubborn.py d idx` — fallback for the 21 types whose primary decomposition did not finish: dimension of the saturated
  ideal, certified real points on random slices, Jacobian rank (`data/stubborn_*.json`).
- `five.py`, `five_brute.py`, `five_id.py`, `gen5.py` — the 34 five-point and 9 six-point three-distance sets.
- `symfam.py` — symmetric constructions (cyclic/dihedral orbit structures + extra coincidences): `data/symfam_*.json`.
- `extend.py`, `examples.py`, `lattice.py` — one-point extensions, regular-polygon subsets, lattice pieces.
- `line.py n` — collinear case (flats of the difference arrangement), used for the table of mixed dimensions.
- `exactify.py` — natural frames and exact coordinate strings for the catalog.
- `lat16.py` — clique search in the triangular lattice for 16-point 7-distance pieces (`data/lattice_16_7.json`, the two
  Erdős–Fishburn examples); the same search with distance set {1,√3,2,√7,3} gives `data/lattice_11_5.json` (Wei's three).
- `data/eight4_shinohara.json` — the 15 eight-point four-distance sets assembled from Shinohara (2008), Theorem 1.2.
- `LITERATURE.md` — what each paper establishes and how the cells cite it.

## Hygiene checks (Oct 2026)
- `reverify.py d` — every unrealizable verdict re-verified independently (Singular `std` = 1 with the Rabinowitsch
  variable; msolve again for d ≥ 5): `data/reverify_{d}.json`, no disagreements (123 + 34 + 5 types).
- `dedup54.py` — exact similarity classes of the (5,4) rigid roots: 235 zero-dimensional systems, 2000 certified real
  roots, orbits of the colour-permuting automorphism group (with rescaling and the four frame symmetries) computed by
  interval arithmetic on the isolating boxes → 468 classes (`data/cert54_rigid_exact.json`, `data/dedup54.json`);
  Sturm-sequence real-root counts of msolve's rational parametrizations agree with the box counts for every system.
  (The first certified figure, 481, came from a dedup that did not rescale.)
- `noreal.py d`, `noreal2.py d` — the 21 components on which random slices found no real point: 12 provably empty
  (a generator that is a positive-definite quadratic), 9 real (hyperplanes through a numerically found point, msolve);
  no family count changes (`data/noreal2_{d}.json`, folded into `data/cert5{d}.json`).
- Hidden rigid configurations (real-isolated points of positive-dimensional types) for d ≥ 5:
  `hidden.py d` — stage 1: equations + all 10×10 Jacobian minors + inequation, msolve; types with an empty
  rank-deficient locus are done (`data/hidden_{d}.json`: 384/443, 216/249, 70/76).
  `hidden_final.py d` — stage 2 for the rest: decompose; per real non-degenerate component eliminate the linearly
  occurring variables (an isomorphism), minimal generating set (mstd), singular locus = component + codim-minors;
  smooth → done; finite → real singular points sphere-tested; positive-dimensional → recurse
  (`data/hidden_final_{d}.json`). `hidden_subsets.py` — smoothness via minors of random generator subsets (a
  superset of the singular locus), used for the two components of d=6 type 242 where Singular crashes and for the
  saturated ideals of the four stubborn d=5 types. `nodetest.py idx` — for the four d=5 types without a decomposition (380, 381, 426, 427): certified local test at each
  candidate point p (interval arithmetic, 300 bits, on msolve's isolating box): generators with provably independent
  gradients cut out a smooth surface S ∋ p containing the curve C; for a further generator g, det of the intrinsic
  Hessian of g|S at p is certified negative, so {g|S = 0} is a smooth arc or an ordinary node with two real branches
  (Morse lemma), and C's purely one-dimensional germ at p lies on real branches. All 20 candidates pass.
  Result: no isolated real point anywhere for d = 5, 6, 7.
  Superseded helpers kept for reference: `hidden2.py`, `hidden2b.py`, `hidden2c.py`, `hidden3.py`, `hidden_242.py`,
  `hidden_stubborn.py`.

**msolve input rule.** msolve accepts only expanded polynomials with integer (or standalone rational) coefficients. It
does not report an error on input such as `(x-1/2)^2` or `3*y/7`: it silently returns "no solutions". Every polynomial
handed to msolve in these scripts is therefore expanded by sympy and, where rational centres occur (sphere tests),
denominators are cleared with `Poly.clear_denoms()`. The first version of the sphere test in `hidden.py` violated this
and reported a non-isolated point as isolated; it was caught by an independent geometric check (the zigzag family of
type 28 at d = 5) and fixed.

**Concurrency rule.** All msolve/Singular scratch files are now process-unique (`tmp/ms_in_<pid>.txt` etc., see
`certify5d.py`). Earlier runs of the hygiene scripts shared `ms_in5.txt`/`ms_out5.txt`; results produced while two such
jobs overlapped were discarded and recomputed by `hygiene_run.sh`, which runs everything strictly sequentially.

## g(7) = 16 review (`g7/`)

`g7/ext11.py` — exact replication of Lemma 7 of Amit–Astra 2026 for the lattice cores (`g7/cores.json` holds the eleven twelve-point lattice cores L0–L10 transcribed from the paper's appendix): every 11-point one-point deletion with ≤ 6 distances, extended by points whose squared distances to the non-extreme anchors lie in the palette plus at most one (two, for five-distance cores) new value; prints the maximum extension size per core (16).
`g7/latt11b.py tri|sq MMAX RATIO2` — independent enumeration of all 11- and 12-point lattice subsets with ≤ 6 distances up to similarity (origin + a minimum-distance partner fixed, DFS with the six-norm constraint; scope: minimum norm ≤ MMAX, squared diameter ≤ RATIO2 × minimum norm), compared with the deletions of the cores. Logs: `g7/latt_tri.log` (11 and 52, matching the paper), `g7/latt_sq.log` (none).

## Marks (provenance of completeness)

Unmarked solid cell: published classification (configurations reconstructed here). `†`: published and independently re-verified by our exact computation. `‡`: our exact computation only (novel). `certify_nd.py n d` is the n-general version of `certify5d.py`; it was run for (4,2..6), (5,2), (5,3) (`data/cert{n}{d}.json`): counts 6 / 9 / 7 / 2 / 1 families-or-rigid, 1, 34 — all agreeing with the literature and the earlier numerics. Note that with the scale fixed by s0 = 1 the variety has the mirror pair x1 = ±1, so component counts per type come in pairs, and rigid roots must be deduplicated up to similarity including rescaling (colour permutations).

## Stratification of the five-point solution sets (`strat*.py`, `degen*.py`)

`sat5b.py d` — saturated ideal (closure of the non-degenerate locus, x1 = 1) of every positive-dimensional type of (5,d), factor-by-factor saturation in Singular (`sat5{d}.json`; `sat5.py` is the slow single-product version).
`strat.py d` — for each one-dimensional type: exact boundary roots (saturated ideal + one coincidence equation, msolve 200 bits), identified as rigid sets of fewer distances / degenerate (two points coincide); the real curve is traced from every root (tangent step + Newton projection with coefficient-normalised generators and an rcond cutoff on the tangent direction; arrival by proximity to the exact roots; escape to infinity when coordinates exceed 300) and from random slice points (closed loops); traversals are grouped into geometric arcs by a permutation- and similarity-invariant fingerprint (sorted distances / diameter) compared over the non-degenerate shape range. Checkpointed, atomic writes (`strat_{d}.json`). `strat_mixed.py` handles the three (5,5) types with both one- and two-dimensional components.
`strat_analyze.py d` — shape-space incidence graph: vertices = rigid sets (identified in the catalogs) and junction points on lower families; edges = arcs; prints valences, connectivity, first Betti number; writes `strat_graph_{d}.json` which `build_data2.py` folds into the cards and cell notes.
Results: (5,4): 35 vertices (34 three-distance sets + pentagon), 108 arcs, 77 degenerate ends, 12 ends at infinity, connected, b1 = 74; the 468 rigid sets are isolated. (5,5): 527 vertices (468 + 34 + 1 rigid, 24 junction points on (5,4) families), 5269 arcs, ~1900 degenerate ends, 105 at infinity, connected, b1 = 4743 (the 3.3 % of arcs with odd multiplicity are almost all half-open degenerate ends; 17 vertex–vertex arcs are uncertain by one). Type 327 has 5 lost traces.
`degen.py d` / `degen_sample.py d` / `degen_analyze2.py d` — degeneration poset: for each positive-dimensional type and each coincidence equation, the components of (saturated ideal + equation) over Q (minAssGTZ, 150 s timeout) with the colouring type of a certified real point, or, as a fallback, the types met by random slices (msolve). Exact for d = 4, 5 (one timeout) and 81 types of d = 6; sampled for the rest (`degen_poset_{d}.json`). Caveat: random small-integer slices hit symmetric special points with positive probability, so sampled lists may include codimension-≥2 strata.

## Three and four points; the summary page

`strat_n.py d n`, `sat_n.py n d`, `degen_n.py n d` are the n-general versions (n = 3, 4 run in seconds). Results: (3,2): one curve, an arc from the equilateral triangle to a collision and an arc to the fold at the collinear triple with the apex at the midpoint. (4,3): 7 vertices (6 two-distance quadruples + the fold at four equally spaced collinear points), 18 arcs, 15 collision ends, 2 at infinity, connected, b1 = 12; `degen_4_4.json` lists the boundary curves of the seven surfaces. Folds (a family meeting its own mirror image, i.e. crossing the collinear/mirror locus) are detected by a palindromic shape sequence; the arc then ends at the mirror-symmetric configuration. Collision ends seen as 'two points merge' and as 'escape to infinity' are now identified by their shape fingerprints over the common range.
`os_inv.py n` — rational Betti numbers of the shape spaces (S_n-invariants of the Orlik–Solomon algebra of the braid arrangement; conjugation acts by (−1)^k).
`summary/` — generator of the summary page "Five-point shape strata" (figures as inline SVG).

## Finiteness of the maximum sets (`g7fin/`)

`g7fin/FINITENESS.md` — the argument that there are only finitely many 16-point seven-distance sets (Wei's core lemma, a pinning lemma, certified rigidity of every 10- and 11-point pattern with ≤ 6 distances), with its state. `g7fin/gen/` — a re-implementation of the Szöllősi–Östergård orderly generation of candidate patterns (nauty canonical forms, Gröbner-basis realizability test mod p); it reproduces Table 9 of their paper at every level (25/1193/537230/14732/1565/635/228 raw; 24/922/20229/4667/1565/635/228 pruned), and the three order-10 patterns outside the known families are real (`g7fin/gen/new3/`) and rigid (`g7fin/pat_rigid_10new.json`). `g7fin/cases/` — Theorem A (no 11-point set has diameter graph C₁₁ and exactly seven distances) and the leaf/cycle lemmas for diameter graphs. `symtest/` — the test of "flexible ⇒ symmetric" (false; 15 counterexamples, all point-deletions of symmetric families).

## The twelvefold set R₁₂ + R₁₂ (`r12/`)

The three new ten-point patterns are subsets of R₁₂ + R₁₂ = {ζ^a + ζ^b : ζ¹² = 1} (73 points). `r12/sweep.py D` enumerates every subset with ≤ D distances up to similarity (DFS with the distance-count bound, 52 s for D = 7); `r12/export_cells.py` writes `data/r12sum.json`, the subsets for the cells (8–12, 5), (9–13, 6), (13–16, 7), which `build_data2.py` merges into the catalogs (label "piece of the twelvefold set", exact coordinates in Q(√3)). Counts of similarity classes in the sweep: (7,4) 34 ⊂ 42, (8,4) 12 ⊂ 15, (9,4) 3 ⊂ 4, (11,5) 3 ⊂ 4, (12,5) 1, (13,6) 2 ⊂ 3, (16,7) 2 = 2 (all consistent with the classified cells); (8,5) 131, (9,5) 48, (10,5) 15, (9,6) 526, (10,6) 195, (11,6) 53, (12,6) 11, (13,7) 299, (14,7) 77, (15,7) 14; lower rows (thousands) are not merged. The sweep exposed a bug in the catalog's similarity test (`build_data.same` compared only the multiset of distances, which identifies homometric sets); it is now a true congruence test (vertex bijection by backtracking). No previously built cell changed.

## Exact classification of the top cells; the two sixteen-point sets (`g7count/`)

With the complete candidate lists (`g7fin/gen/levels/L{7..13}.json`, counts 4635/1565/635/228/62/15/3 = Szöllősi–Östergård's) every pattern of 10–13 points was solved by msolve real root isolation — `g7count/gram_n.py N` uses their Gram-minor system (3×3 minors of the Gram matrix in the ≤ 5 free values, ~1 s per pattern, PSD rank-2 check and coordinates at 100 digits; `realize_all.py` is the coordinate-system route, ~1000× slower, used as a cross-check) — and the realizations were counted up to similarity by congruence (`sets10_by_d.json`, `sets1113_by_d.json`): (10,5) 19, (10,6) 238, (11,5) 4, (11,6) 64, (12,5) 1, (12,6) 14, (13,6) 3. `exact_cells.py` labels them (polygon subsets, lattice pieces) into `data/exact_cells.json`, which `build_data2.py` uses for the solid cells (10,5), (10,6), (11,6), (12,6). `g7count/REPORT.md`: every 16-point seven-distance set contains one of the 257 ten-point ≤6-distance sets or is R₁₄/R₁₅/R₁₅−v plus 1–2 points; `complete10.py`/`extend.py` extend each core by six points (circle intersections and circumcentres as candidates, DFS with a budget of ≤ 2 new values), `polygons.py` settles the polygon cases, `certify16.py` certifies the two resulting patterns: **exactly two sixteen-point seven-distance sets** (Erdős–Fishburn's).
