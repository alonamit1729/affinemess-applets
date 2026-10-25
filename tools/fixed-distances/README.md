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
  saturated ideals of the four stubborn d=5 types. Result: no isolated real point anywhere, except that for d=5
  types 380, 381, 426, 427 (no decomposition) the 20 candidate points are non-isolated only numerically.
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
