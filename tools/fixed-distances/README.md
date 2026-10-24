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
  `data/cert54_rigid.json` (481 rigid), `data/isolated54.json` (empty: no real-isolated points on the curves).
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
