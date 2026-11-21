# construction-tool: build pipeline

The Construction Tool (formerly Distance Hunt; this folder was `tools/distance-hunt/` until 2026-11-20).
`../../construction-tool.html` is `construction-tool.template.html` with `catalog.json` injected:

    python3 build.py                       # writes ../../construction-tool.html

`../../distance-hunt.html` and `../../card-viewer.html` only forward to it (keeping the hash), so old links,
collection links and the homepage tile keep working. The Card Viewer became the tool's *view mode*:
`construction-tool.html#n<n>d<d>c<card>` opens a catalog card as a construction (Show construction, Play,
step, Edit). `SITE_SLUG` in the template names the homepage tile used for links on the site (`distance-hunt`
until the Website Project adds a `construction-tool` tile; the fixed-distances page has the same switch,
`TOOL_SLUG`).

## Catalog constructions (`recipes/`, `pipeline/`)

`recipes/<n>_<d>/<k>.json` holds a construction for every card of the fixed-distances catalog (cards k·1000 to
k·1000+999 of the cell; the page fetches the file it needs): `{cell, labels, items: {num: [kind, dim, level,
equations, unknowns, label index, recipe, verdict?]}}`. Levels: 1 a piece of the triangular or square lattice,
2 straightedge and compass (regular polygons that straightedge and compass can draw count), 3 straightedge,
compass and another regular polygon, 4 sliders / free points pinned by equalities (the equalizer), 0 none found
(the recipe then just places the points). A recipe is the tool's own node list, one short array per node
(`decodeRecipe()` in the template), with the set's points (`m`), the equalities (`q`) and the lattice (`L`).

    cd pipeline
    python3 export_jobs.py /tmp/jobs                    # every card, coordinates refined onto its pattern (35 s)
    node run.js /tmp/jobs /tmp/out 5_4 6_6 ...          # the search (all 46 cells: about 15 min on 2 cores)
    node write_recipes.js /tmp/jobs /tmp/out            # -> ../recipes/ (a file is rewritten only if it changed)

`run.js` merges into what `/tmp/out` already holds, so a cell can be done in pieces (`--from= --to=`,
`--seconds=` for a time budget, `--redo=l4` to search again where only the equalizer was found; `DEEP=1` adds
regular polygons through just two points of the set as starts, slow: about 0.5 s per set).

`construct.js` is the search; `core.js` cuts the geometry functions (`evaluate`, `meet`, `canon`,
`decodeRecipe`, …) out of the page template verbatim, so every construction is built and checked with exactly
the code the browser replays it with. A construction is kept only if replaying it gives the catalog set (same
canonical form) with exactly d distances. The search: starts (a pair of the set's points, a partial lattice,
the catalog's own polygon frame, regular polygons through three or more of the set's points), then points of
the set where two available lines or circles meet (lines through two points, perpendicular bisectors, circles
about a point through another, compass circles of an existing radius); when stuck, a slider on a line or circle
through a missing point (a free point if there is none), and at the end the fewest equalities that pin the
sliders, or leave exactly the family's freedom. Only points of the set and of the start are stepping stones,
so it misses constructions that need helper points.

Results of 2026-11-20 (61,172 cards): rigid sets 4,306 lattice, 18,874 straightedge and compass, 3,716 regular
polygon, 32,483 equalizer, 3 none; families 1,435 straightedge and compass (with free points), 3 polygon, 351
equalizer, 1 none. Against the Galois verdicts (2,298 sets): no set is claimed constructible with tools its
verdict rules out; 66 sets with verdict "straightedge and compass" fall back to the equalizer.

## Fingerprints and verdicts

`catalog.json` holds a similarity fingerprint of every rigid set in the fixed-distances
catalog (59,377 sets on 2026-11-06, reflections identified). Cells with more than 1,000
rigid sets ((6,5), (6,6), (7,6), (8,6)) go to `catalogs/<n>_<d>.json`, which the applet fetches when
first needed. Rebuild when that catalog changes, passing the previous `catalog.json` (with its
`catalogs/` folder beside it) so the constructibility verdicts carry over by shape (card numbers
may shift between versions; compare before publishing):

    cp catalog.json /tmp/catalog.old.json
    python3 canon.py ../fixed-distances/data/applet_data.json catalog.json /tmp/catalog.old.json
(`applet_data.json` is slim since 2026-11-07: `canon.py` restores the big cells' items from `../fixed-distances/data/catalogs/` through `datautil.py`; the result is identical.)

`canon.py` polishes each catalog set by Newton's method on its distance pattern (`refine.py`;
nine (6,3) entries of the source data are off by ~1e-8) and computes the canonical form:
over the ordered pairs at minimum distance, with and without reflection, send the pair to
0 and 1 and sort the points (tolerance 1e-7); keep the lexicographically smallest list.
The game computes the same form in JavaScript (`canon()` in the template) and matches with
tolerance 1e-5.

`verdicts.json` lists, for every catalog card (cell, card number), what it takes to construct
it (null: not analysed; the analysis covers the 2,298 sets of 2026-11-03): `rc` straightedge and compass, `p<m>` a regular m-gon, `tri` angle trisection, `s4` a
quartic with Galois group S4 and two real roots, `S5`/`S6`/`S7` not constructible (non-solvable
Galois group), `cbrt` not constructible (real root of a cubic with one real root). It was made
by the scripts in `analysis/` (run from a folder holding `rigid.json`, the rigid catalog sets
extracted from `../fixed-distances/data/applet_data.json`): `fields.py` refines each set to
1,300 digits and finds the minimal polynomial of every squared-distance ratio by LLL, probing
its Galois group by factoring mod 1,500 primes; `groups.py` identifies degree 4–6 groups with
sympy; `chain.py` tests practical constructibility (lattice or polygon start, then circles of
existing radii, lines, bisectors); `final.py` and `cross.py` combine the two. `canon.py` folds
the verdicts and the explorer's card numbers into `catalog.json`.

`analysis/slide.py` checks which catalog sets one slide reaches (a point of the set on an
existing line or circle, the rest built from it, one equality between two of the set's distances).
Run on the 200 sets ordinary construction cannot reach: 197 need one slide; (5,4) nos. 461, 463
and 466 need several sliders solved together (`analysis/slide_out.json`).
