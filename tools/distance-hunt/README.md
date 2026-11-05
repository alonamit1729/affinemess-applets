# distance-hunt: build pipeline

`../../distance-hunt.html` is `distance-hunt.template.html` with `catalog.json` injected:

    python3 build.py                       # writes ../../distance-hunt.html

`catalog.json` holds a similarity fingerprint of every rigid set in the fixed-distances
catalog (2,298 sets, reflections identified). Rebuild it when that catalog changes:

    python3 canon.py ../fixed-distances/data/applet_data.json catalog.json

`canon.py` polishes each catalog set by Newton's method on its distance pattern (`refine.py`;
nine (6,3) entries of the source data are off by ~1e-8) and computes the canonical form:
over the ordered pairs at minimum distance, with and without reflection, send the pair to
0 and 1 and sort the points (tolerance 1e-7); keep the lexicographically smallest list.
The game computes the same form in JavaScript (`canon()` in the template) and matches with
tolerance 1e-5.

`verdicts.json` lists, for every catalog card (cell, card number), what it takes to construct
it: `rc` straightedge and compass, `p<m>` a regular m-gon, `tri` angle trisection, `s4` a
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
