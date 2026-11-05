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
