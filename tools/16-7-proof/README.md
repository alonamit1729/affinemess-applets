# (16, 7) proof

A runnable proof that, up to similarity, exactly two sets of 16 points in the plane determine at
most seven distinct distances. The web version is at <https://affinemess.com/#applets/16-7-proof>.

## Run it

    python -m venv .venv
    source .venv/bin/activate            # Windows: .venv\Scripts\activate
    python -m pip install -r requirements.txt
    python -m jupyterlab

Open `two-shapes.ipynb` and choose **Run → Run All Cells** (about three and a half minutes). The notebook already
contains its outputs, so it can also simply be read.

Python 3.10 or later. The only library the mathematics needs is SymPy (tested with 1.14).

## The paper

`paper/seven-distances.tex` (and the compiled `paper/seven-distances.pdf`) writes the proof up as a
self-contained paper: g(7) = 16, exactly two 16-point sets with seven distances, and exactly sixteen
15-point sets. It cites this folder for the computations; the notebook follows it, including the
elementary proof of the thirteen-cycle lemma.

## What is here

| file | what it is |
|---|---|
| `two-shapes.ipynb` | the notebook: explanation, code and results |
| `fewdist.py` | all the mathematics: patterns and canonical forms, generation of patterns, exact numbers, realisations, extension search, final classification |
| `figures.py` | SVG drawing helpers (no dependencies) |
| `data/L10.json` | the 228 patterns of ten points with at most six distances, in the numbering the notebook uses |
| `data/solutions.json` | all solutions of the length equations of the 228 patterns (output of `fewdist.solve_many`) |
| `fifteen/` | the classification of fifteen-point sets with seven distances (notebook §8, paper §8): `fifteen.py`, the solutions of the 635 nine-point patterns (`solutions9.json`) and the 15 lattice sets found (`result.json`) |
| `paper/` | the paper, `seven-distances.tex` and `.pdf` |
| `certificates/` | how the 228 are derived: the lists of patterns kept at five to ten points, the program that checks certificates (`check.py`), the program that lists all candidates (`generate.py`), the five-point certificates and samples of the six- and seven-point ones; see `certificates/README.md` |

## Recomputing

Three results are loaded instead of computed when the notebook runs.

* `data/solutions.json`. Setting `RECOMPUTE = True` in the first code cell solves the 228 systems again
  (about four minutes on two cores) and overwrites the file. Every saved solution is re-verified on every
  run by rebuilding its ten points.
* The lists of patterns on six to ten points, and the certificates for the patterns that are discarded
  on the way (section 3 of the notebook). The two large certificate files, `level6.jsonl.gz` and
  `level7.jsonl.gz`, are attached to the release `certificates-v1` of the `affinemess-applets`
  repository; `certificates/README.md` has the links. Download them into `certificates/`, then

      cd certificates
      python generate.py 6 --jobs 2                                 # and 7, 8, 9, 10: 13 minutes in all
      python check.py --jobs 4 level6.jsonl.gz level7.jsonl.gz      # 45 core-minutes

  The first command lists all candidates and checks that each one is kept or has a certificate; the
  second checks every certificate.
* The fifteen-point classification (section 8 of the notebook). The notebook repeats it for a sample of
  the nine-point sets; the full run is

      python fifteen/fifteen.py --jobs 2      # 16 minutes; --resolve also re-solves the 635 patterns
