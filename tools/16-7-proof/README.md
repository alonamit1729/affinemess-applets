# (16, 7) proof

A runnable proof that, up to similarity, exactly two sets of 16 points in the plane determine at
most seven distinct distances. The web version is at <https://affinemess.com/#applets/16-7-proof>.

## Run it

    python -m venv .venv
    source .venv/bin/activate            # Windows: .venv\Scripts\activate
    python -m pip install -r requirements.txt
    python -m jupyterlab

Open `two-shapes.ipynb` and choose **Run → Run All Cells** (two to three minutes). The notebook already
contains its outputs, so it can also simply be read.

Python 3.10 or later. The only library the mathematics needs is SymPy (tested with 1.14).

## What is here

| file | what it is |
|---|---|
| `two-shapes.ipynb` | the notebook: explanation, code and results |
| `fewdist.py` | all the mathematics: patterns and canonical forms, generation of patterns, exact numbers, realisations, extension search, final classification |
| `figures.py` | SVG drawing helpers (no dependencies) |
| `data/L10.json` | the 228 patterns of ten points with at most six distances, in the numbering the notebook uses |
| `data/solutions.json` | all solutions of the length equations of the 228 patterns (output of `fewdist.solve_many`) |
| `certificates/` | how the 228 are derived: the lists of patterns kept at five to ten points, the program that checks certificates (`check.py`), the program that lists all candidates (`generate.py`), the five-point certificates and samples of the six- and seven-point ones; see `certificates/README.md` |

## Recomputing

Two results are loaded instead of computed when the notebook runs.

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
