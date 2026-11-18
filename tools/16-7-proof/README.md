# (16, 7) proof

A runnable proof that, up to similarity, exactly two sets of 16 points in the plane determine at
most seven distinct distances. The web version is at <https://affinemess.com/#applets/16-7-proof>.

## Run it

    python -m venv .venv
    source .venv/bin/activate            # Windows: .venv\Scripts\activate
    python -m pip install -r requirements.txt
    python -m jupyterlab

Open `two-shapes.ipynb` and choose **Run → Run All Cells** (about two minutes). The notebook already
contains its outputs, so it can also simply be read.

Python 3.10 or later. The only library the mathematics needs is SymPy (tested with 1.14).

## What is here

| file | what it is |
|---|---|
| `two-shapes.ipynb` | the notebook: explanation, code and results |
| `fewdist.py` | all the mathematics: patterns and canonical forms, generation of patterns, exact numbers, realisations, extension search, final classification |
| `figures.py` | SVG drawing helpers (no dependencies) |
| `data/L10.json` | the 228 patterns of ten points with at most six distances |
| `data/L9.json` | the 635 nine-point patterns (used for one consistency check) |
| `data/levels_summary.json` | counts for every level of the pattern generation |
| `data/solutions.json` | all solutions of the length equations of the 228 patterns (output of `fewdist.solve_many`) |

## Recomputing

The only result that is loaded instead of computed is `data/solutions.json`. Setting
`RECOMPUTE = True` in the first code cell solves the 228 systems again (about four minutes on two
cores) and overwrites the file. Everything else is computed on every run, and every saved solution
is re-verified on every run by rebuilding its ten points.

`data/L10.json`, `data/L9.json` and `data/levels_summary.json` come from the pattern generator in the
`affinemess-applets` repository (`tools/fixed-distances/g7fin/gen`). Section 3 of the notebook
explains what they are and why the proof needs only the number 228 from the literature.
