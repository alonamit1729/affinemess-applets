# FLINT cross-check of the solutions of the 228 patterns

The paper (section 4) and the notebook (section 4) find all solutions of the 228 ten-point patterns with
SymPy's Gröbner bases. This folder repeats the completeness half with a different implementation, the
Gröbner engine of FLINT (python-flint 0.9), and with more equations: the planarity condition for every
four of the ten points (3x3 Gram minors at every base point), plus `w*q - 1` for "all lengths non-zero and
pairwise different". For each pattern it checks that this system and the saved solution set have the
same reduced Gröbner basis, so they have the same solutions.

| File | What it is |
|---|---|
| `flint_all_origins.py` | the check for all 228 patterns (`python3 flint_all_origins.py .`, about ten minutes on two cores) |
| `flint_all_origins.log` | its output: equality for 224 patterns; 14, 25, 66 and 219 time out |
| `flint_hard.py` | the four remaining patterns with a lighter saturation, still a valid completeness proof (`python3 flint_hard.py . 14 x4 600`, and `none` for 25, 66, 219) |
| `flint_hard.log` | its output: equality for all four; with them, 282 solutions in total |
| `L10.json`, `palettes.log` | the inputs: the 228 patterns (the same list as `../data/L10.json`) and the saved solution sets, one reduced basis per pattern |

These files were written on 15 Nov 2026 to review an earlier version of the proof, made with another
model, whose saved solutions are `palettes.log`; the solution sets are the same as in `../data/solutions.json`
(282 solutions, distributed over the patterns as in Table 2 of the paper). Needs `sympy` and `python-flint`.
