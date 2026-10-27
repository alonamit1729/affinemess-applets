# Reproducing the classification of sixteen-point seven-distance sets

**Claim.** Up to similarity there are exactly two sets of sixteen points in the plane that determine exactly seven distinct distances: the two triangular-lattice pieces of Erdős and Fishburn (the 12-point block of three rows of four, with a row of two above and a row of two below, the short rows on the same side or on opposite sides).

This note explains, for someone who wants to check the claim rather than take it on trust, what has to be believed, what has to be run, and what the runs must print. All code and data are in the folder `tools/fixed-distances/` of the repository `github.com/alonamit1729/affinemess-applets`; paths below are relative to it. The computation has six steps, of which one (the candidate patterns) takes about half a day on one core and the others minutes.

## 0. Shape of the argument

Let X be a sixteen-point set with exactly seven distances.

1. **Reduction (§1, a short proof).** X contains a ten-point subset with at most six distances, or X consists of a regular 14- or 15-gon (or a regular 15-gon minus a vertex) and one or two more points.
2. **Candidate patterns (§3, computation).** The distance pattern of any ten-point set with at most six distances — the edge colouring of K₁₀ by distance, up to permuting points and colours — is one of 228 patterns. This is the orderly generation of Szöllősi and Östergård, re-implemented and re-run; the count 228 is theirs.
3. **Rigidity (§4, computation).** Each of the 228 patterns has finitely many realizations: the polynomial system of each is zero-dimensional over ℚ (Singular).
4. **Realizations (§5, computation).** Solving the 228 systems by real root isolation (msolve) gives all their real solutions: 257 ten-point sets up to similarity.
5. **Extension (§6, computation).** Each of the 257 sets is extended by six points in every way compatible with "seven distances in all"; the extra points are among finitely many candidates (intersections of circles, circumcentres), and a depth-first search finds every valid 6-set. Result: 826 labelled sixteen-point sets, two up to similarity. The polygon cases of step 1 are settled by a direct check (§7).
6. **Certification (§8, computation).** The two patterns found are solved exactly (msolve): each is rigid with a unique realization, the Erdős–Fishburn sets.

Steps 2–4 and 6 are exact (Gröbner bases modulo a prime with confirmation, Singular over ℚ, msolve's certified isolating intervals). Step 5 is carried out in 100-digit floating point; §9 explains why this cannot lose a set and why anything it could add spuriously is caught by step 6.

## 1. The reduction

Notation: D is the diameter of X, X_D the set of diameter endpoints, m = |X_D|, G_D the *diameter graph* (vertices X_D, edges the pairs at distance D), α = α(G_D) its independence number.

**Facts about G_D** (all elementary; proofs in `g7count/REPORT.md` §1). Every diameter endpoint is an extreme point of X; any two diameter segments meet; G_D contains no even cycle and at most one cycle, and every edge meets that cycle. Consequently G_D is a forest or an odd cycle with pendant trees, and α ≥ ⌈m/2⌉ unless G_D is the bare cycle C_m, where α = (m − 1)/2. (This is the diameter-graph lemma of Wei, Electron. J. Combin. 19(4) (2012) P38.)

**A core.** If I ⊂ X_D is independent in G_D, then A = I ∪ (X ∖ X_D) contains no pair at distance D, so A determines at most six distances, and |A| = 16 − m + α.

| m | α | |A| ≥ | conclusion |
|---|---|---|---|
| ≤ 10 | ≥ (m−1)/2 | 11 | A contains a ten-point set with ≤ 6 distances |
| 11 | 5 if G_D = C₁₁, else ≥ 6 | 10 | same |
| 12 | ≥ 6 | 10 | same |
| 13, G_D ≠ C₁₃ | ≥ 7 | 10 | same |
| 13, G_D = C₁₃ | 6 | 9 | X_D is the regular 13-gon (Amit–Astra 2026, Lemma 4: thirteen points with diameter graph C₁₃ and at most seven distances form R₁₃); R₁₃ has six distances, so any ten of its vertices form a ten-point set with ≤ 6 distances |
| 14 | — | — | X_D is a convex 14-gon (diameter endpoints are extreme points) with ≤ 7 distances; a convex n-gon has ≥ ⌊n/2⌋ distances (Altman 1963), so exactly 7, and a convex 14-gon with exactly 7 distances is R₁₄ or R₁₅ minus a vertex (Erdős–Fishburn 1996 and Fishburn, Comput. Geom. 5 (1995); the statement as used is Lemma 4 of Wei 2012: for even n > 8 and t = n/2, R_n or R_{n+1} − 1; for odd n and t = (n − 1)/2, R_n). X = X_D ∪ {q₁, q₂}, and all seven distances already occur in X_D, so every distance from qᵢ to X_D is a chord length of X_D |
| 15 | — | — | X_D = R₁₅ (same theorem, n odd); X = R₁₅ ∪ {q}, same condition on q |
| 16 | — | — | a convex 16-gon has ≥ 8 distances: impossible |

So X contains a ten-point subset with at most six distances, or (m ≥ 14) X is R₁₄, R₁₅ − v or R₁₅ plus one or two points all of whose distances to the polygon are chord lengths of it. The second alternative is checked in §7.

Three cited results enter: Altman's bound, the Erdős–Fishburn/Fishburn classification of convex n-gons with ⌊n/2⌋ distances, and Lemma 4 of Amit–Astra (whose proof is itself computer-assisted, with its own certificate). Nothing else from the literature is used: in particular neither the bound g(7) = 16 nor the finiteness argument of `g7fin/FINITENESS.md` is needed for the count.

## 2. Environment

Python 3.13 with numpy, sympy, networkx, mpmath 1.3 (with gmpy2 for speed), and pynauty 2.8.8 (the generator's canonical forms; `g7fin/gen/ext/` has the sources if the wheel does not build). msolve 0.10.1 (Berthomieu, Eder, Safey El Din) and Singular 4.4 on the path; the scripts look for them at the paths in `g7fin/gen/realize.py`, `g7fin/pat_rigid_d.py` and `g7count/msolve_util.py` (edit the `MS`/`SG` constants if yours differ). Everything was run on one Linux core; timings below are from that machine.

Scripts read and write relative to the folder they are described in, so `cd tools/fixed-distances` first; `g7count/` scripts are run from `g7count/` and find `g7fin/` through a relative path.

## 3. Step 2 — the 228 candidate patterns (`g7fin/gen/`)

**What is computed.** A *pattern* on n points is an edge colouring of K_n with at most six colours, up to permutations of points and of colours. The list L_n is built level by level from n = 4: raw_n is the set of colourings of K_n all of whose (n−1)-point induced sub-colourings lie in L_{n−1} (at n = 4, all 25 colourings), and L_n ⊂ raw_n keeps those that are *complex-realizable with distinct nonzero values*: the ideal generated by the rank-≤2 conditions on the Gram matrix of the squared distances (all 3×3 minors) together with the Rabinowitsch equation t·∏_{a<b}(s_a − s_b)·∏_c s_c − 1 (values pairwise distinct and nonzero) is proper, i.e. its Gröbner basis is not {1}. This is Proposition 4.2 of Szöllősi–Östergård (Electron. J. Combin. 27 (2020) P1.23, with d = 2) applied hereditarily, and it reproduces the two rows of their Table 9 exactly. Since the pattern of a real ten-point set with ≤ 6 distances satisfies both conditions at every level, it lies in L₁₀.

The Gröbner test is done modulo the prime 1073741827 (msolve `-g 2`); a pattern is discarded only when the basis is {1}, and by default each discard is confirmed modulo a second prime (`--no-confirm` skips that; the run reported below used it, for speed, and the agreement with the published counts at every level is the check). Keeping a false candidate is harmless — it has no real realization and drops out in step 4.

**Run.**

    cd tools/fixed-distances/g7fin/gen
    python3 gen.py --to 13 --out levels --workers 1        # resumable: rerun the same command after an interruption

Runtime about 12 h on one core (level 6: 537 230 candidates to test, 2.7 h; level-7 generation 2 h; the rest minutes). The run writes `levels/L{n}.json` (`{"n","count","patterns":[...]}`, each pattern a colour tuple indexed by `itertools.combinations(range(n),2)`), `raw{n}.json`, checkpoints, `summary.json` and `gen.log`.

**Expected output** (`summary.json`; the two rows of Szöllősi–Östergård's Table 9):

| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|
| raw | 25 | 1193 | 537 230 | 14 732 | 1565 | 635 | 228 | 62 | 15 | 3 |
| L_n | 24 | 922 | 20 229 | 4667 | 1565 | 635 | 228 | 62 | 15 | 3 |

`python3 validate.py --levels levels` prints this comparison and also checks L₄, L₅ against the applet's independently certified four- and five-point classifications. The committed `levels/L{7..13}.json` are the lists used below; a verifier who does not want to spend the half day can at least confirm the committed counts and sub-pattern checks with `validate.py`, but regenerating is the real check.

## 4. Step 3 — every pattern is rigid (`g7fin/`)

**What is computed.** For each of the 228 patterns, the coordinate system: p₀ = (0,0), p₁ = (x₁,0), the colour of edge (0,1) has squared value 1, one free value s_c per other colour, one equation |pᵢ − pⱼ|² = s_{colour(ij)} per edge, and the Rabinowitsch equation. The ideal is screened for zero-dimensionality with msolve's Gröbner basis modulo a prime (zero-dimensional iff a pure power of every variable occurs among the leading monomials) and then certified with Singular's `dim(std(I))` over ℚ.

**Run.**

    cd tools/fixed-distances
    python3 g7fin/extract_d.py 10 6          # the 228 patterns of L10 in pat_rigid format -> g7fin/patterns_10d6.json
    python3 g7fin/pat_rigid_d.py 10 6        # -> g7fin/pat_rigid_10d6.json  (about 15 min)

**Expected output.** All 228 records have `status: zerodim` and `dimQ: 0`. (The same scripts with `11 6`, `12 6`, `13 6` certify the 62, 15 and 3 patterns of the next levels; not needed for the claim.)

## 5. Step 4 — all real realizations: 257 ten-point sets (`g7count/`)

**What is computed.** For each pattern, msolve isolates the real solutions of Szöllősi–Östergård's own system — the 3×3 minors of the (n−1)×(n−1) Gram matrix C_ij = (s_0i + s_0j − s_ij)/2 in the ≤ 5 free values, plus the Rabinowitsch equation — to 400 bits (`gram_util.py`, `gram_n.py`). For each real solution the Gram matrix is formed at 100 digits and must be positive semidefinite of rank 2 (eigenvalues: two clearly positive, the rest ≈ 10⁻¹⁰⁰; this is the exact criterion for a Euclidean realization in the plane); coordinates come from its eigendecomposition, and the script asserts that the reconstructed set has exactly the pattern it came from. Because step 3 says the system is zero-dimensional, msolve's output is the complete list of real solutions.

**Run.**

    cd tools/fixed-distances/g7count
    python3 gram_n.py 10                     # -> real10g/pat_<i>.json, 228 files, about 3 min

**Expected output.** Every pattern has `status: real`; the number of solutions per pattern is 1 for 199 patterns, 2 for 20, 3 for 1, 4 for 3, 5 for 2 and 6 for 3 (282 labelled sets); up to similarity these are **257** ten-point sets (214 patterns with one realization, 9 with two, 2 with four, 1 with five, 2 with six). `python3 gram_to_real10.py real10_gram` converts them to the extension's input format and prints `converted 228 patterns`.

**Cross-checks that were run** (optional). `realize_all.py` solves the much larger coordinate system of §4 instead (218 of the 228 patterns in 5554 s; ten are too heavy for a 4 GB machine) and gives congruent realizations class by class; `crosscheck.py` verifies that all 210 ten-point sets with ≤ 6 distances inside the twelvefold set R₁₂ + R₁₂ (`r12/r12_subsets_7.json`) and all 236 ten-point subsets of the 68 eleven-point sets of Amit–Astra's Lemma 6 are among the 257.

## 6. Step 5 — extending the 257 sets to sixteen points (`g7count/`)

**What is computed** (`extend.py`). Let A be one of the 257 sets, Δ(A) its set of squared distances (5 or 6 values; scale so that the smallest is 1), and b = 7 − |Δ(A)| the number of *new* values allowed. A point q of X ∖ A has ten distances to A, all in Δ(A) ∪ N with |N| ≤ b ≤ 2. Either at least two of them are in Δ(A), so q is an intersection point of two circles centred at points of A with radii in Δ(A); or at most one is, so at least nine lie in N, some value occurs at least five times, and q is the circumcentre of three non-collinear points of A. `candidates(A, b)` enumerates both finite families (about 3100 circle points and 120 circumcentres per core), discards points of A, deduplicates, and keeps the points q with |{d(q,a) : a ∈ A} ∖ Δ(A)| ≤ b — at most 16 per core. `extensions(A, 6, b)` is a depth-first search over these survivors in a fixed order that adds a candidate only if its distances to the ones already chosen stay in Δ(A) ∪ N and the running set N of new values never exceeds b; every complete 6-set is recorded. This enumerates every sixteen-point superset of A with ≤ 7 distances.

**Run.**

    cd tools/fixed-distances/g7count
    REAL10DIR=real10_gram EXT10DIR=ext10_gram python3 complete10.py 0 1     # -> ext10_gram/pat_<i>.json, about 1 min
    EXT10DIR=ext10_gram python3 collect.py                                   # -> complete_result.json

**Expected output** (`collect.py`): 257 cores extended; at most 16 surviving candidates and at most 1084 search nodes per core; **826 labelled 6-extensions from 185 cores**, every one with exactly seven distances; **2 similarity classes**: set 0 with 694 hits from 184 cores and set 1 with 132 hits from 89 cores, both reported congruent to the known Erdős–Fishburn sets. The hit counts are themselves a check: 184 and 89 are exactly the numbers of similarity classes of ten-point subsets with ≤ 6 distances of the two sets — every admissible core of each set produced it, and nothing produced anything else.

## 7. The polygon cases (`g7count/polygons.py`, 1 s)

    python3 polygons.py

For X_D ∈ {R₁₄, R₁₅ − v, R₁₅} every extra point has all its distances to X_D among the chord lengths of X_D, hence at least two of them, hence is an intersection of two circles with chord-length radii centred at vertices: `candidates(X_D, 0)` lists all of them. **Expected output:** R₁₄: no point; R₁₅: no point; R₁₅ − v: exactly one point, the missing vertex v (which gives R₁₅, fifteen points). So no sixteen-point seven-distance set has fourteen or more diameter endpoints.

## 8. Step 6 — certification of the two sets (`g7count/certify16.py`, 2 min)

    python3 certify16.py complete_result.json      # -> certify16.json

For each of the two sixteen-point sets found, its colour pattern (121 equations in 36 variables, seven colours) is solved by msolve real root isolation at 120 bits. **Expected output:** each system is zero-dimensional with exactly 4 real solutions (reflection × the sign of x₁), all congruent to the set found. So each pattern is rigid with a unique realization, and the two realizations are the Erdős–Fishburn sets, in triangular coordinates (a, b) ↦ (a + b/2, b√3/2) with rows b = −2, …, 2:

* E₁ (symmetry group of order 2): a ∈ {0, 1}, {−1, 0, 1, 2}, {−2, −1, 0, 1}, {−2, −1, 0, 1}, {−2, −1};
* E₂ (order 6): a ∈ {1, 2}, {−1, 0, 1, 2}, {−2, −1, 0, 1}, {−2, −1, 0, 1}, {−1, 0}.

Both have squared distances {1, 3, 4, 7, 9, 12, 13}; neither contains a thirteen-point six-distance set.

**Optional fast pass** (`cores11.py`, `fast11.py`, 18 s): the same extension from the 68 eleven-point sets with ≤ 6 distances of Amit–Astra's Lemma 6 (52 lattice, 16 circular), by five points, gives the same two sets. It covers every X with m ≤ 11 and G_D ≠ C₁₁ by the table of §1, and is a check on the complete pass rather than part of the proof.

## 9. What is exact, what is numerical, and what is assumed

* **Exact.** The pattern lists (Gröbner bases modulo a prime, each discard confirmable modulo a second prime; the agreement with the published counts at ten levels is a strong check on the implementation). The rigidity certificates (Singular over ℚ). The real solutions (msolve: rational univariate representation and certified isolating intervals; msolve's modular computation is verified probabilistically, as in all current Gröbner engines — a reader who wants more can re-solve the 228 Gram systems with a second tool, they are tiny). The PSD/rank test (an exact criterion, evaluated with eigenvalues separated from zero by more than 0.7).
* **Numerical, but safe in the only direction that matters.** The extension runs in 100-digit arithmetic and treats two squared distances (scaled to minimum 1) as equal when they differ by less than 10⁻⁴⁰. Coordinates are known to about 10⁻⁹⁵ (from isolating intervals of width 2⁻⁴⁰⁰ and well-conditioned eigendecompositions), so two quantities that are exactly equal are never judged unequal: no genuine extension is missed. Two quantities that are exactly different but agree to 40 digits could only *add* candidate points or spuriously merge two new values, i.e. produce extra sixteen-point candidates; every candidate found is then certified exactly in step 6 (and in fact all 826 found are congruent to E₁ or E₂). Spurious coincidences of this size do not occur among algebraic numbers of the degrees and heights appearing here, but the argument does not need that.
* **Cited.** Altman (1963); Erdős–Fishburn (1996) / Fishburn (1995) for convex n-gons with ⌊n/2⌋ distances; Amit–Astra (2026) Lemma 4 for the thirteen-cycle; Szöllősi–Östergård (2020) Proposition 4.2 for the candidate criterion (their published counts are used only as a check on the re-implementation, not as an input).

## 10. File manifest

`g7fin/gen/`: `gen.py` (generator), `patlib.py` (nauty canonical forms), `realize.py` (polynomial systems, msolve drivers), `validate.py`, `levels/L{7..13}.json`, `levels/summary.json`, `levels/gen.log`, `README.md` (definitions, validation at the lower levels). `g7fin/`: `extract_d.py`, `pat_rigid_d.py`, `patterns_10d6.json`, `pat_rigid_10d6.json`. `g7count/`: `gram_util.py`, `gram_n.py`, `gram_to_real10.py`, `real10g/` (the 228 realization files), `extend.py`, `complete10.py`, `collect.py`, `complete_result.json`, `polygons.py`, `certify16.py`, `certify16.json`, `crosscheck.py`, `realize_all.py` (coordinate-route cross-check), `cores11.py`, `fast11.py`, `fast11_result.json`, `REPORT.md` (the full account of the computation as it was first run, with timings).
