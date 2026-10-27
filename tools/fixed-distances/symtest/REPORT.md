# Symmetry of positive-dimensional families of n-point d-distance sets (e(n,d) <= 0)

Data: `applet_data.json` (all family items, n = 3..16, d <= 7), `sat54..sat57.json` (saturated ideals of the
positive-dimensional types for n = 5), `cert5d.json`. Tools in this folder: `symmetry.py` (automorphism group of a
point set by backtracking on the distance matrix, relative tolerance 1e-6, classified as trivial / mirror / C_k / D_k),
`step1_catalog.py`, `step2_sample.py` (random linear slices of the saturated ideals through msolve, exactly as in
`strat.py`), `famcheck.py` (numerical tangent dimension and predictor-corrector trace along a family),
`msolve_check.py` (independent exact check of any catalog family: colouring system + x1 = 1 + random slices),
`superset*.py`, `subsets.py`, `chains.py` (structural analysis of the counterexamples).

e(n,d) = 2n - 4 - C(n,2) + d.

## 1. Catalog samples (one sample configuration per family item)

| cell (n,d) | e(n,d) | families (by dim) | sample has nontrivial symmetry | symmetry types of samples |
|---|---|---|---|---|
| (3,2) | 1 | 1 (1×dim1) | 1 | mirror:1 |
| (3,3) | 2 | 1 (1×dim2) | 0 | trivial:1 |
| (4,3) | 1 | 9 (9×dim1) | 6 | mirror:4, trivial:3, D_2:2 |
| (4,4) | 2 | 7 (7×dim2) | 2 | trivial:5, mirror:1, C_2:1 |
| (4,5) | 3 | 2 (2×dim3) | 0 | trivial:2 |
| (4,6) | 4 | 1 (1×dim4) | 0 | trivial:1 |
| (5,4) | 0 | 25 (25×dim1) | 19 | mirror:18, trivial:6, D_2:1 |
| (5,5) | 1 | 410 (402×dim1, 8×dim2) | 37 | trivial:373, mirror:33, C_2:3, D_2:1 |
| (5,6) | 2 | 245 (244×dim2, 1×dim3) | 8 | trivial:237, mirror:7, C_2:1 |
| (5,7) | 3 | 76 (76×dim3) | 1 | trivial:75, mirror:1 |
| (6,4) | -3 | 4 (4×dim1) | 3 | D_3:2, trivial:1, C_3:1 |
| (6,5) | -2 | 6 (5×dim1, 1×dim2) | 6 | D_2:4, C_2:1, C_3:1 |
| (6,6) | -1 | 10 (9×dim1, 1×dim2) | 8 | D_2:8, trivial:2 |
| (6,7) | 0 | 113 (8×dim1, 105×dim2) | 108 | mirror:92, C_2:16, trivial:5 |
| (7,5) | -6 | 2 (2×dim1) | 2 | mirror:1, D_3:1 |
| (7,6) | -5 | 7 (7×dim1) | 7 | C_3:3, D_2:2, mirror:1, D_3:1 |
| (7,7) | -4 | 10 (9×dim1, 1×dim2) | 9 | D_2:8, trivial:1, C_3:1 |
| (8,6) | -10 | 2 (2×dim1) | 2 | D_4:2 |
| (8,7) | -9 | 4 (4×dim1) | 4 | C_4:2, mirror:1, D_4:1 |
| (9,6) | -16 | 1 (1×dim1) | 1 | D_3:1 |
| (9,7) | -15 | 8 (8×dim1) | 8 | D_3:7, D_4:1 |
| (10,7) | -22 | 2 (2×dim1) | 2 | D_5:2 |

Symmetry types count orientation-reversing elements: "mirror" = one reflection (D_1), C_k = rotations only,
D_k = k reflections + k rotations. The (6,7) column counts 105 "dim 2" items; the catalog's dims are not always
the family's true dimension (see section 4).

Verification that the sample's symmetry is the generic symmetry of its family (all families with e <= 0, n >= 6):
(a) numerical trace along the family from the sample in both tangent directions (`step3_trace.py`): the symmetry type
is constant along every trace (169 families); (b) exact random slices through msolve (`checks_eneg.log`,
`checks_67.log`): every non-degenerate real point found has the catalog sample's symmetry (or, for 13 types that
contain two components, mirror and C_2, both symmetric). 17 families had no non-degenerate real point on the slices
(positive-dimensional degenerate components or no real intersection); for those (a) is the evidence.

## 2. n = 5: generic symmetry from random slices of the saturated ideals (`samples_5d.json`)

4 random slices per type, every real non-degenerate solution classified. No type had samples of mixed symmetry
(trivial and nontrivial), i.e. the symmetry is a property of the family, not of the sample.

| cell | e | positive-dimensional types (sat5d) | generically symmetric | generically trivial |
|---|---|---|---|---|
| (5,4) | 0 | 25 (all dim 1) | 19: 18 mirror, 1 D_2 | 6 (types 62, 84, 92, 245, 325, 326 = catalog #473, #474, #475, #477, #491, #492) |
| (5,5) | 1 | 407 (399 dim 1, 8 dim 2) | dim 1: 27 mirror + 1 D_2; dim 2: 6 mirror | dim 1: 371; dim 2: 2 |
| (5,6) | 2 | 244 (243 dim 2, 1 dim 3) | 7 mirror (6 dim 2, 1 dim 3) | 237 |
| (5,7) | 3 | 76 (dim 3) | 1 mirror | 75 |

In terms of the catalog's 410 (5,5) families: 37 are generically symmetric (33 mirror, 1 D_2, and the 3 C_2 families
#197, #220, #238, which are 1-dimensional components living inside types 226, 250, 270 whose saturated ideals are
2-dimensional, so 2-plane slices never see them; their symmetry was confirmed by tracing) and 373 are generically
trivial. The (5,6) C_2 family #91 is likewise a 2-dim component inside the one 3-dim type. Every catalog sample has
exactly its family's generic symmetry (no catalog sample is "accidentally" symmetric).

## 3. The conjecture

**As stated, the conjecture is false**, in both forms (every member / generic member): there are positive-dimensional
families with e(n,d) <= 0 whose generic member has trivial symmetry, in every cell with e <= 0 that was not populated
exclusively by symmetric constructions:

| cell | e | families | generically symmetric | generically trivial |
|---|---|---|---|---|
| (5,4) | 0 | 25 (complete) | 19 | **6** |
| (6,4) | -3 | 4 | 3 | **1** (#23) |
| (6,5) | -2 | 6 | 6 | 0 |
| (6,6) | -1 | 10 | 8 | **2** (#3, #4) |
| (6,7) | 0 | 113 | 108 | **5** (#0, #1, #2, #3, #8) |
| (7,5) | -6 | 2 | 2 | 0 |
| (7,6) | -5 | 7 | 7 | 0 |
| (7,7) | -4 | 10 | 9 | **1** (#4) |
| (8,6), (8,7), (9,6), (9,7), (10,7) | -10 .. -22 | 17 | 17 | 0 |

Caveat on the cells with 0 counterexamples: the catalogs for n >= 7 are lower bounds found mostly by *symmetric*
searches (their labels read "3-fold rotational symmetry: ...", "dihedral symmetry D3: ..."), and every counterexample
above came from the few non-symmetric seeds (regular-polygon subsets, lattice pieces). So the absence of asymmetric
families there is selection, not evidence. Indeed asymmetric families with e <= 0 exist for every n >= 6 by an
explicit construction (verified numerically for n = 6..12, `chains.py`): the **circular equiangular chain** - points
R·e^(ik·theta), k in a set S of n integers, with theta a free parameter. Its distances are the chords 2R|sin(m·theta/2)|,
m in S - S, so d = |(S - S) ∩ Z_{>0}|, and it is a 1-parameter family for every S; for generic theta its symmetry group
is trivial unless S is symmetric under reversal (S = c - S). With S = {0,1,...,n-2,n}: d = n and e(n,n) = 3n - 4 - C(n,2)
= -1, -4, -8, -13, -19, -26, -34 for n = 6..12. The catalog counterexamples (6,6)#3, (6,7)#0, (7,7)#4 are exactly such
chains (S = {0,1,2,3,4,6}, {0,1,2,3,4,7}, {0,1,2,3,4,5,7}).

**What the data does support (sharpest form found).** Every one of the 15 counterexamples is a sub-configuration of a
positive-dimensional family of configurations with a fixed nontrivial symmetry group: deleting 1-3 points (orbit
points, or the point carrying an extra distance) from a symmetric family with d' >= d distances and e(n',d') < e(n,d)
gives the asymmetric family. The symmetric super-families were found by solving for extra points whose distances to
the configuration lie in its distance set (allowing at most two new values) and verified as genuine 1-parameter
symmetric families (tangent dimension 1, symmetry constant along numerical traces, and for n' <= 8 exact msolve
slices). So the conjecture that survives is:

> Every positive-dimensional family of n-point d-distance sets with e(n,d) <= 0 is the restriction (deletion of
> points) of a positive-dimensional family of configurations with a nontrivial symmetry group that is constant along
> the family; the dependency among the equal-distance equations is inherited from that symmetry. Equivalently: the
> *minimal* (point-deletion-closed) families with e <= 0 are symmetric.

The symmetric super-families have e' <= -2 in every case, and the subset families themselves have e in {0, -1, -3, -4}.
Also supported: for the complete cells with n <= 5, symmetric families are the majority only at e = 0 ((5,4): 19 of 25)
and a small minority once e >= 1 ((5,5): 37 of 410; (5,6): 8 of 245; (5,7): 1 of 76).

## 4. Counterexamples (e <= 0, generic member has trivial symmetry)

All 15 verified three ways: catalog sample has trivial symmetry with tolerance 1e-6 and also 1e-3; numerical trace
along the family keeps trivial symmetry up to the family's boundary in both directions; exact msolve solve of the
colouring system with x1 = 1 and a random linear slice is zero-dimensional with real non-degenerate solutions, all of
trivial symmetry (for (5,4): 6-14 samples per family from the saturated ideals, all trivial).

| family | e | dim (verified) | generic symmetry | structure: sub-configuration of a symmetric family |
|---|---|---|---|---|
| (5,4) #473 | 0 | 1 | trivial | (6,4) #24, C_3 family minus 1 orbit point (same 4 distances) |
| (5,4) #474 | 0 | 1 | trivial | new (6,5) mirror family (e=-2) minus the point carrying the 5th distance |
| (5,4) #475 | 0 | 1 | trivial | new (6,5) mirror family (e=-2) minus 1 point |
| (5,4) #477 | 0 | 1 | trivial | new (6,5) mirror family (e=-2) minus 1 point |
| (5,4) #491 | 0 | 1 | trivial | (6,4) #26, D_3 family minus 1 point (same 4 distances) |
| (5,4) #492 | 0 | 1 | trivial | (6,4) #24, C_3 family minus 1 orbit point (same 4 distances) |
| (6,4) #23 | -3 | 1 | trivial | new (7,5) mirror family (e=-6) minus 1 point |
| (6,6) #3 | -1 | 1 | trivial | (7,6) #58 mirror chain {0..6} minus vertex 5; circular chain S={0,1,2,3,4,6} (points at angles k*theta) |
| (6,6) #4 | -1 | 1 | trivial | new (8,6) mirror family (e=-10) minus 2 points |
| (6,7) #0 | 0 | 1 | trivial | (8,7) #4 mirror chain {0..7} minus vertices 5,6; circular chain S={0,1,2,3,4,7} |
| (6,7) #1 | 0 | 1 | trivial | (7,8) mirror family (e=-3) minus 1 point; 5-chain {0..4} + 1 point |
| (6,7) #2 | 0 | 1 | trivial | (7,8) mirror family (e=-3) minus 1 point; 5-chain {0..4} + 1 point |
| (6,7) #3 | 0 | 1 | trivial | (8,9) mirror family (e=-7) minus 2 points; 4-chain + 2 points |
| (6,7) #8 | 0 | 1 | trivial | (9,9) C_3 family (e=-13) minus 3 points |
| (7,7) #4 | -4 | 1 | trivial | (8,7) #4 mirror chain {0..7} minus vertex 6; circular chain S={0,1,2,3,4,5,7} |

Coordinates (one generic member of each family, shortest distance = 1, point 0 at the origin; `cls` = distance class of each pair in the order (0,1),(0,2),...,(n-2,n-1), classes numbered by increasing length):

- **(5,4) #473**: points (0.000000, 0.000000), (-0.866025, -0.500000), (-0.866025, 0.500000), (0.712069, -0.702110), (-0.614015, 1.467725); distances [1.0, 1.59098, 1.9838, 2.54297]; cls [1, 1, 1, 2, 1, 2, 3, 3, 1, 4]
- **(5,4) #474**: points (0.000000, 0.000000), (-2.075266, 1.198155), (-2.075266, -1.198155), (-1.405571, -1.940792), (-2.383560, -2.149446); distances [1.0, 2.39631, 3.20959, 3.36177]; cls [2, 2, 2, 3, 2, 3, 4, 1, 1, 1]
- **(5,4) #475**: points (0.000000, 0.000000), (-0.866025, -0.500000), (-0.866025, 0.500000), (0.888414, 0.459043), (-0.024275, -1.039867); distances [1.0, 1.04015, 1.75492, 1.99946]; cls [1, 1, 1, 2, 1, 4, 1, 3, 3, 3]
- **(5,4) #477**: points (0.000000, 0.000000), (-0.866025, -0.500000), (-0.866025, 0.500000), (1.192630, -3.065695), (3.251285, 0.500000); distances [1.0, 3.28951, 4.11731, 4.23701]; cls [1, 1, 2, 2, 1, 2, 4, 3, 3, 3]
- **(5,4) #491**: points (0.000000, 0.000000), (6.597009, 0.000000), (3.298505, 5.713178), (0.609361, -0.792893), (2.307159, 5.581902); distances [1.0, 6.03992, 6.59701, 7.03992]; cls [3, 3, 1, 2, 3, 2, 4, 4, 1, 3]
- **(5,4) #492**: points (0.000000, 0.000000), (-1.420077, 0.819882), (-1.420077, -0.819882), (-2.212598, -0.210037), (-0.495675, 1.201303); distances [1.0, 1.29955, 1.63976, 2.22255]; cls [3, 3, 4, 2, 3, 2, 1, 1, 4, 4]
- **(6,4) #23**: points (0.000000, 0.000000), (-1.960435, 1.131857), (-1.960435, -1.131857), (-1.664068, -0.176783), (-0.985133, -0.910982), (-0.678935, 0.734198); distances [1.0, 1.34178, 1.67343, 2.26372]; cls [4, 4, 3, 2, 1, 4, 2, 4, 2, 1, 1, 4, 1, 2, 3]
- **(6,6) #3**: points (0.000000, 0.000000), (2.590247, 0.000000), (0.286110, 1.183355), (1.795123, -0.921937), (1.414594, 1.640206), (2.443338, 0.989150); distances [1.0, 1.21745, 2.01803, 2.16595, 2.59025, 2.63597]; cls [5, 2, 3, 4, 6, 5, 2, 3, 1, 5, 2, 4, 5, 3, 2]
- **(6,6) #4**: points (0.000000, 0.000000), (0.000000, 1.000000), (-0.193867, -0.981028), (-2.006420, -0.803648), (-1.059892, -1.481028), (-1.059892, -0.481028); distances [1.0, 1.16394, 1.82121, 1.99049, 2.16138, 2.69794]; cls [1, 1, 5, 3, 2, 4, 6, 6, 3, 3, 1, 1, 2, 1, 1]
- **(6,7) #0**: points (0.000000, 0.000000), (1.000000, 0.000000), (1.829762, 0.558117), (2.206773, 1.484326), (2.002669, 2.463275), (-0.613401, 2.942073); distances [1.0, 1.91299, 2.65952, 3.00534, 3.17465, 3.35542, 3.41354]; cls [1, 2, 3, 5, 4, 1, 2, 3, 6, 1, 2, 7, 1, 5, 3]
- **(6,7) #1**: points (0.000000, 0.000000), (1.806584, 0.000000), (0.366231, 1.090472), (0.856375, -0.648350), (1.515164, 1.033833), (2.754550, -0.318371); distances [1.0, 1.07412, 1.15033, 1.80658, 1.83427, 1.92664, 2.77289]; cls [4, 3, 2, 5, 7, 4, 3, 2, 1, 4, 3, 7, 4, 6, 5]
- **(6,7) #2**: points (0.000000, 0.000000), (1.388540, 0.000000), (1.744987, 1.342009), (0.539451, 2.031013), (-0.435933, 1.042747), (0.694270, 1.043071); distances [1.0, 1.09241, 1.1302, 1.253, 1.38854, 2.10143, 2.20136]; cls [5, 7, 6, 3, 4, 5, 7, 6, 4, 5, 7, 2, 5, 1, 3]
- **(6,7) #3**: points (0.000000, 0.000000), (2.947693, 0.000000), (-0.054197, 0.998530), (1.751469, -1.599156), (2.306183, -0.767115), (1.861567, 1.562534); distances [1.0, 1.90294, 1.99706, 2.3717, 2.43042, 2.94769, 3.16361]; cls [6, 1, 4, 5, 5, 7, 3, 1, 2, 7, 6, 3, 1, 7, 4]
- **(6,7) #8** (generic msolve sample; the catalog sample is a special lattice point): points (0.000000, 0.000000), (2.098487, 0.000000), (2.259715, -1.765944), (-0.560722, -1.073999), (0.089026, -2.096598), (1.049244, -1.817343); distances [1.0, 1.21156, 1.77329, 2.09849, 2.19573, 2.8679, 2.90408]; cls [4, 6, 2, 4, 4, 3, 6, 7, 4, 7, 5, 2, 2, 3, 1]
- **(7,7) #4**: points (0.000000, 0.000000), (3.336998, 0.000000), (2.370523, -3.193976), (-0.406646, -1.343871), (2.168499, 0.778435), (3.454020, -2.301012), (0.771695, 0.635993); distances [1.0, 1.40405, 2.30398, 2.64296, 3.337, 3.97754, 4.15029]; cls [5, 6, 2, 3, 7, 1, 5, 6, 2, 3, 4, 5, 6, 2, 7, 5, 6, 3, 5, 2, 6]

Symmetric super-configurations (extra points added to the member above, same normalisation; verified: tangent dimension 1, symmetry constant along a numerical trace, and for n <= 8 zero-dimensional random msolve slices with real symmetric solutions):

- **(5,4) #474 + 1 point(s)** -> (6,5), e=-2, symmetry ['mirror']: extra points (0.977989, -2.187656)
- **(5,4) #475 + 1 point(s)** -> (6,5), e=-2, symmetry ['mirror']: extra points (0.841750, -0.539867)
- **(5,4) #477 + 1 point(s)** -> (6,5), e=-2, symmetry ['mirror']: extra points (-2.924681, -3.065695)
- **(6,4) #23 + 1 point(s)** -> (7,5), e=-6, symmetry ['mirror']: extra points (0.296367, 0.955074)
- **(6,6) #4 + 2 point(s)** -> (8,6), e=-10, symmetry ['mirror']: extra points (-0.946528, 0.677379), (-1.812554, 0.177379)
- **(6,7) #1 + 1 point(s)** -> (7,8), e=-3, symmetry ['mirror']: extra points (0.807138, -1.647138)
- **(6,7) #2 + 1 point(s)** -> (7,8), e=-3, symmetry ['mirror']: extra points (0.794144, 0.804173)
- **(6,7) #3 + 2 point(s)** -> (8,9), e=-7, symmetry ['mirror']: extra points (-0.054197, -0.998530), (1.751469, 1.599156)
- **(6,7) #8 + 3 point(s)** -> (9,9), e=-13, symmetry ['C_3']: extra points (0.238267, 0.971200), (1.448738, 1.022599), (2.820437, -0.691945)

## 5. Side findings about the catalog

- Dimension corrections: (6,7) #4, #5, #6, #7 are 2-dimensional families (one random slice leaves a positive-
  dimensional set, two slices give finitely many points), not 1-dimensional; (6,7) #8 is 1-dimensional (two random
  slices miss it), its catalog sample is a lattice point where two branches cross (tangent dimension 2).
- The 3 C_2 families in (5,5) and the C_2 family in (5,6) are lower-dimensional components inside higher-dimensional
  types (sat5d.json records the maximal dimension of a type), consistent with the catalog counts 402 = 399 + 3 and
  245 = 244 + 1.
- New symmetric families not in the catalog (coordinates above, `superset_families.json`): three mirror-symmetric
  1-parameter families in (6,5) (the catalog has only D_2, C_2, C_3 families there), one mirror family in (7,5)
  (its colouring type differs from #17), one mirror family in (8,6). These are the super-families of (5,4) #474,
  #475, #477, (6,4) #23 and (6,6) #4; each contains a point carrying a distance class that no other pair uses.
- Types containing two components of different symmetry (mirror and C_2) on random slices: (6,5) #27, (6,7) #5, #9,
  #12, #14, #17, #18, #19, #22, #23, #24, #27, #30 and (8,7) #4; the catalog shows one of them.

Files: `step1_families.json` (all family items with symmetry), `step1_table.md`, `samples_5{4,5,6,7}.json`,
`step2_summary.json`, `step3_trace.json`, `msolve_generic_summary.json`, `generic_samples.json`,
`superset_families.json`, `counterexamples.txt`, `chains.txt`, `subsets.txt`, `superset*.txt`, logs `*.log`.
