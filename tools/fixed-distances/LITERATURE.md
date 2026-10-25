# Literature notes for "Planar sets with a fixed number of distances"

What each paper establishes, and how the applet's cells use it. g(k) = largest planar k-distance set.

## Einhorn–Schoenberg 1966 (Indag. Math. 28)
- Two-distance sets in the plane: the six 4-point ones; R₅ is the only 5-point one. g(2) = 5.
- Cells (4,2), (5,2).

## Erdős–Fishburn 1996, *Maximum planar sets that determine k distances* (Discrete Math. 160)
- g(2)=5, g(3)=7, g(4)=9, g(5)=12; all maximum sets for k ≤ 4.
- k=3: R₇ and R₆⁺ (hexagon with centre). k=4: R₉ and three more — two pieces of the triangular lattice and "three equilateral triangles with the same centre and a horizontal edge" (their Fig. 1 top right). The last one is the set the applet first found by a mirror-symmetric search (square + apexes + axis point; distances √3−1, √2, 2, √3+1): it has 3-fold symmetry about (0, −1/√3), so the two descriptions coincide.
- k=5: the 12-point lattice piece (rows 3,4,3,2); uniqueness left open there (settled by Shinohara 2008).
- Conjectures: g(6) = 13 with three sets (R₁₃, R₁₂⁺, the 13-point lattice star); for k ≥ 7 every maximum set is a lattice piece. Lattice examples (Fig. 4): (k,n) = (7,16) ×2, (8,19), (9,21), (10,25), (11,27), (13,31). Bao–Yu 2025 point out that the (10,25) example determines 11 distances.
- Lemma 2 (convex polygons, from Altman 1963, Fishburn 1995, EF Geom. Dedicata 1996): a convex n-gon has ≥ ⌊n/2⌋ distances; M₂ₖ₊₁(k) = {R₂ₖ₊₁}; M₂ₖ(k) = {R₂ₖ, R₂ₖ₊₁−1} (k ≥ 4); M₇(4) = {R₈−1, every R₉−2}; M₉(5) = {R₁₀−1, every R₁₁−2}.
- Cells (7,3), (9,4), (16,7).

## Shinohara 2004, *Classification of three-distance sets in two dimensional Euclidean space* (European J. Combin. 25)
- Theorem 1: exactly 34 five-point three-distance sets (figs. 501–534, with a table of the distance pairs (α, β)).
- Theorem 2: no 3-distance set with > 7 points; two 7-point sets; six maximal 6-point sets (figs. 601–606); sixteen maximal 5-point sets. Twenty-four maximal sets with ≥ 5 points in all (Einhorn–Schoenberg had conjectured five).
- Check against the applet: each of our 34 sets was matched to one of figs. 501–534 by its (α, β); all 34 correspond. Apparent misprint: the table gives β = 2 sin 4π/9 for fig. 524, but the construction in the text (∠pt₁p₄ = 7π/9) gives 2 sin 7π/18 — our "18-gon: O + V₀V₄V₇V₁₁". The 9 six-point sets = 6 maximal + the three 6-subsets of R₇, R₆⁺.
- Cells (5,3), (6,3), (7,3).

## Shinohara 2008, *Uniqueness of maximum planar five-distance sets* (Discrete Math. 308)
- Theorem 1.2(a): every 8-point 4-distance set is R₈, R₇⁺, the square with four outward equilateral apexes (Fig. 1(e); distances 1, √2, √(2+√3), 1+√3), or an 8-point subset of a 9-point 4-distance set. Enumerating the subsets gives 15 sets up to congruence (applet cell (8,4)); the earlier extension search had found 14.
- Theorem 1.2(b): the 12-point 5-distance set is unique. Also diameter-graph machinery (Props 3.1–3.4) used later by Wei.
- Cells (8,4), (12,5).

## Wei 2011, *Classification of eleven-point five-distance sets in the plane* (Ars Combin. 102)
- Theorem 13: exactly four 11-point 5-distance sets: R₁₁ and three lattice pieces (all 11-subsets of the 12-point set). Re-derived in the applet by a clique search in the lattice with distance set {1, √3, 2, √7, 3} (`lattice_11_5.json`).
- Cell (11,5).

## Wei 2012, *A proof of Erdős–Fishburn's conjecture for g(6) = 13* (Electron. J. Combin. 19(4) P38); Szöllősi–Östergård 2020, *Constructions of maximum few-distance sets in Euclidean spaces* (Electron. J. Combin. 27(1) P1.23)
- g(6) = 13 with exactly the three conjectured sets. Szöllősi–Östergård classify planar 6-distance sets by computer (isomorph-free generation + Gröbner bases + cliquer).
- Cell (13,6).

## Lan–Wei 2013, *Classification of seven-point four-distance sets in the plane* (Math. Notes 93, 510–522)
- Theorem 8: exactly 42 seven-point four-distance sets, all drawn in Fig. 4 (701–742). Structure given in the proof: R₉−2 (four versions: 701–703, 708), R₈−1 (707), R₇⁺−1 (704), nine 7-subsets of the three-concentric-triangles 9-set (705, 709–716), twenty pieces of the triangular lattice (706, 721–739), two 7-subsets of the square-with-four-apexes 8-set (717, 718), two 7-subsets of the 10-point 5-distance set "R₅ doubled about its centre" (719, 720; it is a regular pentagon with the inner pentagon of its pentagram, ratio φ²), 740 (square-lattice piece, 1:√2:2:√5), 741 (1:√(2+√3):√(4+2√3):√(5+2√3)), 742 (1:√2:√(2+√3):√(4+√3)). Five maximal: 719, 720, 740, 741, 742.
- Method: diameter graphs (Lemmas 6, 7 classify the convex hull cases) plus extension of the 34 five-point three-distance sets as in Shinohara 2008.
- Applet: all 42 reconstructed from these descriptions (`seven4.py`), each verified to have exactly four distances, exact coordinates throughout; cell (7,4) complete (†).

## Bannai–Bannai–Stanton 1983 / Blokhuis 1984
- g(k) ≤ C(k+2, 2): 36 points for k = 7. Row "17–36, 7".

## Budianto 2022; Bao–Yu 2025 (arXiv:2509.00880)
- Partial results toward g(7) = 16; lattice constructions for m ≤ 13 and beyond (16, 19, 21, 24, 27, 27, 31 points for m = 7..13).

## Still wanted
- Wei 2012 (figures), Fishburn 1995 and Erdős–Fishburn Geom. Dedicata 1996 (convex polygons), Budianto 2022.
- Nothing published seems to count the configurations for the non-maximal cells with d ≥ 4 and 6 ≤ n ≤ 10 (other than (7,4), (8,4)); those cells stay as bounds.
