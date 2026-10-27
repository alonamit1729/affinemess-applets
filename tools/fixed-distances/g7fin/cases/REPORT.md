# The three open cases of F(7): what can be proved without the 10-point classification

(Report of the case-analysis agent, 26 Oct 2026; scripts and certificates in this directory.)

Setting as in FINITENESS.md: X(t) a family of 16-point 7-distance sets with constant colouring, D diameter, X_D the m diameter endpoints, G_D diameter graph. Pinning lemma: a rigid subset of ≥ 6 points pins X; every 11-point ≤6-distance set is rigid (certified). Open cases: (a) m=11, G_D=C₁₁; (b) m=12, α=6; (c) m=13, G_D≠C₁₃, α=7.

**Summary.** Case (a) is CLOSED by a new computer-assisted theorem (Theorem A: no 11-point planar set has diameter graph C₁₁ and exactly seven distances). Theorem A also kills every candidate G_D in (b),(c) containing an 11-cycle. Admissible G_D in (b),(c) are enumerated exactly (37 and 175 remain). Still open: forests (33 in (b), 158 in (c)) and suns on C₉/C₇ (4 in (b), 17 in (c)).

## 1. Two structural lemmas for diameter graphs (proved)

**Leaf lemma.** If v ∈ X_D has diameter neighbours a,b,c with b angularly between a and c (neighbours of v lie on a 60° arc of the circle of radius D about v), then b is a leaf of G_D. Proof: |ab|,|bc|<D so another neighbour w of b is distinct from v,a,c; segment bw is vertex-disjoint from va, vc, hence crosses both; but b is interior to the convex cone at v spanned by a,c, and a segment from an interior point of a convex cone leaves through at most one boundary point unless it passes through v (then |bw|>D). ∎ Hence G_D minus leaves has max degree 2: every component is a caterpillar or an odd cycle with pendant leaves (a "sun"); at most one cycle (Wei).

**Cycle lemma.** If G_D contains a cycle C, every edge of G_D meets V(C), so a unicyclic G_D is a connected sun. Proof: an edge xy disjoint from C must cross all 2k+1 edges of C; x,y are extreme points of conv X not in V(C), so outside conv V(C) where the closed polygon C has winding number 0; the interior of xy contains no vertex of C; so the crossing number is even. ∎

Sanity check: with these filters `enum_GD.py 11 5` returns only C₁₁ (Wei's lemma); without the cycle lemma it would also return C₃ ∪ 4K₂ etc.

**Enumeration** (`enum_GD.py`, filters W+L+C, then α; exhaustive, isomorphism-free via networkx): m=12, α=6: 38 graphs (33 caterpillar forests with a perfect matching; suns on C₇ (1), C₉ (3), C₁₁ (1)). m=13, α=7: 181 graphs (158 forests with a 6-matching; suns on C₇ (4), C₉ (13), C₁₁ (6)). Files `GD_m12_alpha6_filtered.json`, `GD_m13_alpha7_filtered.json`. Number of maximum independent sets per graph: 3–64 (m=12), 1–32 (m=13).

## 2. X_D carries all seven classes (proved)

Every 11-subset of X uses all seven classes (else rigid ⇒ pinned). Taking 11 hull vertices: in all three cases X_D contains all seven classes (and every class restricted to X_D has vertex cover ≥ 1,2,3 for m=11,12,13). So X_D is never R₁₁, R₁₂, R₁₃∖v, R₁₃: X_D is a convex m-gon with EXACTLY seven distances; and any 11-point subset Y ⊂ X either has all seven classes or is rigid.

## 3. Theorem A: no 11-point set has diameter graph C₁₁ and exactly seven distances

Y such a set, p₀..p₁₀ convex order, D=1. As in §6.1 of the paper with 13→11: C₁₁ is the step-5 star; the bisector of side p_j p_{j+1} passes through p_{j+6}; hence (4′): |p_i p_{i+1}| < |p_i p_{i+2}| < |p_i p_{i+3}| < |p_i p_{i+4}| < |p_i p_{i+5}| = 1 for all i, both directions. With classes e₁<…<e₆<1, the step-k chord lies in {e_k,e_{k+1},e_{k+2}}; sides take ≤ 3 lengths. The star s_j=p_{5j} is an equilateral closed 11-gon with angle θ_j = 2arcsin(s/2) < π/3 (s the subtended convex side), Σθ=π (rotation index 5). Star r-chords are convex chords: r=2↔step1, r=3↔step4, r=4↔step2, r=5↔step3, squared lengths f(θ₁,θ₂)=3−2cosθ₁−2cosθ₂+2cos(θ₁+θ₂), g(θ₁,θ₂,θ₃)=4−2Σcosθ_i+2cos(θ₁+θ₂)+2cos(θ₂+θ₃)−2cos(θ₁+θ₂+θ₃). Proved: on (0,π/3)², f strictly decreasing in each angle (∂f/∂θ₁ = 2(sinθ₁−sin(θ₁+θ₂))<0); on (0,π/3)³, g strictly decreasing in the middle angle (∂g/∂θ₂ = −4 sin((θ₁+2θ₂+θ₃)/2) sin(θ₁/2) sin(θ₃/2)) and strictly increasing in the outer ones (∂g/∂θ₁ = 2[sinθ₁ + 2sin(θ₃/2)cos(θ₁+θ₂+θ₃/2)] > 0: if the cosine is negative, 2sin(θ₃/2)|cos| < sin(θ₁+θ₂+θ₃/2−π/2) < sinθ₁). The 5-chord is not monotone on the box (checked: sign-definite only when the 4 angles sum ≤ π/2), so only equalities are used for it.

*One side length:* star regular ⇒ Y=R₁₁ ⇒ 5 distances. Excluded.

*Two side lengths:* §6.2–6.4 carry over verbatim (exponents νk−11B_k, F_ν=(z¹¹−ω)/(z−η), r=11⁻¹ mod 2ν exists for all ν≤10). `certify_cycle_odd.py 11 --all-words` (the paper's integer-interval verifier, parametrised in N) proves no nonconstant word closes: 2046 words, 2074 intervals, depth ≤ 2 (`cycle11_allwords_certificate.json`). Validation: for N=13 the script reproduces the published certificate exactly (378/393/depth 3 and 8190/8459/depth 3). PROVED, no distance count needed.

*Three side lengths, seven classes:* `words11.py` treats the cyclic side word (α<β<γ; 8043 words up to dihedral symmetry with all three letters) as a CSP on the 44 chord classes: windows, the 22 chains (4′), same angle-type ⇒ same class, strict inequalities from monotonicity of f,g, and all six classes e₁..e₆ present (§2). Exactly TWO words survive: α⁸βγβ and (αβ)⁴γγβ (both mirror-symmetric, so closure alone leaves a 1-parameter family of closed stars). Each forces ≥ 7 chord coincidences common to all its class assignments (`words11_detail.py`), e.g. for α⁸βγβ: side(β)=g(α,α,α), side(γ)=g(α,α,β)=h(α,α,β,γ)=h(α,β,γ,β), f(β,γ)=g(α,β,γ)=h(α,α,α,β), f(α,β)=h(α,α,α,α), with 8α+2β+γ=π. `words11_certify.py`: 2-D outward-rounded interval branch-and-bound (mpmath iv, 60 bits) over the closed domain 0≤α≤β≤γ≤π/3: for both words every box is disjoint from the domain or has a residual interval excluding 0 (277 and 344 boxes, depth ≤ 19; `words11_certificate.json`). Numerically the residual sum stays ≥ 0.008 (`words11_numeric.py`, heuristic). PROVED modulo mpmath's interval sin/cos (porting to the paper's integer intervals is routine: the equations are finite sums of cosines of integer combinations of α,β,γ).

**Corollary (case (a) closed).** m=11, G_D=C₁₁ is impossible in a family: X_D has ≤6 classes (rigid, pinned) or exactly 7 (Theorem A).

**Corollary (C₁₁-suns).** If G_D contains an 11-cycle (C₁₁+leaf in (b), C₁₁+2 leaves in (c)), the 11 cycle vertices Y are in convex position with G_D[Y]=C₁₁ and ≤7 classes: Y is rigid or contradicts Theorem A. Excluded: 1 graph in (b), 6 in (c).

## 4. What remains and how far the tools reach

Remaining G_D: (b) 33 caterpillar forests with a perfect matching, 3 C₉-suns, 1 C₇-sun; (c) 158 forests with a 6-matching, 13 C₉-suns, 4 C₇-suns (`GD_m12_remaining.json`, `GD_m13_remaining.json`).

**C₉/C₇-suns (reduced, open).** For the cycle vertices the chain has 3 (resp. 2) steps with windows of width 4 (5), and "all classes present" is unavailable (9-point ≤6-distance sets not known rigid). `words_odd.py 9` leaves 12 two-letter, 8 three-letter, 1 four-letter words; N=7 the CSP is weak (16/108/67/2). Two-letter words are rigid whenever the closure polynomial is not identically zero (isolated zeros): `closure_nonzero.py` verifies this exactly for all words except the six period-3 words of N=9 (threefold-symmetric stars close for every θ). Several three-letter N=9 words (α⁶βγβ, (αβ)³γγβ) force no coincidence, so no nonexistence certificate; an isolation (interval-Newton) certificate on closure + assignment equalities would be needed — finite and mechanical, not done.

**Forests (open).** General chain lemma (proved, same bisector argument): for side p_j p_{j+1} with diameter neighbours a of p_j, b of p_{j+1}, the bisector exits the polygon between a and b, so vertices on the arc b→p_j are closer to p_j and those on p_{j+1}→a closer to p_{j+1}. For the perfect matching p_i p_{i+6} (endpoints of 6 pairwise crossing segments are antipodal in convex order) this gives the full chain |p_i p_{i+1}|<…<|p_i p_{i+5}|<1, hence windows step k ∈ {e_k,e_{k+1}}; 87 632 random samples agree (`matching_chain_test.py`, heuristic). Missing: an angle-type parametrisation — the matching 12-gon has 15 dof and no equilateral star, so no "same type ⇒ same length" equalities and no monotonicity; a direct algebraic certificate is 15-dimensional. Practical route stays the 10-point input: for every maximum independent set I, I ∪ (X∖X_D) is a non-rigid 10-set realising one of the three unknown SÖ patterns, as is every 10-subset missing a class with τ=6.

**Idea 1** (τ(G_c) ≥ 6, planarity of min-distance graph) gave no contradiction by itself; its useful consequence is §2. **Idea 4** reduces to τ(G_c) ≤ 5 for some class, which nothing in (b)/(c) forces.

## 5. Files
certify_cycle_odd.py (+ cycle11/13 certificates), enum_GD.py (+ GD_*.json), chord_monotone.py, words11.py, words11_detail.py, words11_numeric.py (heuristic), words11_certify.py (+ words11_certificate.json), words_odd.py (N=7,9,11), words_detail.py, closure_nonzero.py, matching_chain_test.py (heuristic).
