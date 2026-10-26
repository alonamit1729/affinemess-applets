import json
F=json.load(open('summary/figs.json')); F34=json.load(open('summary/figs34.json'))
def esc(s): return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
fig3='''<div class="tw"><table>
<tr><th>from \\ to</th><th>(5,6) families</th><th>(5,5) families</th><th>(5,4) families</th><th>(5,4) rigid</th><th>(5,3) rigid</th><th>(5,2)</th><th>collisions</th></tr>
<tr><td>(5,7) 76 three-folds (sampled)</td><td class="n">1152</td><td class="n">96</td><td class="n">11</td><td class="n">20</td><td class="n">6</td><td class="n">0</td><td class="n">70</td></tr>
<tr><td>(5,6) 245 families (81 exact, rest sampled)</td><td class="n">—</td><td class="n">2887</td><td class="n">314</td><td class="n">202</td><td class="n">111</td><td class="n">0</td><td class="n">224</td></tr>
<tr><td>(5,5) 410 families (exact)</td><td class="n">—</td><td class="n">—</td><td class="n">94</td><td class="n">2239</td><td class="n">596</td><td class="n">12</td><td class="n">376</td></tr>
<tr><td>(5,4) 25 curves (exact)</td><td class="n">—</td><td class="n">—</td><td class="n">—</td><td class="n">—</td><td class="n">83</td><td class="n">7</td><td class="n">25</td></tr>
</table></div>'''
leg1='\n'.join(f'<li><span class="k">{esc(k)}</span>{esc(v)}</li>' for k,v in F['legend1'])
leg4='\n'.join(f'<li><span class="sw" style="background:{c}"></span>{esc(n)}</li>' for c,n in F34['legend4'])
leg4v='\n'.join(f'<li><span class="k">{esc(k)}</span>{esc(v)}</li>' for k,v in F34['legend4v'])
leg2='\n'.join(f'<li><span class="sw" style="background:{c}"></span>{esc(n)}</li>' for n,c in F['legend2'])
html=f'''<title>Five-point shape strata</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo+Narrow:wght@500;600&family=IBM+Plex+Sans:ital,wght@0,400;0,500;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* layout: one reading column, figures breaking a little wider; tokens first */
:root{{ --bg:#f3f5f6; --fg:#1d2428; --fg-soft:#5b6670; --fig-bg:#ffffff; --node-bg:#ffffff; --node-edge:#9aa7b0; --glyph-line:#3a4650; --glyph-dot:#1d2428;
  --arc:#1f6f78; --end-deg:#b5542b; --end-inf:#8a3fa0; --grid:#dde3e6; --line:#d3dade; --accent:#1f6f78;
  --bar-a:#1f6f78; --bar-b:#b5542b; --bar-c:#8a3fa0; --bar-d:#d9a441; --region:#e8eff1;
  --display:"Archivo Narrow","Arial Narrow",system-ui,sans-serif; --body:"IBM Plex Sans",system-ui,sans-serif; --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace; }}
@media (prefers-color-scheme: dark){{ :root:not([data-theme="light"]){{ --bg:#15191c; --fg:#e4e9ec; --fg-soft:#a3adb5; --fig-bg:#1d2327; --node-bg:#1d2327; --node-edge:#5f6c75; --glyph-line:#c6ced4; --glyph-dot:#f0f3f5;
  --arc:#5fb8c2; --end-deg:#e0845b; --end-inf:#c58ad6; --grid:#2d353b; --line:#323b41; --accent:#5fb8c2; --bar-a:#5fb8c2; --bar-b:#e0845b; --bar-c:#c58ad6; --bar-d:#e6c06a; --region:#243036; color-scheme:dark }} }}
:root[data-theme="dark"]{{ --bg:#15191c; --fg:#e4e9ec; --fg-soft:#a3adb5; --fig-bg:#1d2327; --node-bg:#1d2327; --node-edge:#5f6c75; --glyph-line:#c6ced4; --glyph-dot:#f0f3f5;
  --arc:#5fb8c2; --end-deg:#e0845b; --end-inf:#c58ad6; --grid:#2d353b; --line:#323b41; --accent:#5fb8c2; --bar-a:#5fb8c2; --bar-b:#e0845b; --bar-c:#c58ad6; --bar-d:#e6c06a; --region:#243036; color-scheme:dark }}
body{{ background:var(--bg); color:var(--fg); font-family:var(--body); font-size:15.5px; line-height:1.55; padding-block:28px 60px; padding-left:max(16px,env(safe-area-inset-left,0px)); padding-right:max(16px,env(safe-area-inset-right,0px)); }}
main{{ max-width:980px; margin:0 auto; display:grid; gap:34px; }}
header h1{{ font-family:var(--display); font-weight:600; font-size:2.3rem; line-height:1.05; margin:0 0 8px; text-wrap:balance; letter-spacing:.005em; }}
header p.lede{{ max-width:62ch; margin:0; color:var(--fg-soft); font-size:1.02rem; }}
section{{ display:grid; gap:14px; }}
section > h2{{ font-family:var(--display); font-weight:600; font-size:1.45rem; margin:0; padding-top:14px; border-top:1px solid var(--line); text-wrap:balance; }}
section p{{ max-width:66ch; margin:0; }}
.num{{ font-variant-numeric:tabular-nums; }}
figure{{ margin:0; display:grid; gap:8px; }}
figure svg{{ width:100%; height:auto; border:1px solid var(--line); border-radius:3px; }}
figcaption, p.figcap{{ font-size:.9rem; color:var(--fg-soft); max-width:70ch; }}
.stats{{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:10px; }}
.stat{{ border-top:2px solid var(--accent); padding-top:6px; }}
.stat b{{ display:block; font-family:var(--mono); font-size:1.35rem; font-weight:500; font-variant-numeric:tabular-nums; }}
.stat span{{ font-size:.85rem; color:var(--fg-soft); }}
table{{ border-collapse:collapse; width:100%; font-size:.92rem; }}
.tw{{ overflow-x:auto; }}
th,td{{ text-align:left; padding:5px 9px; border-bottom:1px solid var(--line); vertical-align:top; }}
th{{ font-weight:500; color:var(--fg-soft); font-size:.8rem; letter-spacing:.05em; text-transform:uppercase; }}
td.n{{ font-family:var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap; }}
ul.legend{{ list-style:none; margin:0; padding:0; columns:2; column-gap:28px; font-size:.86rem; color:var(--fg-soft); }}
ul.legend li{{ break-inside:avoid; padding:1px 0; }}
ul.legend .k{{ display:inline-block; min-width:2.4em; font-family:var(--mono); color:var(--fg); }}
ul.legend .sw{{ display:inline-block; width:11px; height:11px; margin-right:7px; vertical-align:-1px; border-radius:2px; }}
.caveat{{ border-left:3px solid var(--end-deg); padding-left:12px; color:var(--fg-soft); font-size:.95rem; }}
code{{ font-family:var(--mono); font-size:.88em; }}
@media (max-width:560px){{ ul.legend{{ columns:1; }} header h1{{ font-size:1.9rem; }} }}
</style>
<main>
<header>
<h1>Five-point shape strata</h1>
<p class="lede">What the solution sets behind the five-point row of the fixed-distances table look like: the ambient shape space, the incidence graphs of the one-parameter families, the degeneration poset of all strata, and how each number was computed. Companion to <a href="https://affinemess.com/fixed-distances.html">Planar sets with a fixed number of distances</a>; data and scripts in <code>tools/fixed-distances/</code> of the applets repository.</p>
</header>

<section>
<h2>The ambient space</h2>
<p>Five points in the plane up to similarity have 2·5 − 4 = 6 degrees of freedom. The natural compact model is Kendall's shape space: labelled configurations with coincidences allowed (but not all five points equal), up to orientation-preserving similarity, form ℂP³; the collision loci zᵢ = zⱼ are its ten hyperplanes. Removing them leaves M₀,₆, the moduli space of six distinct points on the sphere with one at ∞. Quotienting by relabelling (S₅) and by reflections (complex conjugation) gives the unlabelled space the table is about. Rational Betti numbers:</p>
<div class="tw"><table>
<tr><th>space</th><th>b₀</th><th>b₁</th><th>b₂</th><th>b₃</th><th>b₄</th><th>b₅</th><th>b₆</th><th>remarks</th></tr>
<tr><td>ℂP³ — labelled, oriented, coincidences allowed</td><td class="n">1</td><td class="n">0</td><td class="n">1</td><td class="n">0</td><td class="n">1</td><td class="n">0</td><td class="n">1</td><td>compact manifold, simply connected</td></tr>
<tr><td>M₀,₆ — labelled, oriented, points distinct</td><td class="n">1</td><td class="n">9</td><td class="n">26</td><td class="n">24</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td>K(π,1), π = pure braid group P₅ / centre; χ = −6</td></tr>
<tr><td>ℂP³/(S₅ × conj) — unlabelled, reflections, coincidences allowed</td><td class="n">1</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td class="n">1</td><td class="n">0</td><td class="n">0</td><td>compact orbifold; conjugation acts by (−1)ᵏ on H²ᵏ</td></tr>
<tr><td>M₀,₆/(S₅ × conj) — unlabelled, reflections, points distinct</td><td class="n">1</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td class="n">0</td><td>rationally a point; orbifold π₁ = B₅/⟨full twist⟩ ⋊ ℤ/2</td></tr>
</table></div>
<p>The invariants were computed on the Orlik–Solomon algebra of the braid arrangement (check: F(ℂ,5)/S₅ has the cohomology of a circle, Arnold's theorem). The orbifold points of the unlabelled space are the configurations with a symmetry; the mirror locus has codimension 3, so the one-parameter families below never fold.</p>
</section>

<section>
<h2>The five-point row</h2>
<div class="tw"><table>
<tr><th>d</th><th>rigid sets</th><th>families</th><th>expected dim 2n−4−C(n,2)+d</th><th>status</th></tr>
<tr><td class="n">2</td><td class="n">1</td><td>—</td><td class="n">−3</td><td>Einhorn–Schoenberg, re-verified</td></tr>
<tr><td class="n">3</td><td class="n">34</td><td>—</td><td class="n">−2</td><td>Shinohara 2004, re-verified</td></tr>
<tr><td class="n">4</td><td class="n">468</td><td class="n">25 curves</td><td class="n">−1</td><td>ours, certified</td></tr>
<tr><td class="n">5</td><td class="n">0</td><td class="n">402 curves, 8 surfaces</td><td class="n">1</td><td>ours, certified</td></tr>
<tr><td class="n">6</td><td class="n">0</td><td class="n">244 surfaces, 1 three-parameter</td><td class="n">2</td><td>ours, certified</td></tr>
<tr><td class="n">7</td><td class="n">0</td><td class="n">76 three-parameter</td><td class="n">3</td><td>ours, certified</td></tr>
</table></div>
<p>Each family is one realizable colouring type of K₅ (which pairs are at equal distance), and its dimension is that of its real solution set modulo similarity. The stratification question is how these pieces are glued: where a family's closure meets the strata with fewer distances.</p>
</section>

<section>
<h2>Four distances: a graph on 35 vertices</h2>
<div class="stats">
<div class="stat"><b class="num">35</b><span>vertices: the 34 three-distance sets and the pentagon</span></div>
<div class="stat"><b class="num">108</b><span>arcs between vertices (25 families)</span></div>
<div class="stat"><b class="num">77 + 12</b><span>half-open arcs: two points merge · one distance → 0</span></div>
<div class="stat"><b class="num">74</b><span>first Betti number; connected</span></div>
</div>
<figure>
{F['fig1']}
<figcaption>Every one-parameter family of five-point four-distance sets is a curve whose ends are three-distance sets (two of its four distance values merge), the pentagon, or a collision. Vertices are drawn as the configurations themselves, numbered as in the (5,3) catalog; teal arcs are the families, parallel arcs are distinct families between the same two sets; dashed rust stubs end with two points merging, violet stubs escape to infinity (the base pair merges after rescaling). The 468 rigid four-distance sets are isolated and not drawn. Layout is force-directed; positions carry no meaning.</figcaption>
</figure>
<ul class="legend">
{leg1}
</ul>
<p>Every three-distance set is reached, with between 4 and 22 arcs each. So the space of five-point sets with at most four distances, modulo similarity, is this graph plus 468 isolated points: H₀ = ℤ⁴⁶⁹, H₁ = ℤ⁷⁴, and the 89 half-open edges run off to the collision strata of ℂP³.</p>
</section>

<section>
<h2>Five distances: 402 curves on 527 vertices</h2>
<div class="stats">
<div class="stat"><b class="num">527</b><span>vertices: 468 + 34 + 1 rigid sets, 24 junction points</span></div>
<div class="stat"><b class="num">5269</b><span>arcs between vertices (402 families; ±17)</span></div>
<div class="stat"><b class="num">≈1900 + 105</b><span>half-open arcs: collisions · infinity</span></div>
<div class="stat"><b class="num">4743</b><span>first Betti number; connected</span></div>
</div>
<p>The same construction one level up. Two new phenomena: 24 curves end at an interior point of a four-distance <em>family</em> rather than at a rigid set (junctions, where a five-distance curve is attached to the middle of an edge of the previous graph), and many curves end at three-distance sets directly, i.e. two coincidences happen at once at a symmetric configuration. Every rigid four-distance set is reached.</p>
<div class="tw"><table>
<tr><th>arc between</th><th>count</th><th>arc between</th><th>count</th></tr>
<tr><td>two rigid four-distance sets</td><td class="n">3376</td><td>four-distance set and a collision</td><td class="n">1528</td></tr>
<tr><td>four-distance and three-distance set</td><td class="n">1653</td><td>three-distance set and a collision</td><td class="n">332</td></tr>
<tr><td>two three-distance sets</td><td class="n">137</td><td>a set and infinity</td><td class="n">105</td></tr>
<tr><td>four-distance set and a junction</td><td class="n">64</td><td>pentagon and a four-distance set</td><td class="n">35</td></tr>
<tr><td>three-distance set / collision and a junction</td><td class="n">4 + 10</td><td>collision to collision, lost traces</td><td class="n">18, 5</td></tr>
</table></div>
<figure>
{F['fig2']}
<figcaption>How many arcs meet at each vertex. Rigid four-distance sets typically carry 11–30 arcs (median 22); the three-distance sets are the hubs, up to 216 arcs at the most symmetric ones (hexagon-with-centre pieces, the 12-gon constructions); junction points carry 1–12.</figcaption>
</figure>
<ul class="legend">
{leg2}
</ul>
<p>The eight two-parameter families of (5,5) are bounded by these curves and by four-distance families, with three-distance sets at the corners; their boundary strata are listed exactly on their cards. The real topology of the surfaces themselves (genus, number of boundary components) has not been computed, so H₁ of the full space with at most five distances is still open; the 1-skeleton is what is above.</p>
</section>


<section>
<h2>The simpler cases: three and four points</h2>
<p>Three points. The shape space of triangles (unlabelled, reflections allowed) is a closed curvilinear triangle: the quotient of Kendall's sphere by the 12 symmetries. Its strata are exactly the table's (3,1), (3,2), (3,3): the equilateral triangle is a corner, the isosceles triangles form two sides meeting there (apex angle running from 60° up to 180° on one side and down to 0° on the other), and the scalene triangles are the interior together with the third side, the collinear triples, which is a mirror of the orbifold rather than a boundary of the stratum. One corner is a collision and not part of the space; the other, the collinear triple with the apex at the midpoint, is where the isosceles family folds: the family meets its own reflection and turns back.</p>
<figure>
{F34['fig3']}
<figcaption>Triangles up to similarity. Teal: the two isosceles arcs (one family, folded at the collinear corner); violet: collinear scalene triples, the mirror edge; the interior is the two-parameter scalene family. For three points the mirror locus has codimension 1, so the reflection quotient has a visible edge; for n points it has codimension n − 2.</figcaption>
</figure>
<p>Four points. The ambient space is four-dimensional (ℂP²/conj is a 4-sphere; its quotient by S₄ is the orbifold the table lives in). The six two-distance quadruples are the rigid points, the nine three-distance families are curves, the seven four-distance families surfaces, and the two five-distance families and the generic quadruples fill the rest. The graph of the curves is small enough to read off completely:</p>
<div class="stats">
<div class="stat"><b class="num">6 + 1</b><span>vertices: the two-distance quadruples, and one fold point</span></div>
<div class="stat"><b class="num">18</b><span>arcs between vertices (9 families)</span></div>
<div class="stat"><b class="num">15 + 2</b><span>half-open arcs: collisions · infinity</span></div>
<div class="stat"><b class="num">12</b><span>first Betti number; connected</span></div>
</div>
<figure>
{F34['fig4']}
<figcaption>The nine one-parameter families of four-point three-distance sets, coloured by family, between the six two-distance quadruples. Dashed stubs end in a collision. The fold vertex is the configuration of four equally spaced collinear points: the trapezoid-with-three-equal-sides family and the parallelogram-with-a-diagonal-equal-to-a-side family (the two branches of one colouring type) both run into it and turn back, so geometrically two arcs meet there.</figcaption>
</figure>
<ul class="legend">
{leg4}
</ul>
<ul class="legend">
{leg4v}
</ul>
<p>The surfaces are bounded by these curves (exact decomposition of each surface's closure): the circumcentre family (A the centre of a circle through B, C, D) by families 1 and 3; "ABC equilateral, D free" by 1 and 4; "AB = AC, AD = BC" by 1, 5, 6, 7; the equilateral path by 1, 2, 5, 6; the two-isosceles-triangles family by 1, 3, 4, 6, 8; kites by 2, 3, 4, 7; parallelograms by 2 (rhombi), 5, 9 (rectangles). Every curve bounds at least one surface, and "D on the unit circle about a vertex of an equilateral triangle" (family 1) bounds five of the seven.</p>
<p>So the pattern that the five-point computation found is already visible here: rigid sets are the vertices, each one-parameter family is a few arcs between them with ends running off to collisions, the graph is connected, and the only new ingredient at n = 5 is the junction (a curve ending in the interior of a lower family), which cannot happen for n ≤ 4 because the four-distance-and-lower families there are the rigid points and the curves of the same picture.</p>
</section>

<section>
<h2>The degeneration poset</h2>
{fig3}
<p class="figcap">Which strata lie in the closure of which: for every family and every coincidence equation (two distance values equal, or one distance zero) the components of the closure were computed and their colouring type identified at a certified real point. Each arrow counts pairs of types. Every lower family is in the boundary of some higher family at every level; the rigid sets of (5,4) and (5,3) appear as corners of the surfaces and three-folds, reached only partially by one-equation boundaries. Exact for d ≤ 5, partly sampled for d = 6, 7.</p>
<p>Reading down: a generic seven-distance three-fold has on its boundary two-parameter six-distance families (1152 incidences over the 76 types), occasionally five-distance curves, four-distance families and even rigid sets at special points. From the two-parameter six-distance families down: 2887 incidences with five-distance curves, 314 with four-distance families, 202 with rigid four-distance sets, 111 with three-distance sets. From the curves of (5,5): 2239 rigid four-distance ends, 94 junction-type ends on four-distance families, 596 three-distance ends, 12 at the pentagon. From the 25 four-distance curves: 83 three-distance ends and 7 at the pentagon.</p>
</section>

<section>
<h2>Method, in one paragraph each</h2>
<p><strong>Families.</strong> For each colouring type the squared-distance equations with the scale fixed, saturated by the inequations (all values distinct and nonzero) in Singular, give the closure of the family; the certification of the counts (msolve real root isolation, primary decomposition, hidden-point and node tests) is described in the applet's About text.</p>
<p><strong>Boundary points and arcs.</strong> Adding one coincidence equation to a curve's ideal gives a zero-dimensional system; msolve isolates its real roots to 200 bits and each root is matched up to similarity against the certified rigid lists. The real curve is then traced numerically from every root (tangent step and Newton projection; arrival detected by proximity to the exact roots; escape when coordinates exceed 300) and from random slice points to catch closed loops (none). Traversals are grouped into geometric arcs by a permutation- and similarity-invariant fingerprint of the shapes along the arc, which collapses the two directions and all symmetric copies of a labelled arc into one.</p>
<p><strong>Poset.</strong> Adding one coincidence equation to a higher-dimensional family's ideal gives its boundary; it was decomposed into irreducible components over ℚ (minAssGTZ, with a timeout) and each component's type read off a certified real point. Where the decomposition did not finish (most of d = 6, all of d = 7) random slices through the boundary supplied certified real points instead.</p>
<div class="caveat">Caveats. Vertex–vertex arc counts are exact up to the few arcs whose two traversals did not pair (17 of 5269 for five distances); collision-end counts are approximate (±40). Type 327 of (5,5) has 5 lost traces. Sampled poset data lists only strata that were hit, and small-integer slices hit symmetric special points with positive probability, so a sampled list can include a corner alongside the genuine boundary faces. The three (5,5) types with both a curve and a surface were traced on their curve components only.</div>
</section>

<section>
<h2>Where things are</h2>
<p><code>tools/fixed-distances/strat.py</code>, <code>strat_analyze.py</code>, <code>degen.py</code>, <code>degen_sample.py</code>, <code>degen_analyze2.py</code>, <code>sat5b.py</code>; data in <code>tools/fixed-distances/data/</code> (<code>strat_graph_4.json</code>, <code>strat_graph_5.json</code>, <code>degen_poset_*.json</code>, <code>strat*_analysis.txt</code>); the README there documents the pipeline and its limits. The cards in the applet's (5,4) and (5,5) catalogs list each family's arcs, and the cards of (5,3), (5,4), (5,6), (5,7) their incidences.</p>
</section>
</main>
'''
open('summary/five-point-strata.html','w').write(html); print(len(html))
