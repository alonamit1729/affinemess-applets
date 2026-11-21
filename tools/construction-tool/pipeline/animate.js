'use strict';
/* node animate.js <out dir> <cell> [<cell> ...]
   For every family (a set that flexes) with a construction in <out dir>/<cell>.json, a loop through the family
   that the catalog card and the tool's view mode both play, written to <out dir>/anim_<cell>.json {num: anim}.
   The unknowns are the recipe's free points (x, y) and sliders (u), in node order.
   - No equalities: every unknown swings on its own, x_j(t) = x_j + a_j sin(w_j t + p_j).    anim = {a, w, p}
   - With equalities: `dim` unknowns swing like that, the others follow by Newton's method; the loop is stored as
     K keyframes of all unknowns (linear in between is close enough at this density).            anim = {k: [[..], ..]}
   Along the whole loop the set must stay a set of n points with exactly d distances, none of them shrinking below
   a third of the shortest one; the amplitude is halved until it does (no animation if even a tiny swing fails). */
const fs = require('fs'), path = require('path');
const core = require('./core.js');
const { evaluate, distClasses, decodeRecipe } = core;

const K = 40, SUB = 6;
const W = [1, 2, 1, 3, 2, 1, 3, 2];
const P = [0, 1.3, 2.1, 0.7, 3.4, 4.2, 5.0, 2.8];

function unknownsOf(nodes) { return nodes.filter(n => n.t === 'free' || n.t === 'glide'); }
function x0Of(U) { return U.flatMap(n => n.t === 'free' ? [n.x, n.y] : [n.u]); }
function ovOf(U, x) {
  const ov = new Map(); let k = 0;
  for (const n of U) { if (n.t === 'free') { ov.set(n.id, { x: x[k], y: x[k + 1] }); k += 2; } else { ov.set(n.id, { u: x[k] }); k += 1; } }
  return ov;
}
/* scale of each unknown: a slider on a circle moves by angle, so divide by the radius */
function scales(nodes, U) {
  const V = evaluate(nodes), s = [];
  for (const n of U) {
    if (n.t === 'free') s.push(1, 1);
    else { const c = V.get(n.on); s.push(c && c.k === 'c' ? 1 / Math.max(0.2, c.r) : 1); }
  }
  return s;
}
function pointsAt(nodes, rec, U, x) {
  const V = evaluate(nodes, ovOf(U, x));
  const P0 = rec.m.map(id => V.get(id));
  return P0.every(p => p && p.x !== undefined && Number.isFinite(p.x) && Number.isFinite(p.y)) ? P0 : null;
}
function shape(P0) {
  let mn = Infinity, R = 0; const c = { x: 0, y: 0 };
  P0.forEach(p => { c.x += p.x / P0.length; c.y += p.y / P0.length; });
  for (let i = 0; i < P0.length; i++) { R = Math.max(R, Math.hypot(P0[i].x - c.x, P0[i].y - c.y)); for (let j = i + 1; j < P0.length; j++) mn = Math.min(mn, Math.hypot(P0[i].x - P0[j].x, P0[i].y - P0[j].y)); }
  return { mn, R };
}
function eqs(rec, V) {
  const out = [];
  for (const g of rec.q || []) {
    const d2 = sg => { const a = V.get(sg[0]), b = V.get(sg[1]); return (a.x - b.x) ** 2 + (a.y - b.y) ** 2; };
    const b0 = d2(g[0]); for (let k = 1; k < g.length; k++) out.push(d2(g[k]) - b0);
  }
  return out;
}
function F(nodes, rec, U, x) { try { const r = eqs(rec, evaluate(nodes, ovOf(U, x))); return r.every(Number.isFinite) ? r : null; } catch (_) { return null; } }
function jac(nodes, rec, U, x, cols) {
  const h = 1e-7, J = [];
  for (const c of cols) {
    const xp = x.slice(), xm = x.slice(); xp[c] += h; xm[c] -= h;
    const fp = F(nodes, rec, U, xp), fm = F(nodes, rec, U, xm); if (!fp || !fm) return null;
    fp.forEach((v, i) => { (J[i] = J[i] || []).push((v - fm[i]) / (2 * h)); });
  }
  return J;
}
function solve(A, b) {   // square, partial pivoting
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c; for (let i = c + 1; i < n; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i;
    if (Math.abs(M[p][c]) < 1e-14) return null;
    [M[c], M[p]] = [M[p], M[c]];
    for (let i = c + 1; i < n; i++) { const f = M[i][c] / M[c][c]; for (let k = c; k <= n; k++) M[i][k] -= f * M[c][k]; }
  }
  const x = new Array(n); for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let k = i + 1; k < n; k++) s -= M[i][k] * x[k]; x[i] = s / M[i][i]; }
  return x;
}
/* Newton on the pivot unknowns (as many as equations), the driven ones held */
function newton(nodes, rec, U, x, piv) {
  x = x.slice();
  for (let it = 0; it < 30; it++) {
    const f = F(nodes, rec, U, x); if (!f) return null;
    if (Math.max(...f.map(Math.abs)) < 1e-12) return x;
    const J = jac(nodes, rec, U, x, piv); if (!J) return null;
    const dx = solve(J, f.map(v => -v)); if (!dx) return null;
    piv.forEach((c, k) => { x[c] += dx[k]; });
  }
  const f = F(nodes, rec, U, x); return f && Math.max(...f.map(Math.abs)) < 1e-10 ? x : null;
}
/* pivot columns: a maximal independent set of columns of the Jacobian */
function pivots(nodes, rec, U, x) {
  const m = x.length, all = [...Array(m).keys()], J = jac(nodes, rec, U, x, all); if (!J) return null;
  const piv = [];
  for (const c of all) { const cols = [...piv, c]; const sub = J.map(r => cols.map(k => r[k])); if (core.rank(sub, cols.length) === cols.length) piv.push(c); if (piv.length === J.length) break; }
  return piv.length === J.length ? piv : null;
}

function animate(res, job, d) {
  const rec = res.rec, nodes = decodeRecipe(rec), U = unknownsOf(nodes); if (!U.length) return null;
  const x0 = x0Of(U), m = x0.length, sc = scales(nodes, U);
  const P0 = pointsAt(nodes, rec, U, x0); if (!P0) return null;
  const s0 = shape(P0);
  const ok = x => {
    const P1 = pointsAt(nodes, rec, U, x); if (!P1) return false;
    const dc = distClasses(P1); if (dc.cls.length !== d || dc.cls[0].d < 1e-9) return false;
    const s = shape(P1); return s.mn > s0.mn / 3 && s.R < s0.R * 2.2;
  };
  const hasEq = (rec.q || []).length > 0;
  const piv = hasEq ? pivots(nodes, rec, U, x0) : [];
  if (piv === null) return null;
  const drv = [...Array(m).keys()].filter(c => !piv.includes(c));
  if (!drv.length) return null;
  for (let A = 0.32; A >= 0.004; A /= 2) {
    const a = x0.map((_, j) => drv.includes(j) ? +(A * sc[j]).toPrecision(4) : 0);
    const at = (t, prev) => {
      const x = x0.map((v, j) => v + a[j] * Math.sin(W[j % 8] * t + P[j % 8]) - a[j] * Math.sin(P[j % 8]));   // starts at x0
      if (!hasEq) return x;
      piv.forEach(c => { x[c] = prev[c]; });
      return newton(nodes, rec, U, x, piv);
    };
    let prev = x0, good = true; const frames = [];
    for (let i = 0; i <= K * SUB && good; i++) {
      const t = 2 * Math.PI * i / (K * SUB), x = at(t, prev);
      if (!x || !ok(x)) { good = false; break; }
      if (i % SUB === 0 && i < K * SUB) frames.push(x.map(v => +v.toFixed(7)));
      prev = x;
    }
    if (!good) continue;
    if (hasEq && Math.max(...prev.map((v, j) => Math.abs(v - x0[j]))) > 1e-6) continue;   // the loop must close
    if (!hasEq) return { a, w: a.map((_, j) => W[j % 8]), p: a.map((_, j) => P[j % 8]) };
    return { k: frames };
  }
  return null;
}

if (require.main === module) {
  const [outDir, ...cells] = process.argv.slice(2);
  for (const cell of cells) {
    const f = path.join(outDir, cell + '.json'); if (!fs.existsSync(f)) continue;
    const O = JSON.parse(fs.readFileSync(f)), [, d] = cell.split('_').map(Number);
    const A = {}; let n = 0, got = 0;
    for (const r of Object.values(O)) {
      if (!r.ok) continue;
      if (!(r.free > 0)) continue;                    // families only (free directions left)
      n++; const a = animate(r, null, d); if (a) { A[r.num] = a; got++; }
    }
    fs.writeFileSync(path.join(outDir, 'anim_' + cell + '.json'), JSON.stringify(A));
    console.log(cell, 'families', n, 'animated', got);
  }
}
module.exports = { animate };
