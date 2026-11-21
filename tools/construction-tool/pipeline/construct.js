'use strict';
/* Find a construction of a catalog set with the construction tool's tools, preferring, in order:
     1  a piece of the triangular or square lattice (no construction needed)
     2  straightedge and compass (a regular polygon counts when straightedge and compass can draw it: 3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30)
     3  straightedge, compass and the regular-polygon tool (7, 9, 11, 13, 14, 18, 36 ...)
     4  all tools, including free points, sliders and the equalizer
   Within a level: fewest equations, then fewest unknowns, then fewest steps.
   Families (sets that flex) always end with free points or sliders: as many free directions as the family's dimension.

   The search: a start (two points / a lattice / a regular polygon with its centre), then, repeatedly, a point of
   the set where two available lines or circles meet: a line through two available points, a perpendicular
   bisector, a circle about an available point through another, or a compass circle whose radius is a distance
   between available points. When that gets stuck, a slider on a line or circle through a missing point (or a free
   point, if there is none) and on again; at the end, the fewest equalities that pin the sliders (or leave exactly
   the family's freedom). Only points of the set and of the start are used as stepping stones. */
const core = require('./core.js');
const { evaluate, meet, curveVal, polyGeom, latVec, distClasses, canon, sameCanon } = core;

const TOL = 1e-8;
const GAUSS = new Set([3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34]);
const MS = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 24, 30, 36];

/* ---------- complex helpers ---------- */
const cx = (x, y) => ({ x, y });
const sub = (a, b) => cx(a.x - b.x, a.y - b.y);
const mul = (a, b) => cx(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
const div = (a, b) => { const d = b.x * b.x + b.y * b.y; return cx((a.x * b.x + a.y * b.y) / d, (a.y * b.x - a.x * b.y) / d); };
const abs = a => Math.hypot(a.x, a.y);
const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const toFrame = (z, A, B) => div(sub(z, A), sub(B, A));

/* ---------- builder: a growing recipe ---------- */
class Builder {
  constructor(T) { this.T = T; this.nodes = []; this.val = new Map(); this.tn = new Array(T.length).fill(-1); this.pts = []; this.curves = []; this.lat = null; }
  clone() {
    const b = new Builder(this.T);
    b.nodes = this.nodes.map(n => ({ ...n })); b.val = new Map(this.val); b.tn = this.tn.slice();
    b.pts = this.pts.slice(); b.curves = this.curves.slice(); b.lat = this.lat; return b;
  }
  P(id) { const v = this.val.get(id); return v && v.x !== undefined ? v : null; }
  add(n) {
    n.id = this.nodes.length; this.nodes.push(n);
    let v = null;
    switch (n.t) {
      case 'base': case 'free': v = cx(n.x, n.y); break;
      case 'glide': { const c = this.val.get(n.on); v = c.k === 'l' ? cx(c.px + n.u * c.dx, c.py + n.u * c.dy) : cx(c.cx + c.r * Math.cos(n.u), c.cy + c.r * Math.sin(n.u)); break; }
      case 'lat': { const a = this.P(n.a), b = this.P(n.b), [ux, uy, vx, vy] = latVec(n.kind, a, b); v = cx(a.x + n.i * ux + n.j * vx, a.y + n.i * uy + n.j * vy); break; }
      case 'poly': v = { poly: polyGeom(this.P(n.a), this.P(n.b), n.m) }; break;
      case 'pv': v = this.val.get(n.poly).poly.verts[n.k]; break;
      case 'pc': { const g = this.val.get(n.poly).poly; v = cx(g.cx, g.cy); break; }
      case 'x': { const m = meet(this.val.get(n.c1), this.val.get(n.c2)); v = m[Math.min(n.k, m.length - 1)]; break; }
      default: v = curveVal(n, id => this.P(id), this.val);
    }
    this.val.set(n.id, v);
    if (v && v.x !== undefined) {
      this.pts.push(n.id);
      this.T.forEach((t, k) => { if (this.tn[k] < 0 && near(v, t)) this.tn[k] = n.id; });
    } else if (v && v.k) this.curves.push(n.id);
    return n.id;
  }
  done() { return this.tn.every(i => i >= 0); }
  reached() { return this.tn.filter(i => i >= 0).length; }
}

/* ---------- curves through a point ---------- */
function curveKey(v) {
  if (v.k === 'c') return 'c' + v.cx.toFixed(7) + ',' + v.cy.toFixed(7) + ',' + v.r.toFixed(7);
  let dx = v.dx, dy = v.dy; if (dy < -1e-12 || (Math.abs(dy) <= 1e-12 && dx < 0)) { dx = -dx; dy = -dy; }
  const off = v.px * dy - v.py * dx;
  return 'l' + dx.toFixed(7) + ',' + dy.toFixed(7) + ',' + off.toFixed(7);
}
const onCurve = (v, p, tol = TOL) => v.k === 'l' ? Math.abs((p.x - v.px) * v.dy - (p.y - v.py) * v.dx) < tol * Math.max(1, abs(p))
  : Math.abs(dist(p, cx(v.cx, v.cy)) - v.r) < tol * Math.max(1, v.r);
/* sine of the angle at which two curves cross at p */
function crossing(a, b, p) {
  const dir = v => v.k === 'l' ? cx(v.dx, v.dy) : (() => { const r = cx(p.x - v.cx, p.y - v.cy), L = abs(r); return cx(-r.y / L, r.x / L); })();
  const u = dir(a), w = dir(b); return Math.abs(u.x * w.y - u.y * w.x);
}
const near = (q, p) => dist(q, p) < 1e-6 * Math.max(1, abs(p));

/* distances between available points: key -> [q, r] */
function distIndex(B) {
  const D = new Map(), P = B.pts;
  for (let i = 0; i < P.length; i++) for (let j = i + 1; j < P.length; j++) {
    const d = dist(B.P(P[i]), B.P(P[j])); if (d < 1e-9) continue;
    const k = d.toFixed(7); if (!D.has(k)) D.set(k, [P[i], P[j]]);
  }
  return D;
}
function curvesThrough(B, p, D) {
  const out = new Map();   // key -> {spec | id, v, cost}
  // a new curve's value is computed exactly as evaluate() will compute it
  const put = (spec, v, cost, id) => {
    if (spec) { v = curveVal(spec, q => B.P(q), B.val); if (!v || !onCurve(v, p, 1e-6)) return; }
    const k = curveKey(v), o = out.get(k); if (!o || cost < o.cost) out.set(k, { spec, v, cost, id });
  };
  for (const id of B.curves) { const v = B.val.get(id); if (onCurve(v, p)) put(null, v, 0, id); }
  const P = B.pts, Q = P.map(id => B.P(id));
  // parallels and perpendiculars through an available point, along a segment between two available points
  const dirs = new Map(), SC = 1e6, KMAX = Math.round(Math.PI * SC);
  const dkey = (dx, dy) => { let t = Math.atan2(dy, dx); if (t < 0) t += Math.PI; if (t >= Math.PI) t -= Math.PI; return Math.round(t * SC) % KMAX; };
  for (let a = 0; a < P.length; a++) for (let b = a + 1; b < P.length; b++) {
    const dx = Q[b].x - Q[a].x, dy = Q[b].y - Q[a].y; if (Math.hypot(dx, dy) < 1e-9) continue;
    const k = dkey(dx, dy); if (!dirs.has(k)) dirs.set(k, [a, b]);
  }
  const findDir = (dx, dy) => { const k = dkey(dx, dy); for (const kk of [k, k - 1, k + 1, (k + KMAX - 1) % KMAX, (k + 1) % KMAX]) if (dirs.has(kk)) return dirs.get(kk); return null; };
  for (let a = 0; a < P.length; a++) {
    const A = Q[a], dx = p.x - A.x, dy = p.y - A.y; if (Math.hypot(dx, dy) < 1e-9) continue;
    const par = findDir(dx, dy); if (par) put({ t: 'par', a: P[par[0]], b: P[par[1]], p: P[a] }, null, 1.05, null);
    const per = findDir(-dy, dx); if (per) put({ t: 'perp', a: P[per[0]], b: P[per[1]], p: P[a] }, null, 1.05, null);
  }
  for (let a = 0; a < P.length; a++) {
    const A = Q[a], r = dist(p, A); if (r < 1e-9) continue;
    // circle about A through another available point at the same distance
    let circ = false;
    for (let b = 0; b < P.length; b++) if (b !== a && Math.abs(dist(Q[b], A) - r) < TOL * Math.max(1, r)) {
      put({ t: 'circ', c: P[a], p: P[b] }, { k: 'c', cx: A.x, cy: A.y, r }, 1, null); circ = true; break;
    }
    if (!circ) { const qr = D.get(r.toFixed(7)); if (qr && Math.abs(dist(B.P(qr[0]), B.P(qr[1])) - r) < TOL * Math.max(1, r)) put({ t: 'comp', p: qr[0], q: qr[1], c: P[a] }, { k: 'c', cx: A.x, cy: A.y, r }, 1.1, null); }
    for (let b = a + 1; b < P.length; b++) {
      const Bp = Q[b], L = dist(A, Bp); if (L < 1e-9) continue;
      const ux = (Bp.x - A.x) / L, uy = (Bp.y - A.y) / L;
      if (Math.abs((p.x - A.x) * uy - (p.y - A.y) * ux) < TOL * Math.max(1, abs(p))) put({ t: 'line', a: P[a], b: P[b] }, { k: 'l', px: A.x, py: A.y, dx: ux, dy: uy }, 1, null);
      if (Math.abs(dist(p, A) - dist(p, Bp)) < TOL * Math.max(1, r)) put({ t: 'bis', a: P[a], b: P[b] }, { k: 'l', px: (A.x + Bp.x) / 2, py: (A.y + Bp.y) / 2, dx: -uy, dy: ux }, 1.05, null);
    }
  }
  return [...out.values()].sort((u, v) => u.cost - v.cost);
}
/* the cheapest pair of curves through p that meet transversally there */
function bestPair(C, p) {
  let best = null;
  for (let i = 0; i < C.length; i++) for (let j = i + 1; j < C.length; j++) {
    const sn = crossing(C[i].v, C[j].v, p); if (sn < 0.02) continue;            // nearly tangent: avoid (numerically fragile)
    const cost = C[i].cost + C[j].cost + (sn < 0.2 ? 0.5 : 0); if (best && cost >= best.cost) continue;
    const m = meet(C[i].v, C[j].v);
    if (!m.length || (m.length === 1 && !(C[i].v.k === 'l' && C[j].v.k === 'l'))) continue;
    const k = m.findIndex(q => near(q, p)); if (k < 0) continue;
    best = { a: C[i], b: C[j], k, cost };
  }
  return best;
}
function curveId(B, c) { return c.id !== null && c.id !== undefined ? c.id : B.add({ ...c.spec }); }

/* add points of the set while two available curves meet at one */
function closure(B) {
  for (;;) {
    const D = distIndex(B);
    let best = null;
    B.T.forEach((p, k) => {
      if (B.tn[k] >= 0) return;
      const C = curvesThrough(B, p, D); if (C.length < 2) return;
      const bp = bestPair(C, p); if (bp && (!best || bp.cost < best.bp.cost)) best = { k, bp };
    });
    if (!best) return B;
    const { a, b, k } = best.bp;
    const ia = curveId(B, a), ib = curveId(B, b);
    // the branch index refers to meet(c1, c2) in this order
    const m = meet(B.val.get(ia), B.val.get(ib)), kk = m.findIndex(q => near(q, B.T[best.k]));
    B.add({ t: 'x', c1: ia, c2: ib, k: kk });
    if (B.tn[best.k] < 0) throw new Error('closure: point not reached');
  }
}

/* ---------- starts ---------- */
function pairStart(T, i, j) {
  const F = T.map(z => toFrame(z, T[i], T[j]));
  const B = new Builder(F); B.add({ t: 'base', x: 0, y: 0 }); B.add({ t: 'base', x: 1, y: 0 });
  return B;
}
const OMEGA = cx(0.5, Math.sqrt(3) / 2);
function latticeCoords(kind, w) {            // w = s + t·v with v = ω (tri) or i (sq)
  if (kind === 'sq') return [w.x, w.y];
  const t = w.y / OMEGA.y; return [w.x - t * OMEGA.x, t];
}
const isInt = v => Math.abs(v - Math.round(v)) < 1e-7;
/* rational approximation with a small denominator */
function rat(v, maxDen = 60) {
  for (let q = 1; q <= maxDen; q++) { const p = Math.round(v * q); if (Math.abs(v * q - p) < 1e-7 * q) return [p, q]; }
  return null;
}
/* gcd in Z[i] / Z[ω] (Euclidean), elements as integer pairs [s, t] meaning s + t·v */
function egcd(kind, a, b) {
  const norm = ([s, t]) => kind === 'sq' ? s * s + t * t : s * s + s * t + t * t;      // ω = e^{iπ/3}: N(s + tω) = s² + st + t²
  const toC = ([s, t]) => kind === 'sq' ? cx(s, t) : cx(s + t * OMEGA.x, t * OMEGA.y);
  const fromC = z => { const [s, t] = latticeCoords(kind, z); return [Math.round(s), Math.round(t)]; };
  while (norm(b) > 0) {
    const q = fromC(div(toC(a), toC(b))), qb = fromC(mul(toC(q), toC(b)));
    const r = [a[0] - qb[0], a[1] - qb[1]];
    if (norm(r) >= norm(b)) { // nearest-integer rounding can miss in rare cases: try the neighbours
      let ok = false;
      for (const [ds, dt] of [[1, 0], [-1, 0], [0, 1], [0, -1], [1, -1], [-1, 1]]) {
        const q2 = [q[0] + ds, q[1] + dt], qb2 = fromC(mul(toC(q2), toC(b))), r2 = [a[0] - qb2[0], a[1] - qb2[1]];
        if (norm(r2) < norm(b)) { a = b; b = r2; ok = true; break; }
      }
      if (!ok) return null;
      continue;
    }
    a = b; b = r;
  }
  return a;
}
/* the whole set on a lattice: returns a level-1 builder */
function latticeFull(T) {
  const n = T.length;
  for (const kind of ['tri', 'sq']) {
    // frame: T[0] -> 0, T[1] -> 1; every point must have rational lattice coordinates
    const W = T.map(z => toFrame(z, T[0], T[1]));
    const R = W.map(w => latticeCoords(kind, w).map(v => rat(v)));
    if (R.some(r => !r[0] || !r[1])) continue;
    let L = 1; const g2 = (a, b) => b ? g2(b, a % b) : Math.abs(a);
    for (const r of R) for (const [, q] of r) L = L / g2(L, q) * q;
    const I = R.map(r => r.map(([p, q]) => p * (L / q)));
    let g = [0, 0]; for (const v of I) { if (v[0] || v[1]) g = g[0] || g[1] ? egcd(kind, g, v) : v; if (!g) break; }
    if (!g) continue;
    // base vector in the W frame: g / L
    const toC = ([s, t]) => kind === 'sq' ? cx(s, t) : cx(s + t * OMEGA.x, t * OMEGA.y);
    const v = div(toC(g), cx(L, 0));
    const F = W.map(w => div(w, v));                   // lattice frame: base points 0 and 1
    if (!F.every(f => latticeCoords(kind, f).every(isInt))) continue;
    // move the origin to a point of the set nearest the centre, so the base pair sits inside the figure
    const B = new Builder(F); B.add({ t: 'base', x: 0, y: 0 }); B.add({ t: 'base', x: 1, y: 0 }); B.lat = kind;
    F.forEach((f, k) => {
      if (B.tn[k] >= 0) return;
      const [i, j] = latticeCoords(kind, f).map(Math.round);
      B.add({ t: 'lat', kind, a: 0, b: 1, i, j });
    });
    if (B.done()) return B;
  }
  return null;
}
/* partial lattice start: base pair (i, j), lattice points of the set */
function latticePartial(T, i, j, kind) {
  const B = pairStart(T, i, j); let hits = 0;
  B.T.forEach((f, k) => {
    if (B.tn[k] >= 0) return;
    const c = latticeCoords(kind, f); if (!c.every(isInt)) return;
    B.add({ t: 'lat', kind, a: 0, b: 1, i: Math.round(c[0]), j: Math.round(c[1]) }); hits++;
  });
  if (hits < 1) return null;
  B.lat = kind; return B;
}
/* regular polygon starts: polygons through two points of the set (k steps apart), deduplicated */
function polygonStarts(T, ms) {
  const n = T.length, seen = new Set(), out = [];
  for (const m of ms) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (i === j) continue;
    for (let k = 1; k <= Math.floor(m / 2); k++) {
      const th = 2 * Math.PI * k / m, w = cx(Math.cos(th), Math.sin(th));
      // centre c with T[j] - c = w (T[i] - c)
      const c = div(sub(T[j], mul(w, T[i])), cx(1 - w.x, -w.y));
      const R = dist(T[i], c); if (R < 1e-9) continue;
      const ang = Math.atan2(T[i].y - c.y, T[i].x - c.x), rot = ((ang % (2 * Math.PI / m)) + 2 * Math.PI / m) % (2 * Math.PI / m);
      const key = m + ':' + c.x.toFixed(6) + ',' + c.y.toFixed(6) + ',' + R.toFixed(6) + ',' + Math.min(rot, 2 * Math.PI / m - rot).toFixed(6);
      if (seen.has(key)) continue; seen.add(key);
      // vertices and the hits among the set's points
      const V = []; for (let t = 0; t < m; t++) { const a = ang + 2 * Math.PI * t / m; V.push(cx(c.x + R * Math.cos(a), c.y + R * Math.sin(a))); }
      const hits = T.filter(p => dist(p, c) < 1e-7 * Math.max(1, R) || V.some(q => dist(p, q) < 1e-7 * Math.max(1, R))).length;
      // base: a vertex that is a point of the set, and the next vertex counterclockwise
      let t0 = V.findIndex(q => T.some(p => dist(p, q) < 1e-7 * Math.max(1, R))); if (t0 < 0) t0 = 0;
      out.push({ m, A: V[t0], B: V[(t0 + 1) % m], hits });
    }
  }
  return out;
}
function polyBuild(T, s) {
  const F = T.map(z => toFrame(z, s.A, s.B));
  const B = new Builder(F); B.add({ t: 'base', x: 0, y: 0 }); B.add({ t: 'base', x: 1, y: 0 });
  const pid = B.add({ t: 'poly', a: 0, b: 1, m: s.m });
  for (let k = 2; k < s.m; k++) B.add({ t: 'pv', poly: pid, k });
  B.add({ t: 'pc', poly: pid });
  B.m = s.m; return B;
}

/* ---------- unknowns and equalities ---------- */
function unknowns(B) { return B.nodes.filter(n => n.t === 'free' || n.t === 'glide'); }
function evalWith(nodes, U, x) {
  const ov = new Map(); let k = 0;
  for (const n of U) { if (n.t === 'free') { ov.set(n.id, { x: x[k], y: x[k + 1] }); k += 2; } else { ov.set(n.id, { u: x[k] }); k += 1; } }
  return evaluate(nodes, ov);
}
const xOf = U => U.flatMap(n => n.t === 'free' ? [n.x, n.y] : [n.u]);
/* rows: for each class, d²(e) − d²(ref) for the class's other pairs; their gradients in the unknowns */
function eqRows(B, pairs, cls) {
  const ref = new Map(), rows = [];
  pairs.forEach((p, e) => { const c = cls[e]; if (!ref.has(c)) ref.set(c, e); else rows.push([ref.get(c), e]); });
  return rows;
}
function jacobian(B, U, x, rows, pairs, h = 1e-6) {
  const cols = x.length;
  const d2 = (val, e) => { const a = val.get(B.tn[pairs[e][0]]), b = val.get(B.tn[pairs[e][1]]); return (a.x - b.x) ** 2 + (a.y - b.y) ** 2; };
  const J = rows.map(() => new Float64Array(cols));
  for (let c = 0; c < cols; c++) {
    const xp = x.slice(), xm = x.slice(); xp[c] += h; xm[c] -= h;
    const vp = evalWith(B.nodes, U, xp), vm = evalWith(B.nodes, U, xm);
    rows.forEach(([r0, e], i) => {
      try { J[i][c] = ((d2(vp, e) - d2(vp, r0)) - (d2(vm, e) - d2(vm, r0))) / (2 * h); } catch (_) { J[i][c] = NaN; }
    });
  }
  return J;
}
function rankOf(rows, m) { return rows.length ? core.rank(rows, m) : 0; }

/* ---------- the stuck step: a slider (or a free point) at a missing point ---------- */
function extend(B, maxUnknowns = 12) {
  let guard = 0;
  while (!B.done() && guard++ < 40) {
    const D = distIndex(B), opts = [];
    B.T.forEach((p, k) => {
      if (B.tn[k] >= 0) return;
      const C = curvesThrough(B, p, D);
      for (const c of C.slice(0, 3)) {
        const b = B.clone(), on = curveId(b, c), v = b.val.get(on);
        const u = v.k === 'l' ? (p.x - v.px) * v.dx + (p.y - v.py) * v.dy : Math.atan2(p.y - v.cy, p.x - v.cx);
        b.add({ t: 'glide', on, u }); closure(b); opts.push({ b, score: b.reached() * 10 - 1 - c.cost * 0.1 });
      }
      if (!C.length) { const b = B.clone(); b.add({ t: 'free', x: p.x, y: p.y }); closure(b); opts.push({ b, score: b.reached() * 10 - 2 }); }
    });
    if (!opts.length) return null;
    opts.sort((u, v) => v.score - u.score || u.b.nodes.length - v.b.nodes.length);
    B = opts[0].b;
    const nu = unknowns(B).reduce((s, n) => s + (n.t === 'free' ? 2 : 1), 0);
    if (nu > maxUnknowns) return null;
  }
  return B.done() ? B : null;
}
/* equalities: greedily the rows that raise the rank, up to (unknowns − dim) */
function chooseEqualities(B, job) {
  const n = job.n, pairs = []; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j]);
  const U = unknowns(B), x = xOf(U), m = x.length;
  if (!m) return { cons: [], eqs: 0, unk: 0, free: 0 };
  const rows = eqRows(B, pairs, job.cls);
  const J = jacobian(B, U, x, rows, pairs);
  // an equation whose gradient is numerical noise holds whatever the unknowns do: it pins nothing
  const okRows = rows.map((r, i) => i).filter(i => J[i].every(Number.isFinite) && Math.max(...J[i].map(Math.abs)) > 1e-5);
  const full = rankOf(okRows.map(i => J[i]), m);
  const want = m - job.dim;
  const chosen = []; let rk = 0;
  // prefer rows whose segments are short (the visible ones) and that reuse a class already in use
  const order = okRows.slice().sort((a, b) => {
    const len = i => { const [r0, e] = rows[i]; const P = q => B.T[q]; return dist(P(pairs[e][0]), P(pairs[e][1])); };
    return len(a) - len(b);
  });
  for (const i of order) {
    if (rk >= want) break;
    const r2 = rankOf([...chosen.map(c => J[c]), J[i]], m);
    if (r2 > rk) { chosen.push(i); rk = r2; }
  }
  // check: the chosen equations, differentiated again with the page's step (1e-7), have full rank
  if (chosen.length) {
    const J2 = jacobian(B, U, x, chosen.map(i => rows[i]), pairs, 1e-7);
    if (rankOf(J2, m) < chosen.length || J2.some(r => Math.max(...r.map(Math.abs)) <= 1e-5)) return { cons: [], eqs: chosen.length, unk: m, free: -1 };
  }
  // groups by class: [ref, e1, e2, ...] as node-id segments
  const groups = new Map();
  for (const i of chosen) {
    const [r0, e] = rows[i], c = job.cls[e];
    if (!groups.has(c)) groups.set(c, [r0]);
    groups.get(c).push(e);
  }
  const seg = e => [B.tn[pairs[e][0]], B.tn[pairs[e][1]]];
  const cons = [...groups.values()].map(g => g.map(seg));
  return { cons, eqs: chosen.length, unk: m, free: m - rk, fullRank: full, want };
}

/* ---------- the whole search ---------- */
function levelOf(B, eqs) {
  if (eqs) return 4;
  if (B.isLattice) return 1;
  const m = B.nodes.filter(n => n.t === 'poly').map(n => n.m);
  if (m.some(x => !GAUSS.has(x))) return 3;
  return 2;
}
function score(c) { return [c.level, c.eqs, c.unk, c.B.nodes.length]; }
const better = (a, b) => { const A = score(a), Bs = score(b); for (let i = 0; i < A.length; i++) if (A[i] !== Bs[i]) return A[i] < Bs[i]; return false; };

function construct(job) {
  const T = []; for (let k = 0; k < job.flat.length; k += 2) T.push(cx(job.flat[k], job.flat[k + 1]));
  const n = T.length; job.n = n;
  const rigid = job.kind === 'rigid';
  const minD = Math.min(...T.flatMap((a, i) => T.slice(i + 1).map(b => dist(a, b))));
  // 1. lattice
  if (rigid) { const L = latticeFull(T); if (L) { L.isLattice = true; return finish(L, job, 1); } }
  // 2./3. straightedge and compass, then polygons
  const starts = [];
  const pairs = []; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j, dist(T[i], T[j])]);
  pairs.sort((a, b) => a[2] - b[2]);
  for (const [i, j] of pairs) {
    starts.push(() => pairStart(T, i, j));
    for (const kind of ['tri', 'sq']) starts.push(() => latticePartial(T, i, j, kind));
  }
  // the catalog's own frame: its coordinates are written relative to a regular m-gon (circumradius 1, centre 0, a vertex at angle 0)
  const fm = /regular (\d+)-gon frame/.exec(job.frame || '');
  if (fm) { const m = +fm[1]; starts.unshift(() => polyBuild(T, { m, A: cx(1, 0), B: cx(Math.cos(2 * Math.PI / m), Math.sin(2 * Math.PI / m)) })); }
  let best = null; const tried = [];
  const consider = B => {
    if (!B) return;
    closure(B); tried.push(B);
    if (B.done()) {
      const c = { B, level: levelOf(B, 0), eqs: 0, unk: 0 };
      if ((!best || better(c, best)) && !verify(encode(B, []), job)) best = c;   // only constructions that replay exactly
    }
  };
  if (rigid) {
    for (const s of starts) consider(s());
    // polygons: those through at least three points of the set, most points first; the first that works is kept
    // then those through two points, whose other vertices and centre serve as stepping stones
    for (const [gauss, minHits] of [[true, 3], [false, 3], ...(process.env.DEEP ? [[true, 2], [false, 2]] : [])]) {
      if (best) break;
      const P = polygonStarts(T, MS.filter(m => GAUSS.has(m) === gauss)).filter(s => minHits === 3 ? s.hits >= 3 : s.hits === 2).sort((a, b) => b.hits - a.hits || a.m - b.m);
      for (const s of P) { consider(polyBuild(T, s)); if (best) break; }
    }
    if (best) return finish(best.B, job, best.level);
  } else {
    for (const s of starts) { const B = s(); if (B) { closure(B); tried.push(B); } }
    const P = polygonStarts(T, MS).filter(s => s.hits >= 3).sort((a, b) => b.hits - a.hits).slice(0, 30);
    for (const s of P) { const B = polyBuild(T, s); closure(B); tried.push(B); }
  }
  // 4. sliders, free points and equalities, from the starts that got furthest
  tried.sort((a, b) => b.reached() - a.reached() || a.nodes.length - b.nodes.length);
  const seen = new Set(); let cands = [];
  for (const B of tried) {
    const key = B.reached() + ':' + B.tn.map(i => i >= 0 ? 1 : 0).join('');
    if (seen.has(key)) continue; seen.add(key); cands.push(B);
    if (cands.length >= 6) break;
  }
  for (const B0 of cands) {
    const B = extend(B0.clone()); if (!B) continue;
    const E = chooseEqualities(B, job);
    if (E.free !== job.dim) continue;
    const c = { B, level: E.eqs ? 4 : levelOf(B, 0), eqs: E.eqs, unk: E.unk, E };
    if ((!best || better(c, best)) && !verify(encode(B, E.cons), job)) best = c;
  }
  if (best) return finish(best.B, job, best.level, best.E);
  return { ok: false, why: 'no construction found' };
}

/* ---------- output: compact recipe, checked by replaying it ---------- */
const R8 = v => +v.toPrecision(14);
function encode(B, cons) {
  const code = n => {
    switch (n.t) {
      case 'base': return ['b', R8(n.x), R8(n.y)];
      case 'free': return ['f', R8(n.x), R8(n.y)];
      case 'glide': return ['g', n.on, R8(n.u)];
      case 'lat': return [n.kind === 'tri' ? 't' : 's', n.a, n.b, n.i, n.j];
      case 'line': return ['l', n.a, n.b];
      case 'bis': return ['bis', n.a, n.b];
      case 'circ': return ['c', n.c, n.p];
      case 'comp': return ['k', n.p, n.q, n.c];
      case 'poly': return ['poly', n.a, n.b, n.m];
      case 'pv': return ['pv', n.poly, n.k];
      case 'pc': return ['pc', n.poly];
      case 'x': return ['x', n.c1, n.c2, n.k];
      case 'par': return ['par', n.a, n.b, n.p];
      case 'perp': return ['perp', n.a, n.b, n.p];
    }
    throw new Error('encode ' + n.t);
  };
  return { n: B.nodes.map(code), m: B.tn.slice(), q: cons || [], ...(B.lat ? { L: B.lat } : {}) };
}
const decode = core.decodeRecipe;
/* replay the recipe with the tool's evaluate() and compare with the catalog set */
function verify(rec, job) {
  const nodes = decode(rec), V = evaluate(nodes);
  const P = rec.m.map(id => V.get(id));
  if (P.some(p => !p || p.x === undefined)) return 'a point of the set is not reached';
  const T = []; for (let k = 0; k < job.flat.length; k += 2) T.push({ x: job.flat[k], y: job.flat[k + 1] });
  if (!sameCanon(canon(P), canon(T))) return 'the replayed set differs from the catalog set';
  const dc = distClasses(P);
  if (dc.cls.length !== job.d) return 'distance count ' + dc.cls.length + ' instead of ' + job.d;
  return null;
}
function finish(B, job, level, E) {
  if (!E) E = { cons: [], eqs: 0, unk: 0, free: 0 };
  const rec = encode(B, E.cons);
  const bad = verify(rec, job);
  if (bad) return { ok: false, why: bad };
  const steps = B.nodes.filter(n => n.t !== 'base' && n.t !== 'pv' && n.t !== 'pc').length;
  return { ok: true, level, eqs: E.eqs, unk: E.unk, free: E.free, steps, rec };
}

module.exports = { construct, decode, verify, encode, GAUSS };
