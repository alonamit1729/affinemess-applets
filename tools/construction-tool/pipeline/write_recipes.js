'use strict';
/* node write_recipes.js <jobs dir> <out dir>
   -> ../recipes/<n>_<d>/<k>.json, cards k·1000 … k·1000+999 of the cell:
      { cell, labels: [...], items: { num: [kind, dim, level, eqs, unknowns, label index, recipe, verdict or 0, animation?] } }
   animation (families): <out dir>/anim_<n>_<d>.json from animate.js
   kind 'rigid' | 'family'; level 1 lattice, 2 straightedge and compass, 3 regular polygon, 4 equalizer,
   0 none found (the recipe then just places the points). A file is rewritten only when its content changes. */
const fs = require('fs'), path = require('path');
const core = require('./core.js');
/* index.json per cell, for the catalog's filters: one 10-character code per card (card k at offset 10·(k−1)):
   0 kind: '0' rigid, '1'..'4' family dimension   1 level '0'..'4'   2 regular polygon used: char(48 + m), '0' none
   3 lattice: 't' triangular, 's' square, '0' none   4 helper points: '1' / '0'
   5 verdict: 'r' straightedge and compass, 'p' a regular polygon, 't' trisection, 'q' S4 (solvable, one real root),
     'c' a cube root, 'u' not solvable, 'b' beyond straightedge and compass (which tool undecided), '?' undecided,
     '-' not analysed
   6 rotations: char(48 + order)   7 mirror: '1' / '0'   8-9 steps (as the tool counts them), two digits */
const VCODE = v => !v ? '-' : v === 'rc' ? 'r' : v[0] === 'p' ? 'p' : v === 'tri' ? 't' : v === 's4' ? 'q' : v === 'cbrt' ? 'c' : v === 'nrc' ? 'b' : v === '?' ? '?' : 'u';
function symmetry(P) {
  const n = P.length, c = P.reduce((s, p) => ({ x: s.x + p.x / n, y: s.y + p.y / n }), { x: 0, y: 0 });
  const Q = P.map(p => ({ x: p.x - c.x, y: p.y - c.y })), R0 = Math.max(...Q.map(q => Math.hypot(q.x, q.y)));
  const tol = 1e-6 * R0, has = q => Q.some(p => Math.hypot(p.x - q.x, p.y - q.y) < tol);
  const p0 = Q.reduce((a, b) => Math.hypot(b.x, b.y) > Math.hypot(a.x, a.y) ? b : a);
  if (Math.hypot(p0.x, p0.y) < tol) return [1, 0];
  let rot = 0, mir = 0;
  for (const q of Q) {
    if (Math.abs(Math.hypot(q.x, q.y) - Math.hypot(p0.x, p0.y)) > tol) continue;
    const t = Math.atan2(q.y, q.x) - Math.atan2(p0.y, p0.x), co = Math.cos(t), si = Math.sin(t);
    if (Q.every(p => has({ x: co * p.x - si * p.y, y: si * p.x + co * p.y }))) rot++;
    const a = (Math.atan2(q.y, q.x) + Math.atan2(p0.y, p0.x)), c2 = Math.cos(a), s2 = Math.sin(a);   // the mirror taking p0 to q
    if (Q.every(p => has({ x: c2 * p.x + s2 * p.y, y: s2 * p.x - c2 * p.y }))) mir = 1;
  }
  return [Math.max(1, rot), mir];
}
function indexCode(item) {
  const [kind, dim, level, eqs, unk, li, rec, v] = item;
  const nodes = core.decodeRecipe(rec), V = core.evaluate(nodes), P = rec.m.map(id => V.get(id));
  const [rot, mir] = P.every(p => p && p.x !== undefined) ? symmetry(P) : [1, 0];
  const poly = nodes.find(n => n.t === 'poly'), marked = new Set(rec.m);
  const helpers = nodes.some(n => n.t === 'x' && !marked.has(n.id));
  const steps = (nodes.some(n => n.t === 'lat') || rec.L && level === 1 ? 1 : 0) + nodes.filter(n => n.t === 'poly' || n.t === 'x' || n.t === 'glide' || n.t === 'free').length + (rec.q || []).length;
  return (kind === 'rigid' ? '0' : String(Math.min(9, dim))) + level + String.fromCharCode(48 + (poly ? poly.m : 0)) + (rec.L ? rec.L[0] : '0') +
    (helpers ? '1' : '0') + VCODE(v) + String.fromCharCode(48 + Math.min(70, rot)) + mir + String(Math.min(99, steps)).padStart(2, '0');
}
const [jobsDir, outDir] = process.argv.slice(2);
const DEST = path.join(__dirname, '..', 'recipes');
const verdicts = {};
for (const [cell, num, v] of JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'verdicts.json')))) if (v) verdicts[cell + '#' + num] = v;
// verdicts from galois_probe.py + classify.py, not yet folded into verdicts.json
const extra = path.join(outDir, 'verdicts_new.json');
if (fs.existsSync(extra)) for (const [cell, num, v] of JSON.parse(fs.readFileSync(extra))) if (v && !verdicts[cell + '#' + num]) verdicts[cell + '#' + num] = v;
let files = 0, written = 0, bytes = 0; const tally = {};
for (const f of fs.readdirSync(jobsDir).filter(f => f.endsWith('.json')).sort()) {
  const cell = f.replace('.json', '').replace('_', ',');
  const jobs = JSON.parse(fs.readFileSync(path.join(jobsDir, f)));
  const res = fs.existsSync(path.join(outDir, f)) ? JSON.parse(fs.readFileSync(path.join(outDir, f))) : {};
  const anim = fs.existsSync(path.join(outDir, 'anim_' + f)) ? JSON.parse(fs.readFileSync(path.join(outDir, 'anim_' + f))) : {};
  const chunks = new Map();
  for (const [num, kind, dim, flat, cls, label] of jobs) {
    const k = Math.floor(num / 1000);
    if (!chunks.has(k)) chunks.set(k, { cell, labels: [], items: {} });
    const C = chunks.get(k);
    let li = C.labels.indexOf(label); if (li < 0) { li = C.labels.length; C.labels.push(label); }
    const r = res[num];
    let item;
    if (r && r.ok) {
      const lat = r.rec.n.find(a => a[0] === 't' || a[0] === 's');      // the lattice a start uses, shown with the construction
      if (lat && !r.rec.L) r.rec.L = lat[0] === 't' ? 'tri' : 'sq';
      item = [kind, dim, r.level, r.eqs, r.unk, li, r.rec];
    }
    else {   // no construction: the points by their coordinates
      const n = [], m = [];
      for (let q = 0; q < flat.length; q += 2) { m.push(n.length); n.push(['b', +flat[q].toPrecision(13), +flat[q + 1].toPrecision(13)]); }
      item = [kind, dim, 0, 0, 0, li, { n, m, q: [] }];
    }
    item.push(verdicts[cell + '#' + num] || 0);
    if (r && r.ok && anim[num]) item.push(anim[num]);
    C.items[num] = item;
    const t = (kind === 'rigid' ? 'R' : 'F') + item[2]; tally[t] = (tally[t] || 0) + 1;
  }
  const dir = path.join(DEST, f.replace('.json', '')); fs.mkdirSync(dir, { recursive: true });
  // the cell's index for the catalog's filters
  const codes = jobs.map(([num]) => indexCode(chunks.get(Math.floor(num / 1000)).items[num]));
  if (jobs.some(([num], i) => num !== i + 1)) throw new Error('card numbers are not 1..N in ' + f);
  { const fn = path.join(dir, 'index.json'), txt = JSON.stringify({ cell, n: jobs.length, s: codes.join('') });
    files++; bytes += txt.length; if (!fs.existsSync(fn) || fs.readFileSync(fn, 'utf8') !== txt) { fs.writeFileSync(fn, txt); written++; } }
  for (const [k, C] of chunks) {
    const fn = path.join(dir, k + '.json'), txt = JSON.stringify(C);
    files++; bytes += txt.length;
    if (!fs.existsSync(fn) || fs.readFileSync(fn, 'utf8') !== txt) { fs.writeFileSync(fn, txt); written++; }
  }
}
console.log(files, 'files,', written, 'written,', (bytes / 1e6).toFixed(1), 'MB;', JSON.stringify(tally));
