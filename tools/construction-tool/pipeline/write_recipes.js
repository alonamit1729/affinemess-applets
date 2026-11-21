'use strict';
/* node write_recipes.js <jobs dir> <out dir>
   -> ../recipes/<n>_<d>/<k>.json, cards k·1000 … k·1000+999 of the cell:
      { cell, labels: [...], items: { num: [kind, dim, level, eqs, unknowns, label index, recipe, verdict] } }
   kind 'rigid' | 'family'; level 1 lattice, 2 straightedge and compass, 3 regular polygon, 4 equalizer,
   0 none found (the recipe then just places the points). A file is rewritten only when its content changes. */
const fs = require('fs'), path = require('path');
const [jobsDir, outDir] = process.argv.slice(2);
const DEST = path.join(__dirname, '..', 'recipes');
const verdicts = {};
for (const [cell, num, v] of JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'verdicts.json')))) if (v) verdicts[cell + '#' + num] = v;
let files = 0, written = 0, bytes = 0; const tally = {};
for (const f of fs.readdirSync(jobsDir).filter(f => f.endsWith('.json')).sort()) {
  const cell = f.replace('.json', '').replace('_', ',');
  const jobs = JSON.parse(fs.readFileSync(path.join(jobsDir, f)));
  const res = fs.existsSync(path.join(outDir, f)) ? JSON.parse(fs.readFileSync(path.join(outDir, f))) : {};
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
    const v = verdicts[cell + '#' + num]; if (v) item.push(v);
    C.items[num] = item;
    const t = (kind === 'rigid' ? 'R' : 'F') + item[2]; tally[t] = (tally[t] || 0) + 1;
  }
  const dir = path.join(DEST, f.replace('.json', '')); fs.mkdirSync(dir, { recursive: true });
  for (const [k, C] of chunks) {
    const fn = path.join(dir, k + '.json'), txt = JSON.stringify(C);
    files++; bytes += txt.length;
    if (!fs.existsSync(fn) || fs.readFileSync(fn, 'utf8') !== txt) { fs.writeFileSync(fn, txt); written++; }
  }
}
console.log(files, 'files,', written, 'written,', (bytes / 1e6).toFixed(1), 'MB;', JSON.stringify(tally));
