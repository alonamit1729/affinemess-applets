'use strict';
/* node run.js <jobs dir> <out dir> <cell> [<cell> ...] [--from=K --to=K] [--threads=2]
   jobs/<n>_<d>.json: [num, kind, dim, flat coords, cls, label] (made by export_jobs.py)
   out/<n>_<d>.json:  {num: result} merged into what is already there (so a cell can be done in pieces) */
const fs = require('fs'), path = require('path');
const { Worker, isMainThread, parentPort, workerData } = require('worker_threads');

if (!isMainThread) {
  const { construct } = require('./construct.js');
  parentPort.on('message', jobs => {
    const out = [];
    for (const j of jobs) {
      const t0 = Date.now(); let r;
      try { r = construct(j); } catch (e) { r = { ok: false, why: 'error: ' + e.message }; }
      r.ms = Date.now() - t0; r.num = j.num; out.push(r);
    }
    parentPort.postMessage(out);
  });
  return;
}

const args = process.argv.slice(2), opt = {};
const pos = args.filter(a => { const m = a.match(/^--(\w+)=(.*)$/); if (m) { opt[m[1]] = m[2]; return false; } return true; });
const [jobsDir, outDir, ...cells] = pos;
const threads = +(opt.threads || 2), from = +(opt.from || 1), to = +(opt.to || 1e9), budget = +(opt.seconds || 1e9) * 1000;
fs.mkdirSync(outDir, { recursive: true });
const T0 = Date.now();

(async () => {
  for (const cell of cells) {
    const [n, d] = cell.split(/[,_]/).map(Number);
    const raw = JSON.parse(fs.readFileSync(path.join(jobsDir, n + '_' + d + '.json')));
    const outFile = path.join(outDir, n + '_' + d + '.json');
    const prev = fs.existsSync(outFile) ? JSON.parse(fs.readFileSync(outFile)) : {};
    const jobs = raw.filter(([num]) => num >= from && num <= to && !(prev[num] && !(opt.redo === 'all' || (opt.redo === 'l4' && (!prev[num].ok || prev[num].level === 4)))))
      .map(([num, kind, dim, flat, cls, label, frame]) => ({ num, kind, dim, flat, cls, d, frame }));
    if (!jobs.length) { console.log(cell, 'nothing to do'); continue; }
    const W = [...Array(threads)].map(() => new Worker(__filename));
    let next = 0, done = 0; const CH = 8;
    const save = () => fs.writeFileSync(outFile, JSON.stringify(prev));
    await new Promise(res => {
      let live = W.length;
      const feed = w => {
        if (next >= jobs.length || Date.now() - T0 > budget) { w.terminate(); if (!--live) res(); return; }
        const part = jobs.slice(next, next + CH); next += CH; w.postMessage(part);
      };
      W.forEach(w => { w.on('message', out => { for (const r of out) { prev[r.num] = r; done++; } if (done % 400 < CH) { save(); console.log(cell, done, '/', jobs.length, ((Date.now() - T0) / 1000).toFixed(0) + 's'); } feed(w); }); feed(w); });
    });
    save();
    const S = {}; for (const r of Object.values(prev)) { const k = r.ok ? 'L' + r.level : 'fail'; S[k] = (S[k] || 0) + 1; }
    console.log(cell, 'done', done, JSON.stringify(S), ((Date.now() - T0) / 1000).toFixed(0) + 's');
    if (Date.now() - T0 > budget) { console.log('time budget reached'); break; }
  }
})();
