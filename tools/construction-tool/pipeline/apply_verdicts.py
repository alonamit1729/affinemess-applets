"""python3 apply_verdicts.py   — fold ../verdicts.json into what the tool's edit mode reads:
the verdicts embedded in ../catalog.json (small cells), and catalogs/<n>_<d>.verdicts.json {num: verdict} for the big
cells, whose fingerprint files (catalogs/<n>_<d>.json, up to 8 MB) are left alone. Files are written only if changed."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
V = {(c, n): v for c, n, v in json.load(open(os.path.join(ROOT, 'verdicts.json')))}
fn = os.path.join(ROOT, 'catalog.json'); C = json.load(open(fn)); k = 0
for it in C['items']:
    v = V.get((it[0], it[3]))
    if v and it[4] != v: it[4] = v; k += 1
if k: open(fn, 'w').write(json.dumps(C, separators=(',', ':')))
print('catalog.json:', k, 'verdicts updated')
for cell, c in C['cells'].items():
    if len(c) > 3 and c[3] == 1:
        D = {str(n): v for (cc, n), v in V.items() if cc == cell and v}
        out = os.path.join(ROOT, 'catalogs', cell.replace(',', '_') + '.verdicts.json'); txt = json.dumps(D, separators=(',', ':'))
        if not os.path.exists(out) or open(out).read() != txt: open(out, 'w').write(txt); print('wrote', out, len(D))
