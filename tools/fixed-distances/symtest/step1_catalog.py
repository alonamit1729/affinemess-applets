"""Step 1: symmetry of the sample configuration of every family item in applet_data.json."""
import json, collections, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry import classify, e_dim
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
cats = D['cats']
rows = []
for key, c in cats.items():
    n, d = map(int, key.split(','))
    for idx, it in enumerate(c['items']):
        if it['kind'] != 'family': continue
        P = [[float(x), float(y)] for x, y in it['coords']]
        r = classify(P)
        rows.append({'cell': key, 'n': n, 'd': d, 'e': e_dim(n, d), 'idx': idx, 'dim': it['dim'], 'label': it['label'],
                     'order': r['order'], 'type': r['type'], 'coords': it['coords'], 'cls': it['cls']})
json.dump(rows, open(os.path.join(HERE, 'step1_families.json'), 'w'), indent=1)
# table by cell
bycell = collections.OrderedDict()
for r in rows:
    b = bycell.setdefault(r['cell'], {'n': r['n'], 'd': r['d'], 'e': r['e'], 'fam': 0, 'sym': 0, 'types': collections.Counter(), 'dims': collections.Counter()})
    b['fam'] += 1; b['dims'][r['dim']] += 1
    if r['order'] > 1: b['sym'] += 1
    b['types'][r['type']] += 1
lines = ['| cell (n,d) | e(n,d) | families (by dim) | sample has nontrivial symmetry | symmetry types of samples |', '|---|---|---|---|---|']
for key, b in sorted(bycell.items(), key=lambda kv: (kv[1]['n'], kv[1]['d'])):
    dims = ', '.join(f'{v}×dim{k}' for k, v in sorted(b['dims'].items()))
    types = ', '.join(f'{t}:{c}' for t, c in sorted(b['types'].items(), key=lambda kv: -kv[1]))
    lines.append(f"| ({b['n']},{b['d']}) | {b['e']} | {b['fam']} ({dims}) | {b['sym']} | {types} |")
open(os.path.join(HERE, 'step1_table.md'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
print()
print('families with e<=0 and trivial sample symmetry:')
for r in rows:
    if r['e'] <= 0 and r['order'] == 1: print(' ', r['cell'], r['idx'], r['dim'], r['label'][:100])
