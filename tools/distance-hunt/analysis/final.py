import json, collections, sys
D = './'  # run from the folder holding rigid.json and the analysis outputs
F = json.load(open(D + 'out_missing.json')) + json.load(open(D + 'out_exact1.json')) + json.load(open(D + 'out_exact2.json'))
G = json.load(open(D + 'groups.json'))
try:
    G.update(json.load(open(D + 'groups_exact.json')))
except Exception:
    pass
CH = {(r['cell'], r['i']): r for r in json.load(open(D + 'chain_out.json'))} if len(sys.argv) > 1 else {}

def pow2(k): return k & (k - 1) == 0

def rkind(q, res):
    if res > 1e-200:
        return 'rc*'   # singular solution; identified separately as golden-ratio values
    if q['deg'] is None:
        return '??'
    if all(pow2(o) for o in q['orders']):
        return 'rc'
    if q['equal_deg']:
        return 'cyclo%d' % q['deg']
    if q['deg'] == 3:
        return 'tri' if q['real_roots'] == 3 else 'IMP:cube-root'
    g = G.get(str(q['poly'])) if q.get('poly') else None
    if g:
        name = g[0].split('.')[-1]
        if name in ('S5', 'S6', 'A5', 'A6', 'PGL2F5', 'PSL2F5'):
            return 'IMP:' + name
        if name == 'S4':
            return 'tri' if q['real_roots'] in (0, 4) else 'S4-'
        if name == 'S4xC2':
            return 'tri?' if q['real_roots'] == 6 else 'S4xC2-'
        return name
    if q['deg'] == 7 and any(o in q['orders'] for o in (10, 12)):
        return 'IMP:S7/A7'
    return 'deg%d?' % q['deg']

def verdict(ks):
    ks = set(ks) - {'rc', 'rc*'}
    if not ks:
        return '1 ruler+compass'
    if any(k.startswith('IMP') for k in ks):
        return '5 impossible (' + ','.join(sorted(k[4:] for k in ks if k.startswith('IMP'))) + ')'
    if any(k in ('S4-', 'S4xC2-') for k in ks):
        return '4 solvable, cubic with one real root (S4)'
    if all(k.startswith('cyclo') for k in ks):
        return '2 cyclotomic (needs m-gon)'
    return '3 trisection'

cells = collections.defaultdict(collections.Counter)
rows = []
for r in F:
    if 'ratios' not in r:
        r['ratios'] = []; r['res'] = 0
        print('ERR', r['cell'], r['i'], r['label'])
    ks = [rkind(q, r['res']) for q in r['ratios']]
    v = verdict(ks)
    ch = CH.get((r['cell'], r['i']))
    rows.append((r['cell'], r['i'], v, ks, ch, r['label']))
    cells[r['cell']][v] += 1

def key(c):
    a, b = c.split(','); return (int(a), int(b))
tot = collections.Counter()
for c in sorted(cells, key=key):
    print(c, dict(sorted(cells[c].items())))
    tot.update(cells[c])
print('TOTAL', dict(sorted(tot.items())))
json.dump(rows, open(D + 'final_rows.json', 'w'))
