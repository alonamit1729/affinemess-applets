"""Match the applet's (n,d) catalog items (applet_data.json, kind 'rigid') to the candidate patterns of
gen/levels/L{n}.json with <= d colours, by nauty canonical form (gen/patlib.py).  Usage: python3 g7fin/catalog_match.py"""
import sys, json, math, itertools, collections
sys.path.insert(0, 'g7fin/gen')
import patlib
A = json.load(open('applet_data.json'))
out = {}
for n, d in [(7, 4), (8, 4), (9, 4), (8, 5), (9, 5), (10, 5)]:
    L = json.load(open(f'g7fin/gen/levels/L{n}.json'))['patterns']
    canon = {}
    for idx, p in enumerate(L):
        if max(p) + 1 <= d:
            canon[patlib.canon(n, tuple(p))] = idx
    items = A['cats'][f'{n},{d}']['items']
    hits = collections.Counter(); miss = []
    for it in items:
        if it['kind'] != 'rigid': continue
        pts = [(float(x), float(y)) for x, y in it['coords']]
        assert len(pts) == n
        vals = []; col = []
        for i, j in itertools.combinations(range(n), 2):
            x = math.dist(pts[i], pts[j])
            for c, v in enumerate(vals):
                if abs(x - v) < 1e-6 * max(1, x): col.append(c); break
            else: vals.append(x); col.append(len(vals) - 1)
        k = len(vals)
        c = patlib.canon(n, tuple(patlib.rgf(col)))
        if c in canon: hits[canon[c]] += 1
        else: miss.append((it['label'], k))
    covered = sorted(hits)
    print(f'({n},{d}): {len(items)} catalog items, {len(canon)} patterns; items matched {sum(hits.values())}, '
          f'distinct patterns hit {len(hits)}, unmatched items {miss}')
    multi = {f'L{n}#{i}': c for i, c in hits.items() if c > 1}
    if multi: print('   patterns with several catalog realizations:', multi)
    nothit = [f'L{n}#{i}' for i in canon.values() if i not in hits]
    print('   patterns with no catalog item:', nothit)
    out[f'{n},{d}'] = {'items': len(items), 'patterns': len(canon), 'matched': sum(hits.values()),
                       'distinct_hit': len(hits), 'unmatched': miss, 'multi': multi, 'nothit': nothit}
json.dump(out, open('g7fin/catalog_match.json', 'w'), indent=1)
