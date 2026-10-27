"""Extract from gen/levels/L{n}.json the candidate patterns with at most k colours into g7fin/patterns_{n}d{k}.json
(pat_rigid input format).  Usage: python3 g7fin/extract_d.py n k"""
import sys, json, itertools, collections
n, k = int(sys.argv[1]), int(sys.argv[2])
L = json.load(open(f'g7fin/gen/levels/L{n}.json'))
assert L['n'] == n and L['count'] == len(L['patterns'])
pairs = list(itertools.combinations(range(n), 2))
out = []
for idx, p in enumerate(L['patterns']):
    assert len(p) == len(pairs)
    kk = max(p) + 1
    if kk > k: continue
    assert sorted(set(p)) == list(range(kk)), (idx, p)   # rgf: colours 0..kk-1 all used
    col = {f'{a},{b}': c for (a, b), c in zip(pairs, p)}
    out.append({'name': f'L{n}#{idx}', 'k': kk, 'sources': [f'L{n}.json'], 'col': col})
fn = f'g7fin/patterns_{n}d{k}.json'
json.dump(out, open(fn, 'w'))
print(fn, len(out), collections.Counter(o['k'] for o in out))
