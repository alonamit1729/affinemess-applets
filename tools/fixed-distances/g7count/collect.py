"""Collect the complete pass: statistics and the similarity classes of all 16-point sets found."""
import json, glob, os
from extend import *
SCR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
recs = [json.load(open(f)) for f in sorted(glob.glob('ext10/pat_*.json'), key=lambda f: int(f.split('_')[-1][:-5]))]
print('patterns processed:', len(recs))
from collections import Counter
print('status:', Counter(r['status'] for r in recs))
print('k:', Counter(r['k'] for r in recs))
print('labelled real solutions: total', sum(r['nsol'] for r in recs), 'distribution', sorted(Counter(r['nsol'] for r in recs).items()))
print('distinct realizations (similarity classes): total', sum(r['nreal'] for r in recs), 'per pattern', sorted(Counter(r['nreal'] for r in recs).items()))
print('patterns with no real realization:', [r['idx'] for r in recs if r['nreal'] == 0])
cores = [c for r in recs for c in r['cores']]
print('cores extended:', len(cores), ' max survivors', max(c['ncand'] for c in cores), ' max dfs nodes', max(c['dfs_nodes'] for c in cores),
      ' total labelled 6-extensions', sum(c['nsol'] for c in cores), ' cores with >=1 extension', sum(1 for c in cores if c['nsol']))
print('time: msolve', round(sum(json.load(open(f))['time'] for f in glob.glob('real10/pat_*.json'))), 's; extension', round(sum(r['time'] for r in recs)), 's')
C = Classes()
K = json.load(open(os.path.join(SCR, 'r12', 'r12_subsets_7.json')))['16,7']
for r in recs:
    for ci, c in enumerate(r['cores']):
        for s in c['sets']:
            X = [(mp.mpf(x), mp.mpf(y)) for x, y in s['pts']]
            assert s['k'] == 7
            idx, new = C.add(X, {'first_pattern': r['idx'], 'core': ci})
            C.reps[idx][2] += s['copies'] - 1
print('DISTINCT 16-POINT 7-DISTANCE SETS:', len(C.reps))
out = []
for i, (X, info, hits) in enumerate(C.reps):
    known = [j for j, k in enumerate(K) if congruent(X, [(mp.mpf(x), mp.mpf(y)) for x, y in k['pts']], tol=mp.mpf('1e-7'))]
    ncores = sum(1 for r in recs for c in r['cores'] for s in c['sets'] if congruent(X, [(mp.mpf(x), mp.mpf(y)) for x, y in s['pts']]))
    print(' set', i, 'labelled hits', hits, 'from', ncores, 'cores; congruent to known set(s)', known, '; pattern colours', max(pattern(X)) + 1)
    out.append({'pts': to_str(X), 'pts_float': to_floats(X), 'hits': hits, 'ncores': ncores, 'known': known, 'pattern': pattern(X)})
json.dump({'npatterns': len(recs), 'sets': out}, open('complete_result.json', 'w'), indent=1)
