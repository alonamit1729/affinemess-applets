"""Sanity checks on the complete pass: (1) every 10-point set of the R12+R12 catalog with <= 6 distances and every
10-point subset of the 68 eleven-point cores is congruent to one of the msolve realizations; (2) counts."""
import json, glob, itertools, os
from extend import *
from cores11 import all_sets
SCR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = Classes()
npat = 0
for f in sorted(glob.glob(os.environ.get('EXT10DIR','ext10')+'/pat_*.json')):
    d = json.load(open(f)); npat += 1
    for c in d['cores']:
        A = [(mp.mpf(x), mp.mpf(y)) for x, y in c['A']]
        idx, new = C.add(A, d['idx'])
        assert new, ('duplicate realization across patterns', d['idx'], C.reps[idx][1])
print('patterns processed', npat, ' distinct 10-point realizations', len(C.reps))
def lookup(P, tol):
    k = inv_key(P)
    # tolerant lookup: compare with all reps having same number of distances (slow but fine)
    for r in C.reps:
        if congruent(r[0], P, tol=tol): return True
    return False
r12 = json.load(open(os.path.join(SCR, 'r12', 'r12_subsets_7.json')))
miss = 0; tot = 0
for key in ('10,5', '10,6'):
    for s in r12[key]:
        P = [(mp.mpf(x), mp.mpf(y)) for x, y in s['pts']]; tot += 1
        if not lookup(P, mp.mpf('1e-6')): miss += 1; print('  r12 set missing', key)
print('r12 catalog 10-point sets:', tot, 'missing from realizations:', miss)
S = Classes(); miss = 0
for A, info in all_sets():
    for i in range(11):
        P = [p for j, p in enumerate(A) if j != i]
        idx, new = S.add(P)
        if new and not lookup(P, mp.mpf('1e-30')): miss += 1; print('  10-subset of core missing', info, i)
print('10-subsets of the 68 eleven-point cores:', len(S.reps), 'missing:', miss)
