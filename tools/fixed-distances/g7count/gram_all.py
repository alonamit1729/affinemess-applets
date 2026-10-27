"""Gram route for all 228 patterns (cross-check of the coordinate route, and the source for the patterns whose
coordinate system is too heavy for msolve): real10g/pat_<idx>.json with the s-values, eigenvalues, PSD flags and
reconstructed coordinates (100 digits).  Usage: python3 gram_all.py [idx ...]"""
import json, os, sys, time
from gram_util import *
L = json.load(open(os.path.join(SCR,'g7fin','gen','levels','L10.json')))['patterns']
sys.path.insert(0, os.path.join(SCR, 'g7fin', 'gen')); import patlib
idxs = [int(a) for a in sys.argv[1:]] or range(len(L))
for idx in idxs:
    fn = f'real10g/pat_{idx}.json'
    if os.path.exists(fn): continue
    col = tuple(L[idx])
    V, sols, dt = gram_realize(10, col)
    rec = {'idx': idx, 'col': list(col), 'vars': V, 'time': dt, 'source': 'gram'}
    if isinstance(sols, str):
        rec['status'] = sols; rec['sets'] = []
    else:
        rec['status'] = 'real' if sols else 'noreal'
        R = sets_from_gram(10, col, V, sols)
        rec['nsol'] = len(sols); rec['npsd'] = sum(r['psd'] for r in R)
        rec['eig'] = [[mp.nstr(x, 30) for x in r['eig'][-3:]] for r in R]
        sets = []
        for r in R:
            if not r['psd']: continue
            P = coords_from_gram(r['G'], 10)
            assert patlib.canon(10, pattern(P)) == patlib.canon(10, col)
            assert pattern(P) == col  # same labelling: edge colours must agree exactly
            sets.append([[mp.nstr(p[0], 100), mp.nstr(p[1], 100)] for p in P])
        rec['sets'] = sets
    tmp = fn + f'.tmp{os.getpid()}'
    json.dump(rec, open(tmp, 'w')); os.replace(tmp, fn)
    print(idx, rec['status'], rec.get('nsol'), 'psd', rec.get('npsd'), round(dt, 1), flush=True)
