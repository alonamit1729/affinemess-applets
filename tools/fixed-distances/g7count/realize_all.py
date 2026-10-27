"""Real realizations of all 228 ten-point patterns (msolve real root isolation, -p 400).
Usage: python3 realize_all.py <worker> <nworkers>.  Checkpoint: real10/pat_<idx>.json (atomic)."""
import sys, json, os, time
from msolve_util import *
L = json.load(open(os.path.join(SCR,'g7fin','gen','levels','L10.json')))['patterns']
w, W = int(sys.argv[1]), int(sys.argv[2])
outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'real10')
for idx in range(w, len(L), W):
    fn = os.path.join(outdir, f'pat_{idx}.json')
    if os.path.exists(fn): continue
    col = tuple(L[idx])
    V, sols, dt = realize_pattern(10, col, precision=400, timeout=14000)
    rec = {'idx': idx, 'col': list(col), 'vars': V, 'time': dt}
    if isinstance(sols, str):
        rec['status'] = sols
    else:
        rec['status'] = 'real' if sols else 'noreal'
        rec['sols'] = [[[str(a), str(b)] for a, b in sol] for sol in sols]
    tmp = fn + f'.tmp{os.getpid()}'
    json.dump(rec, open(tmp, 'w')); os.replace(tmp, fn)
    print(idx, rec['status'], len(sols) if isinstance(sols, list) else '', round(dt,1), flush=True)
print('worker', w, 'done')
