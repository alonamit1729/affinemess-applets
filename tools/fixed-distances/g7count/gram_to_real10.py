"""Convert the Gram-route realizations (real10g/pat_<i>.json, 100-digit coordinates with p0 at the origin) into the
input format of complete10.py (real10/pat_<i>.json with 'vars' and 'sols' = one interval [a, a] per variable), rotating
each set so that p1 lies on the positive x-axis.  Usage: python3 gram_to_real10.py [outdir]   (default real10)"""
import json, glob, os, sys
import mpmath as mp
mp.mp.dps = 110
out = sys.argv[1] if len(sys.argv) > 1 else 'real10'
os.makedirs(out, exist_ok=True)
VARS = ['x1'] + [v for i in range(2, 10) for v in (f'x{i}', f'y{i}')] + [f's{c}' for c in range(1, 6)] + ['t']
for f in sorted(glob.glob('real10g/pat_*.json')):
    d = json.load(open(f))
    sols = []
    for S in d['sets']:
        P = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in S]
        z = P[1] - P[0]; rot = abs(z) / z
        Q = [(p - P[0]) * rot for p in P]
        vals = {'x1': Q[1].real}
        for i in range(2, 10): vals[f'x{i}'] = Q[i].real; vals[f'y{i}'] = Q[i].imag
        for c in range(1, 6): vals[f's{c}'] = mp.mpf(0)   # colour values are not used downstream
        vals['t'] = mp.mpf(0)
        sols.append([[mp.nstr(vals[v], 100), mp.nstr(vals[v], 100)] for v in VARS])
    rec = {'idx': d['idx'], 'col': d['col'], 'vars': VARS, 'status': d['status'], 'source': 'gram', 'sols': sols}
    fn = os.path.join(out, f"pat_{d['idx']}.json"); tmp = fn + f'.tmp{os.getpid()}'
    json.dump(rec, open(tmp, 'w')); os.replace(tmp, fn)
print('converted', len(glob.glob('real10g/pat_*.json')), 'patterns into', out)
