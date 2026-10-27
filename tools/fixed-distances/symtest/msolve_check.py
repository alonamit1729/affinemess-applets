"""independent check of an n-point family item: solve the colouring system with x1=1 and dim random linear slices
with msolve; list real solutions with distinct nonzero distance values and their symmetry."""
import json, sys, os, re, random, itertools, subprocess
from os import getpid
from fractions import Fraction
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify
from step2_sample_lib import run_msolve, parse_real
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
cell = sys.argv[1]; idx = int(sys.argv[2]) if sys.argv[1] != 'super' else 0; nsl = int(sys.argv[3]) if len(sys.argv) > 3 else 3
if cell == 'super':
    S = json.load(open(os.path.join(HERE, 'superset_families.json'))); it = S[sys.argv[2]]; idx = sys.argv[2]; it['dim'] = 1
else: it = D['cats'][cell]['items'][idx]
n = it['n']; cls0 = it['cls']; d = max(cls0) + 1
rl = {}
for c in cls0: rl.setdefault(c, len(rl))
cls = [rl[c] for c in cls0]  # colour of pair (0,1) becomes class 0, whose squared length is normalised to 1
E = list(itertools.combinations(range(n), 2))
XV = ['x1'] + [f'{c}{i}' for i in range(2, n) for c in ('x', 'y')]; SV = [f's{k}' for k in range(1, d)]; V = XV + SV
SYM = {v: sp.Symbol(v) for v in V}
def expand(expr): return str(sp.expand(sp.sympify(expr, locals=SYM))).replace('**', '^').replace(' ', '')
def coord(i):
    if i == 0: return ('0', '0')
    if i == 1: return ('x1', '0')
    return (f'x{i}', f'y{i}')
Sv = ['1'] + SV
eqs = []
for k, (i, j) in enumerate(E):
    (a, b), (c, e) = coord(i), coord(j)
    eqs.append(expand(f'({a}-({c}))^2+({b}-({e}))^2-({Sv[cls[k]]})'))
eqs.append('x1-1')
def pts(v):
    P = [[0, 0], [v[0], 0]]
    for i in range(2, n): P.append([v[1 + 2 * (i - 2)], v[2 + 2 * (i - 2)]])
    return np.array(P)
def margin(v):
    s = [1.0] + list(v[len(XV):]); sc = max(abs(x) for x in s)
    return min(min(abs(a - b) for a, b in itertools.combinations(s, 2)), min(s[1:])) / sc
random.seed(hash(str(idx)) % 1000)
DIM = int(sys.argv[4]) if len(sys.argv) > 4 else it['dim']
for t in range(nsl):
    lins = ['+'.join(f'({random.randint(-5, 5)})*{x}' for x in XV) + f'-({random.randint(-3, 3)})' for _ in range(DIM)]
    st, md = parse_real(run_msolve(eqs + [expand(l) for l in lins], V))
    good = [(classify(pts(w))['type'], round(margin(w), 4)) for w in md if margin(w) > 1e-6]
    print(f'{cell} #{idx} slice {t}: status {st}, {len(md)} real solutions, non-degenerate: {good}', flush=True)
    if os.environ.get('SAVE'):
        fn = os.path.join(HERE, 'generic_samples.json'); G = json.load(open(fn)) if os.path.exists(fn) else {}
        G.setdefault(f'{cell}:{idx}', []).extend([{'pts': pts(w).tolist(), 'type': classify(pts(w))['type'], 'margin': margin(w)} for w in md if margin(w) > 1e-6])
        json.dump(G, open(fn, 'w'))
