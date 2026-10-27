"""Step 2: random points on every positive-dimensional family of 5-point d-distance sets (saturated ideals in
sat5d.json), via random linear slices solved with msolve exactly as strat.py does; symmetry at each point.
Usage: python3 step2_sample.py d [nslices]   -> writes symtest/samples_5d.json (incremental, resumable)"""
import json, sys, os, re, math, random, subprocess, itertools, time, collections
from os import getpid
from fractions import Fraction
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
from symmetry import classify
PID = str(getpid())
MS = '/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
d = int(sys.argv[1]); NSL = int(sys.argv[2]) if len(sys.argv) > 2 else 4; n = 5
XV = ['x1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4']; SV = [f's{k}' for k in range(1, d)]; V = XV + SV
SYM = {v: sp.Symbol(v) for v in V}
TMP = os.path.join(HERE, 'tmp'); os.makedirs(TMP, exist_ok=True)

def run_msolve(eqs, vars_):
    fi, fo = os.path.join(TMP, f'ms_in_{PID}.txt'), os.path.join(TMP, f'ms_out_{PID}.txt')
    open(fi, 'w').write(','.join(vars_) + '\n0\n' + ',\n'.join(eqs) + '\n')
    subprocess.run([MS, '-f', fi, '-o', fo, '-p', '200'], capture_output=True, text=True, timeout=900)
    return open(fo).read()

def parse_real(out):
    s = out.replace('\n', '').replace(' ', '')
    if s.startswith('[-1]'): return 'empty', []
    if s.startswith('[1,'): return 'posdim', []
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s.rstrip(':'))
    if not tokens or tokens[0] != '[': return 'error', []
    def build(i):
        if tokens[i] == '[':
            lst = []; i += 1
            while tokens[i] != ']':
                v, i = build(i); lst.append(v)
                if tokens[i] == ',': i += 1
            return lst, i + 1
        tk = tokens[i]
        if '^' in tk:
            num, den = tk.split('/'); base, ex = den.split('^'); return Fraction(int(num), int(base) ** int(ex)), i + 1
        return Fraction(tk), i + 1
    tree, _ = build(0)
    sols = tree[1][1] if isinstance(tree[1], list) and len(tree[1]) > 1 else []
    return 'zerodim', [[float((lo + hi) / 2) for lo, hi in sol] for sol in sols]

def expand(expr): return str(sp.expand(sp.sympify(expr, locals=SYM))).replace('**', '^').replace(' ', '')
def pts_of(v): return np.array([[0, 0], [v[0], 0], [v[1], v[2]], [v[3], v[4]], [v[5], v[6]]], float)
def svals(v): return [1.0] + list(v[7:7 + d - 1])
def margin(v):
    s = svals(v); sc = max(abs(x) for x in s)
    return min(min(abs(a - b) for a, b in itertools.combinations(s, 2)), min(s[1:])) / sc

SAT = json.load(open(os.path.join(ROOT, f'sat5{d}.json')))
R = json.load(open(os.path.join(ROOT, f'cert5{d}.json'))); cols = {r['idx']: r['col'] for r in R}
outfn = os.path.join(HERE, f'samples_5{d}.json')
out = json.load(open(outfn)) if os.path.exists(outfn) else {}
t0 = time.time()
for k, v in sorted(SAT.items(), key=lambda kv: int(kv[0])):
    if v.get('timeout') or k in out: continue
    idx = int(k); dim = v['dim']; gens = list(v['gens'])
    if 'x1-1' not in gens: gens.append('x1-1')
    random.seed(1000 * d + idx)
    rec = {'idx': idx, 'col': cols[idx], 'dim': dim, 'samples': [], 'slices': 0, 'empty_slices': 0}
    tries = 0
    while tries < 3 * NSL and (rec['slices'] < NSL or len(rec['samples']) < 3) and tries < 12:
        tries += 1
        lins = ['+'.join(f'({random.randint(-5, 5)})*{x}' for x in XV) + f'-({random.randint(-3, 3)})' for _ in range(dim)]
        st, md = parse_real(run_msolve(gens + [expand(l) for l in lins], V))
        if st != 'zerodim': rec.setdefault('odd', []).append(st); continue
        rec['slices'] += 1
        good = 0
        for w in md:
            if margin(w) <= 1e-6: continue
            P = pts_of(w); c = classify(P)
            rec['samples'].append({'v': w, 'type': c['type'], 'order': c['order'], 'margin': margin(w)}); good += 1
        if good == 0: rec['empty_slices'] += 1
    out[k] = rec
    json.dump(out, open(outfn + '.tmp', 'w')); os.replace(outfn + '.tmp', outfn)
    types = collections.Counter(s['type'] for s in rec['samples'])
    print(f'd={d} type {idx} dim {dim}: {len(rec["samples"])} samples from {rec["slices"]} slices {dict(types)} {time.time()-t0:.0f}s', flush=True)
print('done', len(out))
