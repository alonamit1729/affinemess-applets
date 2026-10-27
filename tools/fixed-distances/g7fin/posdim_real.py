"""Real points on the positive-dimensional patterns found by pat_rigid_d.py: random hyperplane sections of the
coordinate system (gen/realize.py conventions) solved by msolve real root isolation.  Prints, per pattern and slice,
the number of real solutions and one real configuration (and the distance values), to see whether the
one-dimensional component is real.  Usage: python3 g7fin/posdim_real.py"""
import sys, os, json, random, itertools, math
sys.path.insert(0, 'g7fin/gen')
import realize
from realize import system_coords, msolve_real, real_status
from realcheck import random_linear

def parse_real(out):
    """msolve -p output -> list of solutions (midpoints of the isolating boxes), using realize.real_status's grammar."""
    import re
    from fractions import Fraction
    s = out.replace('\n', '').replace(' ', '').rstrip(':')
    if s.startswith('[-1'): return 'empty'
    if s.startswith('[1,'): return 'posdim'
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s)
    def build(i):
        if tokens[i] == '[':
            lst = []; i += 1
            while tokens[i] != ']':
                v, i = build(i); lst.append(v)
                if tokens[i] == ',': i += 1
            return lst, i + 1
        tk = tokens[i]
        if '^' in tk:
            num, den = tk.split('/'); base, ex = den.split('^')
            return Fraction(int(num), int(base) ** int(ex)), i + 1
        return Fraction(tk), i + 1
    tree, _ = build(0)
    sols = tree[1][1]
    return [[float((a + b) / 2) for a, b in sol] for sol in sols]

targets = {'9d5': ['L9#124'], '8d5': ['L8#398', 'L8#1081']}
res = {}
for f, names in targets.items():
    n = int(f[0])
    P = {p['name']: p for p in json.load(open(f'g7fin/patterns_{f}.json'))}
    pairs = list(itertools.combinations(range(n), 2))
    for nm in names:
        col = tuple(P[nm]['col'][f'{a},{b}'] for a, b in pairs)
        k = max(col) + 1
        V, polys = system_coords(n, col)
        XV = [v for v in V if v[0] in 'xy']
        print('==', nm, 'k', k, 'vars', V)
        random.seed(1)
        for sl in range(6):
            lin = random_linear(XV)
            out = msolve_real(V, polys + [lin], precision=120, timeout=600, tag=f'_{nm.replace("#","_")}_{sl}')
            st = real_status(out)
            print('  slice', sl, lin, '->', st)
            if st == 'real':
                sols = parse_real(out)
                print('   #real solutions in slice:', len(sols))
                for sol in sols[:2]:
                    d = dict(zip(V, sol))
                    pts = [(0.0, 0.0), (d['x1'], 0.0)] + [(d[f'x{i}'], d[f'y{i}']) for i in range(2, n)]
                    svals = [1.0] + [d[f's{c}'] for c in range(1, k)]
                    dist = {c: [] for c in range(k)}
                    for (a, b), c in zip(pairs, col):
                        dist[c].append(math.dist(pts[a], pts[b]))
                    print('   pts', [(round(x, 5), round(y, 5)) for x, y in pts])
                    print('   sq values', [round(v, 6) for v in svals], 'check spread per class',
                          [round(max(dist[c]) - min(dist[c]), 9) for c in range(k)])
                res.setdefault(nm, []).append({'slice': lin, 'nreal': len(sols), 'sol': sols[0] if sols else None})
json.dump(res, open('g7fin/posdim_real.json', 'w'), indent=1)
