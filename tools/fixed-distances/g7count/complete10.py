"""COMPLETE PASS: for every real realization A of each of the 228 ten-point patterns (real10/pat_<idx>.json from
realize_all.py), enumerate all 16-point supersets with <= 7 distances (extend.py, 6 extra points).
Usage: python3 complete10.py <worker> <nworkers>   (loops until every pattern has been processed; checkpoint
ext10/pat_<idx>.json)."""
import sys, os, json, time, glob, itertools
from fractions import Fraction
from extend import *
SCR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SCR, 'g7fin', 'gen'))
import patlib
HERE = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(SCR, 'g7fin', 'gen', 'levels', 'L10.json')))['patterns']
w, W = int(sys.argv[1]), int(sys.argv[2])
ONCE = len(sys.argv) > 3
os.makedirs(os.path.join(HERE, 'ext10'), exist_ok=True)


def realizations(rec):
    """distinct (up to similarity) 10-point sets from the msolve solutions; checks each one's pattern."""
    V = rec['vars']; col = tuple(rec['col'])
    canon = patlib.canon(10, col)
    C = Classes()
    for sol in rec.get('sols', []):
        vv = {}
        for name, (a, b) in zip(V, sol):
            fa, fb = Fraction(a), Fraction(b)
            m = (fa + fb) / 2
            vv[name] = mp.mpf(m.numerator) / mp.mpf(m.denominator)
        P = [(mp.mpf(0), mp.mpf(0)), (vv['x1'], mp.mpf(0))] + [(vv[f'x{i}'], vv[f'y{i}']) for i in range(2, 10)]
        pat = pattern(P)
        assert patlib.canon(10, pat) == canon, ('pattern mismatch', rec['idx'])
        C.add(P)
    return [r[0] for r in C.reps]


while True:
    todo = []
    for idx in range(w, len(L), W):
        fo = os.path.join(HERE, 'ext10', f'pat_{idx}.json')
        fi = os.path.join(HERE, 'real10', f'pat_{idx}.json')
        if os.path.exists(fo): continue
        if not os.path.exists(fi): todo.append(idx); continue
        rec = json.load(open(fi))
        t0 = time.time()
        reals = realizations(rec)
        k = max(rec['col']) + 1
        out = {'idx': idx, 'k': k, 'status': rec['status'], 'nsol': len(rec.get('sols', [])), 'nreal': len(reals), 'cores': []}
        for A in reals:
            A = normalise(A)
            assert len(dist_values(A)) == k
            Qs, st = extensions(A, 6, 7 - k)
            sets = []
            Cx = Classes()
            for Q in Qs:
                X = A + Q
                kX = len(dist_values(X)); assert kX <= 7
                Cx.add(X, kX)
            out['cores'].append({'A': to_str(A), **st, 'sets': [{'pts': to_str(r[0]), 'k': r[1], 'copies': r[2]} for r in Cx.reps]})
        out['time'] = round(time.time() - t0, 1)
        tmp = fo + f'.tmp{os.getpid()}'
        json.dump(out, open(tmp, 'w')); os.replace(tmp, fo)
        print(idx, 'k', k, 'reals', len(reals), 'sets', [len(c['sets']) for c in out['cores']],
              'cands', [c['ncand'] for c in out['cores']], round(out['time'], 1), flush=True)
    if not todo or ONCE: break
    time.sleep(60)
print('worker', w, 'done')
