"""Certify the Jacobian rank at a real point of each "stubborn" five-point type (stubborn.py), so that the
certificate of those types no longer has a numerical step.

For the 21 colouring types of K5 whose primary decomposition never finished, stubborn.py showed that the saturated
ideal (Rabinowitsch variable t enforcing "values distinct and nonzero") has dimension exactly e, found real points
on random slices, and reported the Jacobian rank there numerically. Here, for each type:

  1. the first slice of stubborn.py with real solutions is solved again by msolve (same random seed, so the same
     slice), which returns isolating boxes of width 2^-200, each containing exactly one real solution p;
  2. at the first such box, the 10 x |V| Jacobian of the ten distance equations in V = (x1, x2, y2, x3, y3, x4, y4,
     s1, ..., s_{d-1}) is evaluated in interval arithmetic (mpmath.iv, 300 bits), and a 10 x 10 minor is shown to
     exclude 0 on the whole box (interval Gaussian elimination with columns chosen by a pivoted QR of the
     midpoint matrix).

Then rank J(p) = 10 = |V| - e. The ten equations cut out a real-analytic manifold of dimension e near p, which
contains the variety; since every component of the variety through p has complex dimension e (step of stubborn.py),
the variety coincides with that manifold near p. So the real points of the type form an e-dimensional family near
p: the type is a genuine e-parameter family, as counted. p has t != 0, so it is non-degenerate.

    cd tools/fixed-distances && python3 certify_stubborn.py      (msolve; about a minute)
Writes data/certify_stubborn.json.
"""
import sys, os, json, glob, re, random, subprocess
from fractions import Fraction
import numpy as np, scipy.linalg
from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
os.makedirs('tmp', exist_ok=True)
iv.prec = 300


def prelude(d):
    """certify5d.py's definitions (as in stubborn.py) for d colours"""
    src = open('certify5d.py').read().split('import os')[0]
    argv = sys.argv; sys.argv = [argv[0], str(d)]
    G = {'__name__': 'certify5d'}
    exec(src, G)
    sys.argv = argv
    return G


def parse_boxes(out):
    """msolve output -> list of real solutions, each a list of (lo, hi) Fractions"""
    s = out.replace('\n', '').replace(' ', '')
    if s.startswith('[-1]') or s.startswith('[1,'): return []
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s.rstrip(':'))
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
    return tree[1][1] if isinstance(tree[1], list) and len(tree[1]) > 1 else []


def interval_det_nonzero(M):
    """interval Gaussian elimination on a square list-of-lists of iv.mpf; True if det provably != 0"""
    A = [row[:] for row in M]; k = len(A)
    for c in range(k):
        piv = max(range(c, k), key=lambda r: float(abs(A[r][c]).a))   # largest lower bound
        if 0 in A[piv][c]: return False
        A[c], A[piv] = A[piv], A[c]
        for r in range(c + 1, k):
            f = A[r][c] / A[c][c]
            A[r] = [A[r][j] - f * A[c][j] for j in range(k)]
    return True


def certify(path):
    sb = json.load(open(path)); idx = sb['idx']; d = int(os.path.basename(path).split('_')[1])
    G = prelude(d)
    col = G['classes'][idx]; assert list(col) == sb['col']
    V, XV, expand, system, ineq = G['V'], G['XV'], G['expand'], G['system'], G['ineq']
    e = 2 * 5 - 4 - 10 + d
    eqs = system(col) + [expand(ineq)]; VT = V + ['t']
    random.seed(idx)
    for k in range(6):   # the slices of stubborn.py, in the same order
        lins = ['+'.join(f'({random.randint(-5, 5)})*{v}' for v in XV) + f'-({random.randint(-3, 3)})' for _ in range(e)]
        boxes = parse_boxes(G['run_msolve'](eqs + [expand(l) for l in lins], VT))
        if boxes: break
    else:
        return {'idx': idx, 'd': d, 'certified': False, 'reason': 'no real slice point'}
    import sympy as sp
    syms = [sp.Symbol(v) for v in V]
    F = [sp.sympify(q.replace('^', '**'), locals={v: s for v, s in zip(V, syms)}) for q in system(col)]
    Jsym = sp.Matrix(F).jacobian(syms)
    Jf = [[sp.lambdify(syms, Jsym[r, c], modules=[{}]) for c in range(len(V))] for r in range(10)]   # plain arithmetic
    for b, box in enumerate(boxes):
        q = lambda r: iv.mpf(r.numerator) / r.denominator          # an interval enclosing the rational r
        X = [iv.mpf([q(lo).a, q(hi).b]) for lo, hi in box[:len(V)]]
        J = [[iv.mpf(f(*X)) for f in row] for row in Jf]
        mid = np.array([[float(x.mid) for x in row] for row in J])
        _, _, perm = scipy.linalg.qr(mid, pivoting=True)
        cols = sorted(perm[:10])
        if interval_det_nonzero([[J[r][c] for c in cols] for r in range(10)]):
            width = max(float(hi - lo) for lo, hi in box)
            return {'idx': idx, 'd': d, 'e': e, 'certified': True, 'slice': k, 'box': b, 'boxes_on_slice': len(boxes),
                    'box_width': width, 'minor_columns': [V[c] for c in cols],
                    'point': [float((lo + hi) / 2) for lo, hi in box[:len(V)]]}
    return {'idx': idx, 'd': d, 'certified': False, 'reason': 'no minor certified'}


def main():
    res = [certify(p) for p in sorted(glob.glob('data/stubborn_*.json'))]
    for r in res: print(r['d'], r['idx'], r['certified'], r.get('minor_columns', r.get('reason')), flush=True)
    ok = all(r['certified'] for r in res)
    json.dump({'all_certified': ok, 'types': res}, open('data/certify_stubborn.json', 'w'), indent=1)
    print(len(res), 'types;', 'all certified' if ok else 'NOT all certified')


if __name__ == '__main__':
    main()
