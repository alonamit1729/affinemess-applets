#!/usr/bin/env python3
"""Real realizability of the patterns of a level (after gen.py): msolve over Q with real root isolation on the
coordinate system (p0=(0,0), p1=(x1,0), colour of (0,1) = 1, Rabinowitsch variable t for distinct nonzero values),
exactly the pat_rigid.py / certify5d.py conventions.

Per pattern: 'empty' (no complex solution at all: the mod-p screen kept a false positive), 'rigid_real' (zero-
dimensional with real points: a real realization exists and there are finitely many), 'rigid_noreal'
(zero-dimensional, no real point: complex-only), 'posdim' (positive-dimensional; then --slices random hyperplane
sections are tried: 'posdim_real' if a section has a real point with distinct nonzero values, else 'posdim_unknown').

Usage: python3 realcheck.py n [--levels DIR] [--slices 4] [--precision 200]
Writes DIR/real{n}.json  [{"col":..., "k":..., "status":..., "nreal":...}] (checkpointed, resumable).
"""
import argparse, json, os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import realize
from realize import system_coords, msolve_real, real_status


def random_linear(XV):
    """Random integer linear form in the coordinate variables, msolve syntax (no zero coefficients: msolve
    silently misparses '0*x')."""
    while True:
        cs = [random.randint(-5, 5) for _ in XV]
        if any(cs):
            break
    terms = [f'{c}*{v}' for c, v in zip(cs, XV) if c]
    terms.append(str(random.randint(-3, 3)) or '0')
    return '+'.join(t for t in terms if t != '0').replace('+-', '-')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int)
    ap.add_argument('--levels', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'levels'))
    ap.add_argument('--file', default=None, help='pattern file to check (default DIR/L{n}.json)')
    ap.add_argument('--slices', type=int, default=4)
    ap.add_argument('--precision', type=int, default=200)
    ap.add_argument('--timeout', type=int, default=3600)
    args = ap.parse_args()
    n = args.n
    pats = [tuple(c) for c in json.load(open(args.file or os.path.join(args.levels, f'L{n}.json')))['patterns']]
    fn = os.path.join(args.levels, f'real{n}.json')
    done = json.load(open(fn)) if os.path.exists(fn) else []
    t0 = time.time()
    for i, col in enumerate(pats):
        if i < len(done):
            continue
        k = max(col) + 1
        V, P = system_coords(n, col)
        try:
            out = msolve_real(V, P, precision=args.precision, timeout=args.timeout)
            st = real_status(out)
        except Exception as ex:
            st = 'error'
        rec = {'col': list(col), 'k': k}
        if st == 'empty':
            rec['status'] = 'empty'
        elif st == 'real':
            rec['status'] = 'rigid_real'
        elif st == 'zerodim_noreal':
            rec['status'] = 'rigid_noreal'
        elif st == 'posdim':
            random.seed(i)
            XV = [v for v in V if v[0] in 'xy']
            found = False
            # expected dimension 2n+k-4-C(n,2) (free coordinates + free colour values - equations);
            # start there and add hyperplanes while the section stays positive-dimensional
            d0 = max(1, 2 * n + k - 4 - n * (n - 1) // 2)
            for d in range(d0, len(XV) + 1):       # cut with d random hyperplanes until zero-dimensional
                for _ in range(args.slices):
                    lins = [random_linear(XV) for _ in range(d)]
                    try:
                        s2 = real_status(msolve_real(V, P + lins, precision=args.precision, timeout=args.timeout))
                    except Exception:
                        s2 = 'error'
                    if s2 == 'real':
                        found = True
                        break
                    if s2 == 'posdim':
                        break          # need more hyperplanes
                if found or s2 != 'posdim':
                    break
            rec['status'] = 'posdim_real' if found else 'posdim_unknown'
        else:
            rec['status'] = st
        done.append(rec)
        json.dump(done, open(fn, 'w'))
        print(i, '/', len(pats), 'k', k, rec['status'], '%.0fs' % (time.time() - t0), flush=True)
    import collections
    print(collections.Counter(r['status'] for r in done))


if __name__ == '__main__':
    main()
