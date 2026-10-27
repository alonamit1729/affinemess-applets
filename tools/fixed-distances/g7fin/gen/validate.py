#!/usr/bin/env python3
"""Validation of the generator output against Szollosi-Ostergard (Table 9) and the certified small cases.

Usage: python3 validate.py [--levels DIR]
Checks, for every level present in DIR:
  * raw / pruned counts against Table 9 of arXiv:1804.06040;
  * n=4: pruned list == real-realizable K4 types of scratchpad/cert4{2..6}.json;
  * n=5: pruned list contains every real-realizable K5 type of scratchpad/cert5{2..6}.json (and lists the
         complex-only ones);
  * every n: every n-point sub-pattern of the 225 known realizable 10-point patterns (g7fin/patterns_10.json)
    is in the raw list and in the pruned list;
  * n=10: comparison with patterns_10.json (which known patterns are found, which candidates are new).
"""
import argparse, collections, itertools, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import patlib
from patlib import canon, induced, dict_to_col

SCRATCH = os.path.abspath(os.path.join(HERE, '..', '..'))
SO = {4: (25, 24), 5: (1193, 922), 6: (537230, 20229), 7: (14732, 4667), 8: (1565, 1565),
      9: (635, None), 10: (228, None), 11: (62, None), 12: (15, None), 13: (3, None), 14: (0, None)}


def load(fn):
    return [tuple(c) for c in json.load(open(fn))['patterns']] if os.path.exists(fn) else None


def cert_real(n, want='real'):
    """Canonical forms of the certified types: want='real' (a rigid realization, a prime with real points or
    real slice hits), 'unresolved' (decomposition timed out: status unknown), 'all'."""
    out = set()
    for k in range(2, 7):
        fn = os.path.join(SCRATCH, f'cert{n}{k}.json')
        if not os.path.exists(fn):
            continue
        for r in json.load(open(fn)):
            isreal = bool(r.get('rigid')) or any(p.get('real') or p.get('real_slices') for p in r['primes'])
            unres = bool(r.get('unresolved'))
            if want == 'all' or (want == 'real' and isreal) or (want == 'unresolved' and unres and not isreal):
                out.add(canon(n, tuple(r['col'])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--levels', default=os.path.join(HERE, 'levels'))
    args = ap.parse_args()
    known10 = [dict_to_col(10, p['col']) for p in json.load(open(os.path.join(SCRATCH, 'g7fin', 'patterns_10.json')))]
    print(f'{len(known10)} known realizable 10-point patterns loaded')
    for n in range(4, 15):
        raw = load(os.path.join(args.levels, f'raw{n}.json'))
        L = load(os.path.join(args.levels, f'L{n}.json'))
        if raw is None and L is None:
            continue
        so_raw, so_pr = SO.get(n, (None, None))
        print(f'\n=== level {n} ===')
        if raw is not None:
            print(f'raw candidates: {len(raw)}   (S-O row 1: {so_raw})  {"OK" if so_raw == len(raw) else "MISMATCH"}')
        if L is not None:
            kc = collections.Counter(max(c) + 1 for c in L)
            print(f'pruned (complex-realizable): {len(L)}   (S-O row 2: {so_pr})  '
                  f'{"OK" if so_pr == len(L) else ("n/a" if so_pr is None else "MISMATCH")}   by #colours {dict(sorted(kc.items()))}')
        if n <= 10:
            subs = set()
            for col in known10:
                for U in itertools.combinations(range(10), n):
                    subs.add(canon(n, induced(10, col, list(U))))
            if raw is not None:
                miss = [c for c in subs if c not in set(raw)]
                print(f'{len(subs)} distinct {n}-subpatterns of the known 10-point patterns; missing from raw: {len(miss)}')
            if L is not None:
                miss = [c for c in subs if c not in set(L)]
                print(f'   missing from pruned: {len(miss)}')
        if n in (4, 5) and L is not None:
            real = cert_real(n)
            unres = cert_real(n, 'unresolved')
            Ls = set(L)
            print(f'certified real-realizable {n}-types: {len(real)} (+{len(unres)} unresolved); in pruned list: {len(real & Ls)}; '
                  f'real types NOT in pruned list: {len(real - Ls)}; unresolved in pruned list: {len(unres & Ls)}; '
                  f'pruned but neither (complex-only): {len(Ls - real - unres)} '
                  f'by #colours {dict(sorted(collections.Counter(max(c)+1 for c in Ls-real-unres).items()))}')
        rfn = os.path.join(args.levels, f'real{n}.json')
        if os.path.exists(rfn):
            rc = json.load(open(rfn))
            st = collections.Counter(r['status'] for r in rc)
            print(f'realcheck ({len(rc)} patterns): {dict(st)}')
            if n in (4, 5):
                real = cert_real(n)
                unres = cert_real(n, 'unresolved')
                rreal = {canon(n, tuple(r['col'])) for r in rc if r['status'] in ('rigid_real', 'posdim_real')}
                rnot = {canon(n, tuple(r['col'])) for r in rc if r['status'] not in ('rigid_real', 'posdim_real')}
                checked = rreal | rnot
                print(f'   real by realcheck: {len(rreal)}; of these certified real: {len(rreal & real)}, '
                      f'certified-unresolved: {len(rreal & unres)}, certified NOT real (conflict): {len(rreal - real - unres)}; '
                      f'certified real but not found real by realcheck: {len(real & rnot)} '
                      f'{dict(collections.Counter(r["status"] for r in rc if canon(n, tuple(r["col"])) in real & rnot))}')
        if n == 10 and L is not None:
            Ls = set(L)
            found = [c for c in known10 if canon(10, c) in Ls]
            new = [c for c in L if c not in {canon(10, k) for k in known10}]
            print(f'known patterns found in L10: {len(found)}/{len(known10)}; candidates not among the known: {len(new)}')
            for c in new:
                print('   new candidate:', c)


if __name__ == '__main__':
    main()
