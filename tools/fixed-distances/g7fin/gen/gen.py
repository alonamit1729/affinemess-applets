#!/usr/bin/env python3
"""Hereditary (orderly) generation of candidate colouring patterns of n-point planar sets with <= 6 distances.

Level n is built from the pruned list L_{n-1}: every parent pattern is extended by one vertex in all ways
(colours of the n-1 new edges, at most 6 colours in total, new colours introduced in increasing order);
an extension is kept only if every induced sub-pattern that contains the new vertex is in the pruned list
of its order (checked incrementally during the extension: 4- and 5-subsets by table lookup, larger ones by
nauty canonical form).  Survivors are canonicalised and deduplicated -> raw list (Szollosi-Ostergard
Table 9, row 1).  Then each raw pattern is tested for complex realizability (Groebner basis of the
rank<=2 / distance system with the Rabinowitsch inequation, mod a prime) -> pruned list L_n (row 2).

Files (in --out, default ./levels):
  raw{n}.json   {"n","count","patterns":[rgf tuples]}      canonical raw candidates
  test{n}.json  {"method","prime","status":[1/0/-1 per raw pattern]}  checkpoint of the realizability tests
  L{n}.json     {"n","count","patterns":[...]}            pruned list
  summary.json  counts per level
Resuming: rerun the same command; finished levels are skipped, a half-tested level continues.

Usage: python3 gen.py --to 10 [--out DIR] [--workers 1] [--method auto|gram|coords] [--no-confirm] [--gen-only]
"""
import argparse, itertools, json, os, sys, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import patlib
from patlib import EDGES, EIDX, rgf, canon, all_rgf

MAXCOL = 6


# ------------------------------------------------------------------ level data
def load_level(outdir, n):
    fn = os.path.join(outdir, f'L{n}.json')
    if not os.path.exists(fn):
        return None
    return [tuple(c) for c in json.load(open(fn))['patterns']]


def save_json(fn, obj):
    tmp = fn + f'.tmp{os.getpid()}'
    with open(tmp, 'w') as f:
        json.dump(obj, f)
    os.replace(tmp, fn)


def lookup_table(n, Lset):
    """Set of all rgf-labelled colourings of K_n whose canonical form is in Lset (n = 4, 5)."""
    return {col for col in all_rgf(len(EDGES[n])) if canon(n, col) in Lset}


# ------------------------------------------------------------------ extension of one parent
class Extender:
    def __init__(self, n, tables, Lsets):
        """tables: {4: set, 5: set} of labelled rgf colourings; Lsets: {j: set of canonical} for j>=6."""
        self.n = n
        self.tables = tables
        self.Lsets = Lsets
        v = n - 1
        self.v = v
        # for each step i, the list of (j, seq) sub-pattern extraction recipes for subsets containing i and v
        self.checks = [[] for _ in range(v)]
        eidx = EIDX[n - 1]
        for i in range(v):
            for j in range(4, min(i + 2, n - 1) + 1):
                if j == 4 and i >= 3 and n - 1 >= 5:
                    continue   # every K4 containing i and v lies in a K5 {a,b,c,i,v} checked at this step
                for T in itertools.combinations(range(i), j - 2):
                    U = sorted(T + (i, v))
                    seq = []
                    for a, b in itertools.combinations(U, 2):
                        seq.append(-1 - a if b == v else eidx[(a, b)])
                    self.checks[i].append((j, tuple(seq)))
        # assembling the full colouring of K_n in EDGES[n] order
        self.full = []
        for a, b in EDGES[n]:
            self.full.append(-1 - a if b == v else eidx[(a, b)])

    def ok(self, j, seq, P, c):
        col = rgf([P[x] if x >= 0 else c[-1 - x] for x in seq])
        if j <= 5:
            return col in self.tables[j]
        return canon(j, col) in self.Lsets[j]

    def extend(self, P):
        """Yield all labelled extensions (full rgf colourings of K_n) of parent P passing hereditary checks."""
        n, v = self.n, self.v
        k = max(P) + 1
        c = [0] * v
        checks, tables, Lsets, full = self.checks, self.tables, self.Lsets, self.full
        out = []

        def rec(i, kk):
            if i == v:
                out.append(rgf([P[x] if x >= 0 else c[-1 - x] for x in full]))
                return
            for col in range(min(MAXCOL, kk + 1)):
                c[i] = col
                good = True
                for j, seq in checks[i]:
                    sub = rgf([P[x] if x >= 0 else c[-1 - x] for x in seq])
                    if j <= 5:
                        if sub not in tables[j]:
                            good = False
                            break
                    elif canon(j, sub) not in Lsets[j]:
                        good = False
                        break
                if good:
                    rec(i + 1, max(kk, col + 1))
        rec(0, k)
        return out


_EXT = None


def _init_worker(n, tables, Lsets):
    global _EXT
    _EXT = Extender(n, tables, Lsets)


def _work_parents(parents):
    res = {}
    for P in parents:
        for col in _EXT.extend(P):
            cf = canon(_EXT.n, col)
            res[cf] = 1
    return list(res)


def generate_raw(n, Lprev, tables, Lsets, workers, log, part_file=None):
    """part_file: checkpoint of the generation ({'done': chunks processed, 'raw': [...]}), written every 50 chunks
    (atomically); chunks are processed in order so a resumed run skips the first `done` chunks."""
    parents = Lprev
    t0 = time.time()
    raw = set()
    chunks = [parents[i:i + 8] for i in range(0, len(parents), 8)]
    done = 0
    if part_file and os.path.exists(part_file):
        with open(part_file) as f:
            part = json.load(f)
        done = part['done']; raw = set(tuple(c) for c in part['raw'])
        log(f'  level {n}: resuming raw generation at chunk {done}/{len(chunks)} with {len(raw)} raw')

    def ckpt():
        if part_file:
            tmp = part_file + '.tmp'
            with open(tmp, 'w') as f:
                json.dump({'done': done, 'raw': [list(c) for c in raw]}, f)
            os.replace(tmp, part_file)
    todo = chunks[done:]
    if workers > 1:
        with Pool(workers, initializer=_init_worker, initargs=(n, tables, Lsets)) as pool:
            for res in pool.imap(_work_parents, todo):
                raw.update(res)
                done += 1
                if done % 50 == 0:
                    log(f'  level {n}: {done}/{len(chunks)} parent chunks, {len(raw)} raw so far, {time.time()-t0:.0f}s'); ckpt()
    else:
        _init_worker(n, tables, Lsets)
        for ch in todo:
            raw.update(_work_parents(ch))
            done += 1
            if done % 50 == 0:
                log(f'  level {n}: {done}/{len(chunks)} parent chunks, {len(raw)} raw so far, {time.time()-t0:.0f}s'); ckpt()
    log(f'  level {n}: raw generation done, {len(raw)} candidates, {time.time()-t0:.0f}s')
    if part_file and os.path.exists(part_file):
        os.remove(part_file)
    return sorted(raw)


# ------------------------------------------------------------------ realizability tests
_TEST_CFG = {}


def _init_test(cfg):
    global _TEST_CFG
    _TEST_CFG = cfg
    import realize  # noqa


def _test_one(item):
    import realize
    idx, n, col = item
    cfg = _TEST_CFG
    method = cfg['method']
    if method == 'auto':
        method = 'gram' if n <= 7 else 'coords'
    try:
        ok = realize.complex_realizable(n, col, method, prime=realize.PRIMES[0])
        if not ok and cfg['confirm']:
            ok2 = realize.complex_realizable(n, col, method, prime=realize.PRIMES[1])
            if ok2:
                return idx, 2   # disagreement between primes: keep, flag
        return idx, int(ok)
    except Exception as ex:  # timeout etc.: keep as candidate, flag
        return idx, 3


def test_level(n, raw, outdir, cfg, workers, log):
    fn = os.path.join(outdir, f'test{n}.json')
    if os.path.exists(fn):
        st = json.load(open(fn))
        status = st['status']
        if len(status) != len(raw) or st.get('method') != cfg['method']:
            log(f'  level {n}: stale test file, restarting tests')
            status = [-1] * len(raw)
    else:
        status = [-1] * len(raw)
    todo = [(i, n, raw[i]) for i in range(len(raw)) if status[i] == -1]
    log(f'  level {n}: testing {len(todo)} of {len(raw)} candidates ({cfg["method"]}, {workers} workers)')
    t0 = time.time()
    ndone = 0

    def save():
        save_json(fn, {'n': n, 'method': cfg['method'], 'prime': cfg['prime'], 'confirm': cfg['confirm'], 'status': status})
    if workers > 1:
        with Pool(workers, initializer=_init_test, initargs=(cfg,)) as pool:
            for idx, s in pool.imap_unordered(_test_one, todo, chunksize=4):
                status[idx] = s
                ndone += 1
                if ndone % 500 == 0:
                    save()
                    el = time.time() - t0
                    log(f'  level {n}: tested {ndone}/{len(todo)}, {el:.0f}s, {el/ndone*1000:.1f} ms each, '
                        f'eta {el/ndone*(len(todo)-ndone)/60:.1f} min, kept so far {sum(1 for s in status if s>0)}')
    else:
        _init_test(cfg)
        for item in todo:
            idx, s = _test_one(item)
            status[idx] = s
            ndone += 1
            if ndone % 500 == 0:
                save()
                el = time.time() - t0
                log(f'  level {n}: tested {ndone}/{len(todo)}, {el:.0f}s, {el/ndone*1000:.1f} ms each')
    save()
    log(f'  level {n}: tests done in {time.time()-t0:.0f}s; '
        f'kept {sum(1 for s in status if s>0)} (flagged prime-disagreement {status.count(2)}, errors {status.count(3)})')
    return status


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--to', type=int, default=10)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'levels'))
    ap.add_argument('--workers', type=int, default=1, help='processes (this container has one effective core)')
    ap.add_argument('--method', default='auto', choices=['auto', 'gram', 'coords'],
                    help='auto = gram (Szollosi-Ostergard minors) for n<=7, coords for n>=8 (fewer, cheaper polynomials)')
    ap.add_argument('--no-confirm', action='store_true',
                    help='do not re-test discards with a second prime (default: every discard is confirmed mod a 2nd prime)')
    ap.add_argument('--no-prune', action='store_true', help='skip realizability tests (L_n := raw_n)')
    ap.add_argument('--gen-only', action='store_true', help='stop after generating the raw list of the last level')
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    logf = open(args.log or os.path.join(args.out, 'gen.log'), 'a')

    def log(msg):
        line = time.strftime('%H:%M:%S ') + msg
        print(line, flush=True)
        logf.write(line + '\n')
        logf.flush()
    import realize
    cfg = {'method': args.method, 'confirm': not args.no_confirm, 'prime': realize.PRIMES[0]}
    summary_fn = os.path.join(args.out, 'summary.json')
    summary = json.load(open(summary_fn)) if os.path.exists(summary_fn) else {}

    # level 4: all colourings of K4 with <= 6 colours
    L = {}
    for n in range(4, args.to + 1):
        Ln = load_level(args.out, n)
        if Ln is not None:
            L[n] = Ln
            log(f'level {n}: loaded {len(Ln)} pruned patterns')
            continue
        log(f'level {n}: start')
        rawfn = os.path.join(args.out, f'raw{n}.json')
        if os.path.exists(rawfn):
            raw = [tuple(c) for c in json.load(open(rawfn))['patterns']]
            log(f'  level {n}: loaded {len(raw)} raw candidates')
        else:
            if n == 4:
                raw = sorted({canon(4, col) for col in all_rgf(6)})
            else:
                tables = {}
                Lsets = {}
                for j in range(4, n):
                    if j <= 5:
                        tables[j] = lookup_table(j, set(L[j]))
                    else:
                        Lsets[j] = set(L[j])
                raw = generate_raw(n, L[n - 1], tables, Lsets, args.workers, log, part_file=os.path.join(args.out, f'rawpart{n}.json'))
            save_json(rawfn, {'n': n, 'count': len(raw), 'patterns': [list(c) for c in raw]})
        if args.gen_only and n == args.to:
            log(f'level {n}: raw {len(raw)} (gen-only, stopping)')
            break
        if args.no_prune:
            status = [1] * len(raw)
        else:
            status = test_level(n, raw, args.out, cfg, args.workers, log)
        Ln = [raw[i] for i in range(len(raw)) if status[i] > 0]
        save_json(os.path.join(args.out, f'L{n}.json'), {'n': n, 'count': len(Ln), 'patterns': [list(c) for c in Ln]})
        L[n] = Ln
        kc = {}
        for c in Ln:
            kc[max(c) + 1] = kc.get(max(c) + 1, 0) + 1
        summary[str(n)] = {'raw': len(raw), 'pruned': len(Ln), 'by_colours': kc,
                           'flagged': status.count(2), 'errors': status.count(3)}
        save_json(summary_fn, summary)
        log(f'level {n}: raw {len(raw)}  pruned {len(Ln)}  by #colours {kc}')
    log('done')


if __name__ == '__main__':
    main()
