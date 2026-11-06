"""Canonical form under similarity (reflections included), shared with the JS game."""
import json, sys, collections, refine
from functools import cmp_to_key
TOL = 1e-7

def cmpc(a, b):
    if abs(a.real - b.real) > TOL: return -1 if a.real < b.real else 1
    if abs(a.imag - b.imag) > TOL: return -1 if a.imag < b.imag else 1
    return 0

def cmpl(A, B):
    for a, b in zip(A, B):
        c = cmpc(a, b)
        if c: return c
    return 0

def canon(z):
    n = len(z)
    ds = [abs(z[i] - z[j]) for i in range(n) for j in range(i + 1, n)]
    dmin = min(ds)
    best = None
    for i in range(n):
        for j in range(n):
            if i == j or abs(z[i] - z[j]) > dmin * (1 + 1e-7): continue
            for refl in (False, True):
                w = [(q - z[i]) / (z[j] - z[i]) for q in z]
                if refl: w = [q.conjugate() for q in w]
                w.sort(key=cmp_to_key(cmpc))
                if best is None or cmpl(w, best) < 0: best = w
    return best

def classes(z):
    ds = sorted(abs(z[i] - z[j]) for i in range(len(z)) for j in range(i + 1, len(z)))
    cl = [ds[0]]
    for d in ds[1:]:
        if d - cl[-1] > 1e-7 * d: cl.append(d)
    return len(cl)

def fingerprint(job):
    cell, num, lab, coords, cls, n, d = job
    z, res = refine.refine(coords, cls)
    if len(z) != n or classes(z) != d: return (cell, num, lab, None, res)
    cf = canon(z)
    return (cell, num, lab, [v for q in cf for v in (round(q.real, 7), round(q.imag, 7))], res)


BIG = 1000   # cells with more rigid sets than this are written to catalogs/<n>_<d>.json and loaded on demand

if __name__ == '__main__':
    # python3 canon.py <fixed-distances applet_data.json> catalog.json [old catalog.json for verdicts]
    import os
    from multiprocessing import Pool
    here = os.path.dirname(os.path.abspath(__file__))
    D = json.load(open(sys.argv[1]))
    old = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else None
    oldv = collections.defaultdict(list)          # cell -> [(flat, verdict)]
    if old:
        for it in old['items']:
            if len(it) > 4: oldv[it[0]].append((it[2], it[4]))
    elif os.path.exists(os.path.join(here, 'verdicts.json')):
        print('no old catalog given: verdicts are matched by shape only when one is given')
    jobs = []
    out = {'cells': {}, 'items': []}
    for cell, c in D['cats'].items():
        n, d = map(int, cell.split(','))
        rig = [(k + 1, it) for k, it in enumerate(c['items']) if it['kind'] == 'rigid']
        fam = [it for it in c['items'] if it['kind'] != 'rigid']
        out['cells'][cell] = [len(rig), 1 if c.get('complete') else 0, len(fam)] + ([1] if len(rig) > BIG else [])
        jobs += [(cell, num, it['label'], it['coords'], it['cls'], n, d) for num, it in rig]
    with Pool(2) as P: res = P.map(fingerprint, jobs, chunksize=64)
    per = collections.defaultdict(list); bad = 0; dup = 0; ver = []
    for cell, num, lab, flat, r in res:
        if flat is None: bad += 1; print('SKIP (distance count off after refining)', cell, num, lab[:50]); continue
        v = None
        for f2, v2 in oldv.get(cell, []):
            if all(abs(a - b) < 1e-5 for a, b in zip(f2, flat)): v = v2; break
        per[cell].append((num, lab, flat, v)); ver.append([cell, num, v])
    for cell, L in per.items():
        seen = set()
        for num, lab, flat, v in L:
            key = tuple(round(x, 4) for x in flat)
            if key in seen: dup += 1
            seen.add(key)
    os.makedirs(os.path.join(here, 'catalogs'), exist_ok=True)
    labels = []; lab_idx = {}
    for cell, L in sorted(per.items(), key=lambda a: tuple(map(int, a[0].split(',')))):
        if out['cells'][cell][0] > BIG:
            bl = []; bi = {}; items = []
            for num, lab, flat, v in L:
                if lab not in bi: bi[lab] = len(bl); bl.append(lab)
                items.append([bi[lab], [round(x, 6) for x in flat], num, v])
            fn = os.path.join(here, 'catalogs', cell.replace(',', '_') + '.json')
            open(fn, 'w').write(json.dumps({'cell': cell, 'labels': bl, 'items': items}, separators=(',', ':')))
            print('wrote', fn, os.path.getsize(fn))
            continue
        for num, lab, flat, v in L:
            if lab not in lab_idx: lab_idx[lab] = len(labels); labels.append(lab)
            out['items'].append([cell, lab_idx[lab], flat, num, v])
    out['labels'] = labels
    json.dump(ver, open(os.path.join(here, 'verdicts.json'), 'w'), separators=(',', ':'))
    known = sum(1 for _, _, v in ver if v)
    print('sets', len(res), 'embedded', len(out['items']), 'skipped', bad, 'duplicates', dup, 'with verdict', known)
    s = json.dumps(out, separators=(',', ':'))
    open(sys.argv[2], 'w').write(s)
    print('bytes', len(s))
