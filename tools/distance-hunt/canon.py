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

if __name__ == '__main__':
    D = json.load(open(sys.argv[1]))
    out = {'cells': {}, 'items': []}
    labels = []; lab_idx = {}
    dup = 0
    for cell, c in D['cats'].items():
        n, d = map(int, cell.split(','))
        rig = [it for it in c['items'] if it['kind'] == 'rigid']
        fam = [it for it in c['items'] if it['kind'] != 'rigid']
        out['cells'][cell] = [len(rig), 1 if c.get('complete') else 0, len(fam)]
        seen = []
        for it in rig:
            z0 = [complex(float(x), float(y)) for x, y in it['coords']]
            z, res = refine.refine(it['coords'], it['cls'])
            mv = max(abs(a-b) for a,b in zip(z,z0))
            if mv > 1e-8 or res > 1e-10: print('REFINED', cell, it['label'][:50], '%.1e %.1e' % (mv, res))
            assert len(z) == n and classes(z) == d, (cell, it['label'], classes(z))
            cf = canon(z)
            for s in seen:
                if all(abs(a - b) < 1e-5 for a, b in zip(s, cf)):
                    dup += 1; print('DUP', cell, it['label'][:60]); break
            seen.append(cf)
            lab = it['label']
            if lab not in lab_idx: lab_idx[lab] = len(labels); labels.append(lab)
            out['items'].append([cell, lab_idx[lab], [v for q in cf for v in (round(q.real, 7), round(q.imag, 7))]])
    out['labels'] = labels
    print('items', len(out['items']), 'labels', len(labels), 'dups', dup)
    s = json.dumps(out, separators=(',', ':'))
    open(sys.argv[2], 'w').write(s)
    print('bytes', len(s))
