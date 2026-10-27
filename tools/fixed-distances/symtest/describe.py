import json, os, sys, itertools, collections, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify, dmat
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
def col_auts(n, cls, allow_colour_perm):
    E = list(itertools.combinations(range(n), 2)); idx = {e: k for k, e in enumerate(E)}
    auts = []
    for p in itertools.permutations(range(n)):
        m = {}; ok = True
        for k, (i, j) in enumerate(E):
            a, b = sorted((p[i], p[j])); c2 = cls[idx[(a, b)]]
            if allow_colour_perm:
                if m.setdefault(cls[k], c2) != c2: ok = False; break
            elif c2 != cls[k]: ok = False; break
        if ok and (not allow_colour_perm or len(set(m.values())) == len(m)): auts.append(p)
    return auts
for cell, idx in [(a.split(':')[0], int(a.split(':')[1])) for a in sys.argv[1:]]:
    it = D['cats'][cell]['items'][idx]; n = it['n']; cls = it['cls']; d = max(cls) + 1
    P = np.array([[float(x), float(y)] for x, y in it['coords']]); Dm = dmat(P)
    E = list(itertools.combinations(range(n), 2))
    vals = [np.mean([Dm[i, j] for k, (i, j) in enumerate(E) if cls[k] == c]) for c in range(d)]
    order = sorted(range(d), key=lambda c: vals[c]); rank = {c: r for r, c in enumerate(order)}
    m = vals[order[0]]
    print(f'=== ({cell}) #{idx}  {it["label"]}  catalog dim {it["dim"]}')
    print('  distance values (shortest = 1):', [round(vals[c] / m, 5) for c in order])
    print('  class sizes:', [sum(1 for k in range(len(E)) if cls[k] == c) for c in order])
    M = [['' if i == j else str(rank[cls[E.index((min(i, j), max(i, j)))]] + 1) for j in range(n)] for i in range(n)]
    print('  distance-class matrix (1 = shortest):'); [print('    ', ' '.join(f'{x:>1}' for x in row)) for row in M]
    print('  sample coords (shortest distance 1):', [[round(x, 5), round(y, 5)] for x, y in (P - P[0]) / m])
    print('  symmetry of sample:', classify(P)['type'], '; colouring automorphisms (colour-preserving):', len(col_auts(n, cls, False)), '; (up to colour permutation):', len(col_auts(n, cls, True)))
    if it.get('strat'): print('  strat:', it['strat'].replace('\n', ' | '))
