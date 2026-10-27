"""Is every counterexample a sub-configuration (delete 1 or 2 points) of a symmetric catalog family?
Match by canonical edge colouring type (up to relabelling points and colours), then confirm numerically that the
colouring of the sub-configuration of the symmetric sample equals that of the counterexample (same type)."""
import json, os, sys, itertools, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import dmat, e_dim
rows = json.load(open(os.path.join(HERE, 'step1_families.json')))
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
def canon_col(n, cls):
    E = list(itertools.combinations(range(n), 2)); idx = {e: k for k, e in enumerate(E)}; best = None
    # cheap invariant first: sorted class sizes
    for p in itertools.permutations(range(n)):
        c = [None] * len(E)
        for k, (i, j) in enumerate(E):
            a, b = sorted((p[i], p[j])); c[idx[(a, b)]] = cls[k]
        rl = {}; out = []
        for x in c:
            if x not in rl: rl[x] = len(rl)
            out.append(rl[x])
        t = tuple(out)
        if best is None or t < best: best = t
    return best
def invariant(n, cls):
    E = list(itertools.combinations(range(n), 2))
    size = collections.Counter(cls)
    vd = []
    for v in range(n):
        vd.append(tuple(sorted(size[cls[k]] for k, (i, j) in enumerate(E) if v in (i, j))))
    return (tuple(sorted(size.values())), tuple(sorted(vd)))
def colouring(P):
    n = len(P); Dm = dmat(P); E = list(itertools.combinations(range(n), 2)); sc = Dm.max()
    vals = []
    cls = []
    for i, j in E:
        v = Dm[i, j]; k = next((t for t, u in enumerate(vals) if abs(u - v) < 1e-6 * sc), None)
        if k is None: vals.append(v); k = len(vals) - 1
        cls.append(k)
    return cls, len(vals)
if __name__ == '__main__':
  targets = ['5,4:473', '5,4:474', '5,4:475', '5,4:477', '5,4:491', '5,4:492', '6,4:23', '6,6:3', '6,6:4', '6,7:0', '6,7:1', '6,7:2', '6,7:3', '6,7:8', '7,7:4']
  tcanon = {}
  for a in targets:
      cell, idx = a.split(':'); it = D['cats'][cell]['items'][int(idx)]
      tcanon[a] = (it['n'], canon_col(it['n'], it['cls']), invariant(it['n'], it['cls']), max(it['cls']) + 1)
  # sources: every family or rigid item in the catalog with nontrivial symmetry and n <= 8 (deletions of 1 or 2 points)
  found = collections.defaultdict(list)
  for key, c in D['cats'].items():
      n, d = map(int, key.split(','))
      if n > 9 or n < 6: continue
      for j, it in enumerate(c['items']):
          P = np.array([[float(x), float(y)] for x, y in it['coords']])
          from symmetry import classify
          sym = classify(P)['type']
          if sym == 'trivial' or it['kind'] != 'family': continue
          for k in (1, 2, 3):
              if n - k < 5: continue
              for dele in itertools.combinations(range(n), k):
                  keep = [i for i in range(n) if i not in dele]; Q = P[keep]
                  cls, dq = colouring(Q)
                  inv = None
                  for a, (tn, tc, tinv, td) in tcanon.items():
                      if tn != len(keep) or td != dq: continue
                      if inv is None: inv = invariant(len(keep), cls)
                      if inv != tinv: continue
                      if canon_col(len(keep), cls) == tc:
                          found[a].append((key, j, it['kind'], it['dim'], sym, dele, dq, it['label'][:60]))
  for a in targets:
      print(a, 'e =', e_dim(*map(int, a.split(':')[0].split(','))), ':')
      seen = set()
      for f in found[a]:
          tag = (f[0], f[1], len(f[5]))
          if tag in seen: continue
          seen.add(tag)
          print(f'    subset of ({f[0]}) #{f[1]} [{f[2]}, dim {f[3]}, symmetry {f[4]}] deleting {len(f[5])} point(s) {f[5]} -> {f[6]} distances; label: {f[7]}')
      if not found[a]: print('    not found as a 1-, 2- or 3-point deletion of a symmetric catalog configuration with n<=9')
