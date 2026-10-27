"""FAST PASS: all 16-point 7-distance supersets of the 68 eleven-point <=6-distance sets."""
import time, json, os
from extend import *
from cores11 import all_sets
sets = all_sets()
C = Classes()
log = []
t00 = time.time()
for i, (A, info) in enumerate(sets):
    k = len(dist_values(A))
    t0 = time.time()
    Qs, st = extensions(A, 5, 7 - k)
    found = []
    for Q in Qs:
        X = A + Q
        kX = len(dist_values(X))
        assert kX <= 7
        if kX < 7: print('  !! set with', kX, 'distances', info)
        idx, new = C.add(X, {'first_core': info})
        found.append(idx)
    rec = {'core': info, 'k': k, **st, 'time': round(time.time() - t0, 1), 'classes': sorted(set(found))}
    log.append(rec); print(i, rec, flush=True)
print('total time', round(time.time() - t00, 1))
print('distinct 16-point 7-distance sets:', len(C.reps))
out = [{'pts': to_str(r[0]), 'pts_float': to_floats(r[0]), 'info': r[1], 'hits': r[2], 'pattern': pattern(r[0])} for r in C.reps]
json.dump({'log': log, 'sets': out}, open('fast11_result.json', 'w'), indent=1)
