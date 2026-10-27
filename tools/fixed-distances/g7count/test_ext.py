import time, sys
from extend import *
from cores11 import all_sets
sets = all_sets()
i = int(sys.argv[1])
A, info = sets[i]
k = len(dist_values(A))
t0=time.time()
Qs, st = extensions(A, 5, 7-k)
print(info, 'k=',k, st, round(time.time()-t0,1))
C = Classes()
for Q in Qs:
    X = A+Q
    kX = len(dist_values(X))
    C.add(X, kX)
print('distinct 16-sets:', len(C.reps), [r[1] for r in C.reps])
