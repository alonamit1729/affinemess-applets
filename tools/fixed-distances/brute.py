import math, itertools, random
import numpy as np
from scipy.optimize import fsolve
random.seed(3); np.random.seed(3)
V=[0,1,2,3]; E=list(itertools.combinations(V,2))
# enumerate the 9 colourings
seen=set(); types=[]
for col in itertools.product(range(3),repeat=6):
    if len(set(col))<3: continue
    key=None
    for p in itertools.permutations(V):
        for q in itertools.permutations(range(3)):
            c=tuple(q[col[E.index(tuple(sorted((p[a],p[b]))))]] for a,b in E)
            if key is None or c<key: key=c
    if key not in seen: seen.add(key); types.append(key)
def pts(z): return [(0.0,0.0),(1.0,0.0),(z[0],z[1]),(z[2],z[3])]
def dm(P): return [[math.dist(P[i],P[j]) for j in range(4)] for i in range(4)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q); s=1.0
    for p in itertools.permutations(range(4)):
        s=D2[p[0]][p[1]]/D1[0][1]
        if all(abs(D1[i][j]*s-D2[p[i]][p[j]])<1e-5 for i in range(4) for j in range(4)): return True
    return False
total=0
for ti,col in enumerate(types):
    cls=[[E[i] for i in range(6) if col[i]==k] for k in range(3)]
    sols=[]
    for (p,q) in itertools.permutations(range(3),2):
        def F(z):
            P=pts(z); d=lambda e: math.dist(P[e[0]],P[e[1]])
            eqs=[]
            for c in cls:
                for e in c[1:]: eqs.append(d(e)**2-d(c[0])**2)
            eqs.append(d(cls[q][0])**2-4*d(cls[p][0])**2)
            return eqs
        for _ in range(400):
            z0=np.random.uniform(-3,3,4)
            z,info,ier,msg=fsolve(F,z0,full_output=True,xtol=1e-13)
            if ier!=1 or max(abs(v) for v in F(z))>1e-9: continue
            P=pts(z)
            ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(4),2))
            if ds[0]<1e-6: continue
            vals=[]
            for dd in ds:
                if not vals or dd-vals[-1]>1e-6: vals.append(dd)
            if len(vals)!=3: continue
            if any(congruent(P,Q) for Q in sols): continue
            sols.append(P)
    sizes=sorted(len(c) for c in cls)
    print('type',ti+1,sizes,cls,'->',len(sols))
    total+=len(sols)
print('total',total)
