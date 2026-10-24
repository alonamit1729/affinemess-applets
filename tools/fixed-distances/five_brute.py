import math, itertools, json, sys
import numpy as np
from scipy.optimize import least_squares
np.random.seed(7)
E5=list(itertools.combinations(range(5),2)); E4=list(itertools.combinations(range(4),2))
# all colourings
cols=np.array(list(itertools.product(range(3),repeat=10)),dtype=np.int64)
cols=cols[(cols==0).any(1)&(cols==1).any(1)&(cols==2).any(1)]
pw=3**np.arange(10)[::-1]
best=None
for p in itertools.permutations(range(5)):
    eperm=[E5.index(tuple(sorted((p[a],p[b])))) for a,b in E5]
    permuted=cols[:,eperm]
    for q in itertools.permutations(range(3)):
        qa=np.array(q)
        code=(qa[permuted]*pw).sum(1)
        best=code if best is None else np.minimum(best,code)
canon=np.unique(best)
print('classes',len(canon))
# realizable K4 colouring types (canonical codes of 4-vertex colourings, using colours up to relabel)
def canon4(col):
    key=None
    for p in itertools.permutations(range(4)):
        for q in itertools.permutations(range(3)):
            c=tuple(q[col[E4.index(tuple(sorted((p[a],p[b]))))]] for a,b in E4)
            if key is None or c<key: key=c
    return key
ok4=set()
# 3-colour types: all 9 are realizable
for col in itertools.product(range(3),repeat=6):
    if len(set(col))==3: ok4.add(canon4(col))
# 2-colour realizable types: square C4+matching, rhombus 5+1, tri+centre K3+K13, trapezoid P4+P4, kite (K3+pendant)+P3
def col_from(edges1):
    return tuple(0 if e in edges1 else 1 for e in E4)
two_types=[ [(0,1),(1,2),(2,3),(0,3)], [(0,1),(0,2),(0,3),(1,2),(1,3)], [(0,1),(0,2),(1,2)], [(0,1),(1,2),(2,3)], [(0,1),(0,2),(1,2),(0,3)] ]
for t in two_types: ok4.add(canon4(col_from(t)))
def decode(code):
    return [(code//3**(9-i))%3 for i in range(10)]
survivors=[]
for code in canon:
    col=decode(int(code))
    good=True
    for S in itertools.combinations(range(5),4):
        sub=tuple(col[E5.index((S[a],S[b]))] for a,b in E4)
        if canon4(sub) not in ok4: good=False;break
    if good: survivors.append(col)
print('survivors',len(survivors))
def pdist(P): return [math.dist(P[i],P[j]) for i,j in E5]
def dm(P): return [[math.dist(P[i],P[j]) for j in range(5)] for i in range(5)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(5)):
        s=D2[p[0]][p[1]]/D1[0][1]
        if all(abs(D1[i][j]*s-D2[p[i]][p[j]])<1e-6 for i in range(5) for j in range(i+1,5)): return True
    return False
found=[]
for col in survivors:
    def res(z):
        Q=[(0,0),(1,0)]+[(z[2*i],z[2*i+1]) for i in range(3)]
        v=[1.0,z[6],z[7]]  # colour 0 edges... scale: set one class to 1? better: free scale
        return [math.dist(Q[i],Q[j])-z[6+col[k]] for k,(i,j) in enumerate(E5)]
    def res2(z):
        Q=[(0,0),(z[0],0)]+[(z[1+2*i],z[2+2*i]) for i in range(3)]
        return [math.dist(Q[i],Q[j])-z[7+col[k]] for k,(i,j) in enumerate(E5)]+[z[7]-1.0]
    for _ in range(150):
        z0=np.concatenate([np.random.uniform(-2,2,7),np.random.uniform(0.3,2.5,3)])
        r=least_squares(res2,z0,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=400)
        if max(abs(r.fun))>1e-9: continue
        z=r.x; Q=[(0,0),(z[0],0)]+[(z[1+2*i],z[2+2*i]) for i in range(3)]
        ds=pdist(Q)
        if min(ds)<1e-5: continue
        vals=sorted(z[7:10])
        if vals[1]-vals[0]<1e-5 or vals[2]-vals[1]<1e-5: continue
        if any(congruent(Q,P) for P,_ in found): continue
        found.append((Q,col)); print('NEW',len(found),[round(v,5) for v in sorted(z[7:10])],col); sys.stdout.flush()
print('TOTAL',len(found))
json.dump([{'pts':Q,'col':c} for Q,c in found],open('five_brute.json','w'))
