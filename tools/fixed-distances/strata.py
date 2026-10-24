"""Enumerate combinatorial types (set partitions of the edges of K_n up to S_n) for n points with
exactly d distances, test realizability in the plane by random-start least squares, record the
dimension of the realization space at each solution and cluster rigid solutions by congruence.
Usage: python3 strata.py n d seed starts"""
import itertools, math, sys, json, random
import numpy as np
from scipy.optimize import least_squares
n=int(sys.argv[1]); d=int(sys.argv[2]); seed=int(sys.argv[3]); STARTS=int(sys.argv[4])
rng=np.random.default_rng(seed)
E=list(itertools.combinations(range(n),2)); m=len(E)
# --- set partitions of E into d blocks, canonical under S_n ---
def partitions(items,k):
    # restricted growth strings
    def rec(i,cur,mx):
        if i==len(items):
            if mx==k: yield tuple(cur)
            return
        for b in range(min(mx+1,k)):
            if b==mx and k-mx > len(items)-i: pass
            cur.append(b); yield from rec(i+1,cur,max(mx,b+1)); cur.pop()
    yield from rec(0,[],0)
perms=[]
for p in itertools.permutations(range(n)):
    perms.append([E.index(tuple(sorted((p[a],p[b])))) for a,b in E])
def canon(rgs):
    best=None
    for pe in perms:
        c=[rgs[pe[i]] for i in range(m)]
        # relabel to restricted growth
        mp={}; out=[]
        for x in c:
            if x not in mp: mp[x]=len(mp)
            out.append(mp[x])
        t=tuple(out)
        if best is None or t<best: best=t
    return best
classes=set()
for rgs in partitions(E,d):
    classes.add(canon(rgs))
classes=sorted(classes)
print('n',n,'d',d,'classes',len(classes),file=sys.stderr)
# --- realizability ---
NV=2*n-3  # coordinates: x1; (x_i,y_i) i>=2
def unpack(z):
    P=[(0.0,0.0),(z[0],0.0)]
    for i in range(2,n): P.append((z[2*i-3],z[2*i-2]))
    return P
def make_res(col):
    def res(z):
        P=unpack(z); v=z[NV:NV+d]
        r=[math.dist(P[i],P[j])-v[col[k]] for k,(i,j) in enumerate(E)]
        r.append(v[0]-1.0)
        return r
    return res
def jac_rank(res,z):
    h=1e-6; f0=np.array(res(z)); J=np.zeros((len(f0),len(z)))
    for i in range(len(z)):
        zz=z.copy(); zz[i]+=h; J[:,i]=(np.array(res(zz))-f0)/h
    s=np.linalg.svd(J,compute_uv=False)
    return int((s>1e-5*s[0]).sum())
def dm(P): return [[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(n)):
        s=D2[p[0]][p[1]]/D1[0][1]
        ok=True
        for i in range(n):
            for j in range(i+1,n):
                if abs(D1[i][j]*s-D2[p[i]][p[j]])>1e-6: ok=False;break
            if not ok: break
        if ok: return True
    return False
def dist_jac_rank(z):
    # rank of the Jacobian of the distance map itself (collinear sets are singular points of it)
    h=1e-6; P=unpack(z); D0=np.array([math.dist(P[i],P[j]) for i,j in E]); J=np.zeros((m,NV))
    for i in range(NV):
        zz=z.copy(); zz[i]+=h; Q=unpack(zz); J[:,i]=(np.array([math.dist(Q[a],Q[b]) for a,b in E])-D0)/h
    s=np.linalg.svd(J,compute_uv=False); return int((s>1e-6*s[0]).sum())
def local_dim(res,z):
    """dimension of the solution set near z, by perturb-and-resolve (used at singular points)"""
    disp=[]
    for _ in range(14):
        z0=z+rng.normal(0,2e-3,len(z))
        r=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=400)
        if max(abs(r.fun))>1e-11: continue
        disp.append(r.x-z)
    if not disp: return None
    A=np.array(disp)
    if np.abs(A).max()<1e-4: return 0
    s=np.linalg.svd(A,compute_uv=False)
    return int((s>3e-4).sum())
results=[]
for ci,col in enumerate(classes):
    res=make_res(col)
    rigid=[]; dims={}
    sample={}
    for _ in range(STARTS):
        z0=np.concatenate([rng.uniform(-2,2,NV),rng.uniform(0.3,2.5,d)])
        r=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=400)
        if max(abs(r.fun))>1e-11: continue
        z=r.x; P=unpack(z); v=sorted(z[NV:NV+d])
        dd=[math.dist(P[i],P[j]) for i,j in E]
        if min(dd)<1e-3*max(dd): continue
        if any(v[k+1]-v[k]<1e-3*v[-1] for k in range(d-1)): continue
        if dist_jac_rank(z)<NV:
            dim=local_dim(res,z)
            if dim is None: continue
        else:
            dim=NV+d-jac_rank(res,z)
        dims[dim]=dims.get(dim,0)+1
        if dim==0:
            if not any(congruent(P,Q) for Q in rigid): rigid.append(P)
        else:
            if dim not in sample: sample[dim]=P
    if dims:
        results.append({'col':list(col),'dims':dims,'rigid':[[list(p) for p in P] for P in rigid],'sample':{str(k):[list(p) for p in P] for k,P in sample.items()}})
        print(ci,col,dims,'rigid',len(rigid),file=sys.stderr)
json.dump({'n':n,'d':d,'classes':len(classes),'realizable':results},open(f'strata_{n}_{d}_{seed}.json','w'))
tot_rigid=sum(len(r['rigid']) for r in results)
print('n',n,'d',d,'classes',len(classes),'realizable',len(results),'rigid configs',tot_rigid,'dims',sorted(set(k for r in results for k in r['dims'])))
