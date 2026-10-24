"""Recompute the dimension at the stored sample points of a strata file with a robust perturbation test
(the original run misjudged collinear configurations). Rewrites the file in place."""
import json, math, itertools, sys
import numpy as np
from scipy.optimize import least_squares
fn=sys.argv[1]
D=json.load(open(fn)); n=D['n']; d=D['d']
E=list(itertools.combinations(range(n),2)); NV=2*n-3
rng=np.random.default_rng(123)
def unpack(z):
    P=[(0.0,0.0),(z[0],0.0)]
    for i in range(2,n): P.append((z[2*i-3],z[2*i-2]))
    return P
def pack(P):
    # P[0] at origin, P[1] on x-axis assumed
    z=[P[1][0]]
    for i in range(2,n): z+= [P[i][0],P[i][1]]
    return z
def make_res(col):
    def res(z):
        P=unpack(z); v=z[NV:NV+d]
        r=[math.dist(P[i],P[j])-v[col[k]] for k,(i,j) in enumerate(E)]; r.append(v[0]-1.0); return r
    return res
def jac_rank(res,z):
    h=1e-6; f0=np.array(res(z)); J=np.zeros((len(f0),len(z)))
    for i in range(len(z)):
        zz=z.copy(); zz[i]+=h; J[:,i]=(np.array(res(zz))-f0)/h
    s=np.linalg.svd(J,compute_uv=False); return int((s>1e-5*s[0]).sum())
def dist_jac_rank(z):
    h=1e-6; P=unpack(z); D0=np.array([math.dist(P[i],P[j]) for i,j in E]); J=np.zeros((len(E),NV))
    for i in range(NV):
        zz=z.copy(); zz[i]+=h; Q=unpack(zz); J[:,i]=(np.array([math.dist(Q[a],Q[b]) for a,b in E])-D0)/h
    s=np.linalg.svd(J,compute_uv=False); return int((s>1e-6*s[0]).sum())
def local_dim(res,z):
    disp=[]
    for _ in range(16):
        z0=z+rng.normal(0,2e-3,len(z))
        r=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=400)
        if max(abs(r.fun))>1e-11: continue
        disp.append(r.x-z)
    if not disp: return None
    A=np.array(disp)
    if np.abs(A).max()<1e-4: return 0
    s=np.linalg.svd(A,compute_uv=False)
    return int((s>3e-4).sum())
def dim_at(res,z):
    if dist_jac_rank(z)<NV: return local_dim(res,z)
    return NV+d-jac_rank(res,z)
def dm(P): return [[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(n)):
        s=D2[p[0]][p[1]]/D1[0][1]
        if all(abs(D1[i][j]*s-D2[p[i]][p[j]])<1e-6 for i in range(n) for j in range(i+1,n)): return True
    return False
changed=0
for r in D['realizable']:
    col=r['col']; res=make_res(col)
    newsample={}
    for k,P in list(r['sample'].items()):
        P=[tuple(p) for p in P]
        # values
        vals=[0.0]*d
        for kk,(i,j) in enumerate(E): vals[col[kk]]=math.dist(P[i],P[j])
        z=np.array(pack(P)+vals,float)
        # normalise scale: v0=1
        s=vals[0]; z[:NV]/=s; z[NV:]/=s
        dim=dim_at(res,z)
        if dim is None: dim=int(k)
        if dim!=int(k): changed+=1
        if dim==0:
            if not any(congruent(P,Q) for Q in r['rigid']): r['rigid'].append([list(p) for p in P])
        else:
            if str(dim) not in newsample: newsample[str(dim)]=[list(p) for p in P]
    r['sample']=newsample
    r['dims']={**({'0':len(r['rigid'])} if r['rigid'] else {}),**{k:1 for k in newsample}}
json.dump(D,open(fn,'w'))
tot=sum(len(r['rigid']) for r in D['realizable'])
fams={}
for r in D['realizable']:
    for k in r['sample']: fams[k]=fams.get(k,0)+1
print(fn,'changed',changed,'rigid',tot,'family types by dim',fams)
