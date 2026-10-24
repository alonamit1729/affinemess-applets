"""Systematic search for components of S(n,d) (n-point sets with exactly d distances, mod similarity) via
symmetric parametrizations: cyclic C_m and dihedral D_m orbit structures, optionally with up to two extra
coincidences imposed. Each solution's true local dimension is measured in the full configuration space.
Output: symfam.json  {"n,d": [ {pts, dim, hash, struct, d}, ... ] }"""
import math, itertools, json, sys
import numpy as np
import networkx as nx
from scipy.optimize import least_squares
rng=np.random.default_rng(11)
NMAX=16; DMAX=7; PMAX=5; JMAX=2
def pairs(n): return list(itertools.combinations(range(n),2))
def distinct(vals,tol):
    v=sorted(vals); out=[v[0]]
    for x in v[1:]:
        if x-out[-1]>tol: out.append(x)
    return out
# ---------------- structures ----------------
def structures(n):
    out=[]
    # C_m: c + k*m = n, k orbits of regular m-gons
    for m in range(2,n+1):
        for c in (0,1):
            if (n-c)%m: continue
            k=(n-c)//m
            if k<1: continue
            p=2*k-2
            if p<=PMAX: out.append(('C',m,c,k,0,0,p))
    # D_m, m>=2: c + a*m + b*m + 2m*g = n ; params a+b+2g-1
    for m in range(2,n+1):
        for c in (0,1):
            for g in range(0,n//(2*m)+1):
                rest=n-c-2*m*g
                if rest<0 or rest%m: continue
                ab=rest//m
                for a in range(ab+1):
                    b=ab-a
                    if a+b+g==0: continue
                    p=a+b+2*g-1
                    if p<=PMAX: out.append(('D',m,c,a,b,g,p))
    # D_1 (one mirror): f axis points + p pairs ; params f+2p-2
    for f in range(0,n+1):
        if (n-f)%2: continue
        pp=(n-f)//2
        if pp==0: continue
        p=f+2*pp-2
        if p<=PMAX: out.append(('M',1,f,pp,0,0,p))
    return out
def build(S,z):
    kind,m,c,a,b,g,p=S
    P=[]
    if kind=='C':
        k=a
        if c: P.append((0.0,0.0))
        radii=[1.0]+list(z[0:k-1]); angs=[0.0]+list(z[k-1:2*k-2])
        for r,t in zip(radii,angs):
            for j in range(m): P.append((r*math.cos(t+2*math.pi*j/m),r*math.sin(t+2*math.pi*j/m)))
    elif kind=='D':
        if c: P.append((0.0,0.0))
        vals=[1.0]+list(z)  # first radius normalised
        idx=0
        for _ in range(a):
            r=vals[idx]; idx+=1
            for j in range(m): P.append((r*math.cos(2*math.pi*j/m),r*math.sin(2*math.pi*j/m)))
        for _ in range(b):
            r=vals[idx]; idx+=1
            for j in range(m): P.append((r*math.cos(2*math.pi*j/m+math.pi/m),r*math.sin(2*math.pi*j/m+math.pi/m)))
        for _ in range(g):
            r=vals[idx]; th=vals[idx+1]; idx+=2
            for j in range(m):
                for s in (1,-1): P.append((r*math.cos(2*math.pi*j/m+s*th),r*math.sin(2*math.pi*j/m+s*th)))
    else:  # mirror: f axis points (0,y_i) with y_1=0 ; pairs (x_j,y_j) with x_1=1
        f=c; pp=a
        ys=[0.0]+list(z[0:f-1]) if f>=1 else []
        rest=list(z[max(f-1,0):])
        for y in ys: P.append((0.0,y))
        xs=[1.0]; yy=[]
        # rest = [y_1, x_2, y_2, ..., x_pp, y_pp]  (f>=1) ; if f==0: scale fixed by x_1=1, translation by y_1=0
        if f>=1:
            yy=[rest[0]]; rest=rest[1:]
        else:
            yy=[0.0]
        for j in range(pp-1): xs.append(rest[2*j]); yy.append(rest[2*j+1])
        for x,y in zip(xs,yy): P+=[(x,y),(-x,y)]
    return P
def nparams(S): return S[6]
# ---------------- analysis ----------------
def dists(P): return [math.dist(P[i],P[j]) for i,j in pairs(len(P))]
def colouring(P,tol):
    ds=dists(P); vals=distinct(ds,tol)
    return [min(range(len(vals)),key=lambda k:abs(vals[k]-x)) for x in ds], vals
def type_hash(P,col):
    n=len(P); G=nx.Graph()
    mult={}
    for cidx in col: mult[cidx]=mult.get(cidx,0)+1
    # relabel classes by (multiplicity, sorted degree sequence) to be label-invariant
    E=pairs(n)
    deg={}
    for (i,j),cidx in zip(E,col):
        deg.setdefault(cidx,[0]*n); deg[cidx][i]+=1; deg[cidx][j]+=1
    key={cidx:(mult[cidx],tuple(sorted(deg[cidx]))) for cidx in mult}
    order=sorted(set(key.values()))
    for (i,j),cidx in zip(E,col): G.add_edge(i,j,c=str(order.index(key[cidx])))
    return nx.weisfeiler_lehman_graph_hash(G,edge_attr='c',iterations=4)
def generic_point_and_dim(P,col,d):
    """move to a generic (non-singular) point of the component through P if P is singular; return (point, local dim).
    dim 0 means P is isolated in its stratum."""
    n=len(P); E=pairs(n); NV=2*n-3
    ang=math.atan2(P[1][1]-P[0][1],P[1][0]-P[0][0]); co,si=math.cos(-ang),math.sin(-ang)
    Q=[((x-P[0][0])*co-(y-P[0][1])*si,(x-P[0][0])*si+(y-P[0][1])*co) for x,y in P]
    z=[Q[1][0]]
    for i in range(2,n): z+=[Q[i][0],Q[i][1]]
    vals=[0.0]*d
    for k,(i,j) in enumerate(E): vals[col[k]]=math.dist(Q[i],Q[j])
    z=np.array(z+vals); v0=vals[0]
    def unpack(z):
        R=[(0.0,0.0),(z[0],0.0)]
        for i in range(2,n): R.append((z[2*i-3],z[2*i-2]))
        return R
    def res(z):
        R=unpack(z); v=z[NV:]
        return [math.dist(R[i],R[j])-v[col[k]] for k,(i,j) in enumerate(E)]+[v[0]-v0]
    def jac(z):
        f0=np.array(res(z)); h=1e-6; J=np.zeros((len(f0),len(z)))
        for i in range(len(z)):
            zz=z.copy(); zz[i]+=h; J[:,i]=(np.array(res(zz))-f0)/h
        return J
    def is_singular(z):
        J=jac(z); Jd=J[:len(E),:NV]; sd=np.linalg.svd(Jd,compute_uv=False)
        return int((sd>1e-6*sd[0]).sum())<NV, J
    sing,J=is_singular(z)
    if not sing:
        s=np.linalg.svd(J,compute_uv=False); rank=int((s>1e-6*s[0]).sum())
        return P, len(z)-rank
    # singular: look for nearby solutions
    moved=[]
    for _ in range(16):
        z0=z+rng.normal(0,2e-3,len(z))
        r=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=300)
        if max(abs(r.fun))>1e-11: continue
        if np.abs(r.x-z).max()>1e-4: moved.append(r.x)
    if not moved: return P, 0
    for zz in moved:
        sing2,J2=is_singular(zz)
        if not sing2:
            s=np.linalg.svd(J2,compute_uv=False); rank=int((s>1e-6*s[0]).sum())
            return unpack(zz), len(zz)-rank
    # all nearby solutions singular too: estimate from displacement span (may overestimate at crossings)
    A=np.array([zz-z for zz in moved]); return P, int((np.linalg.svd(A,compute_uv=False)>3e-4).sum())
def sim_key(P):
    ds=sorted(dists(P)); m=ds[0]; return tuple(round(x/m,6) for x in ds)
results={}
def record(n,d,P,S,extra=''):
    col,vals=colouring(P,1e-7*max(dists(P)))
    if len(vals)!=d: return False
    h0=type_hash(P,col); key=f'{n},{d}'; L=results.setdefault(key,[])
    sk0=sim_key(P)
    for e in L:
        if e['hash']==h0 and (e['dim']>=1 or e.get('sk')==sk0): return False
    G,dim=generic_point_and_dim(P,col,d)
    colG,valsG=colouring(G,1e-7*max(dists(G)))
    if len(valsG)!=d or not ok_config(G): G,colG=P,col
    h=type_hash(G,colG)
    key=f'{n},{d}'
    L=results.setdefault(key,[])
    if dim==0:
        sk=sim_key(P)
        if any(e['dim']==0 and e.get('sk')==sk for e in L): return False
        L.append({'pts':[list(p) for p in P],'dim':0,'hash':h,'sk':sk,'struct':str(S)+extra,'d':d})
    else:
        if any(e['hash']==h and e['dim']==dim for e in L): return False
        L.append({'pts':[list(p) for p in G],'dim':dim,'hash':h,'struct':str(S)+extra,'d':d})
    return True
def ok_config(P):
    ds=dists(P); return min(ds)>1e-3*max(ds)
def process(n,S):
    p=nparams(S)
    # generic sample(s)
    best=None
    for _ in range(6):
        z=rng.uniform(0.3,2.0,p) if p else np.zeros(0)
        # angles for D generic orbits must be in (0, pi/m); C angles in (0, 2pi/m)
        z=randomize(S,z)
        P=build(S,z)
        if not ok_config(P): continue
        col,vals=colouring(P,1e-7*max(dists(P)))
        if best is None or len(vals)>best[0]: best=(len(vals),z,P,col)
    if best is None: return
    dgen,z,P,col=best
    # j extra coincidences, j=0..JMAX, with p-j>=0 and dgen-j<=DMAX
    for j in range(0,min(JMAX,p)+1):
        d=dgen-j
        if d>DMAX or d<1: continue
        if j==0:
            record(n,d,P,S); continue
        # choose j disjoint-or-chained coincidence pairs among the dgen classes
        for pat in itertools.combinations(itertools.combinations(range(dgen),2),j):
            # merged class mapping
            parent=list(range(dgen))
            def find(x):
                while parent[x]!=x: x=parent[x]
                return x
            for a,b in pat: parent[find(a)]=find(b)
            roots=sorted(set(find(c) for c in range(dgen)))
            if len(roots)!=d: continue
            newcol=[roots.index(find(c)) for c in col]
            def res(w):
                Q=build(S,w[:p]); v=w[p:]
                return [math.dist(Q[i],Q[j_])-v[newcol[k]] for k,(i,j_) in enumerate(pairs(n))]
            found=0
            for _ in range(6):
                w0=np.concatenate([randomize(S,rng.uniform(0.2,2.5,p)),rng.uniform(0.3,3,d)])
                r=least_squares(res,w0,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=300)
                if max(abs(r.fun))>1e-11: continue
                Q=build(S,r.x[:p])
                if not ok_config(Q): continue
                colQ,valsQ=colouring(Q,1e-7*max(dists(Q)))
                if len(valsQ)!=d: continue
                if record(n,d,Q,S,extra=' +%d coincidence(s)'%j): found+=1
                if found>=4: break
def randomize(S,z):
    kind,m,c,a,b,g,p=S
    z=np.array(z,float)
    if kind=='C':
        k=a
        for i in range(k-1): z[k-1+i]=rng.uniform(0.05,2*math.pi/m-0.05)
    elif kind=='D':
        idx=a+b-1
        for _ in range(g):
            if idx>=0: z[idx]=rng.uniform(0.3,2.0)
            z[idx+1]=rng.uniform(0.05,math.pi/m-0.05); idx+=2
    else:
        z=rng.uniform(-2,2,p)
    return z
if __name__=='__main__':
    ns=[int(x) for x in sys.argv[1:]] or list(range(4,NMAX+1))
    for n in ns:
        SS=structures(n)
        print('n',n,'structures',len(SS),file=sys.stderr)
        for S in SS:
            try: process(n,S)
            except Exception as e: print('ERR',S,e,file=sys.stderr)
        json.dump(results,open('symfam_%d.json'%n,'w'))
        summ={k:{} for k in results if k.startswith(f'{n},')}
        for k in summ:
            for e in results[k]: summ[k][e['dim']]=summ[k].get(e['dim'],0)+1
        summ={k:dict(sorted(v.items())) for k,v in summ.items()}
        print(n,summ)
