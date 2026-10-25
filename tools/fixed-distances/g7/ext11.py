"""Independent replication of Lemma 7 (lattice cores): every 11-point lattice core A with <=6 distances, extended by
points whose distances to A lie in S plus at most one (six-distance cores) or two (five-distance cores) new values,
with the largest new value the diameter: maximum total size.  Exact rational / sympy arithmetic."""
import json, itertools, math
from fractions import Fraction as F
import sympy as sp
cores12=json.load(open('cores.json'))
Q=lambda p: p[0]*p[0]+p[0]*p[1]+p[1]*p[1]
d=lambda p,q: Q((p[0]-q[0],p[1]-q[1]))
def eu(p): return (p[0]+p[1]/2, p[1]*math.sqrt(3)/2)
def key(P):
    import numpy as np
    D=np.array([[math.dist(eu(a),eu(b)) for b in P] for a in P]); D/=D[D>0].min()
    return tuple(sorted(tuple(np.round(sorted(r),6)) for r in D))
cores={}
for L in cores12:
    for i in range(12):
        A=[tuple(p) for j,p in enumerate(L) if j!=i]
        S=sorted({d(a,b) for a,b in itertools.combinations(A,2)})
        if len(S)<=6:
            k=key(A)
            if k not in cores: cores[k]=(A,S)
print('11-point lattice cores (similarity classes):',len(cores),'with 5 distances:',sum(1 for A,S in cores.values() if len(S)==5))
def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(P):
    P=sorted(set(P)); lo=[]
    for p in P:
        while len(lo)>1 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    up=[]
    for p in reversed(P):
        while len(up)>1 and cross(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return set(lo[:-1]+up[:-1])
def maxclique(V,adj):
    best=0
    def go(cand,size):
        nonlocal best
        if size+len(cand)<=best: return
        if not cand: best=max(best,size); return
        c=cand[0]
        go([x for x in cand[1:] if adj(c,x)],size+1); go(cand[1:],size)
    go(list(V),0); return best
# exact arithmetic for q = (x,y) in triangular coordinates; distance squared Q(q-p) is a rational (or u+v*sqrt(D))
results=[]
for k,(A,S) in cores.items():
    H=hull(A); I=[p for p in A if p not in H]; Sset=set(S); mS=max(S)
    if len(S)==6:
        anchors=next(t for t in itertools.combinations(I,3) if cross(*t)!=0)
        A0,B0,C0=anchors; w=(B0[0]-A0[0],B0[1]-A0[1]); v=(C0[0]-A0[0],C0[1]-A0[1])
        # Q(q-A0)=a, Q(q-B0)=b, Q(q-C0)=c; with r=q-A0: Q(r)-Q(r-w)=a-b -> 2 r.w (in Q-form) - Q(w) = a-b
        # Q-form bilinear: <r,w> = r0 w0 + (r0 w1 + r1 w0)/2 + r1 w1 ; Q(r-w)=Q(r)-2<r,w>+Q(w)
        cands=set()
        for a,b,c in itertools.product(S,repeat=3):
            # 2<r,w> = a-b+Q(w); 2<r,v> = a-c+Q(v): linear system in r=(x,y)
            M=[[2*w[0]+w[1], w[0]+2*w[1]],[2*v[0]+v[1], v[0]+2*v[1]]]; rhs=[a-b+Q(w),a-c+Q(v)]
            det=M[0][0]*M[1][1]-M[0][1]*M[1][0]
            x=F(rhs[0]*M[1][1]-M[0][1]*rhs[1],det); y=F(M[0][0]*rhs[1]-rhs[0]*M[1][0],det)
            q=(A0[0]+x,A0[1]+y)
            if d(q,A0)==a and d(q,B0)==b and d(q,C0)==c: cands.add(q)
        cands=[q for q in cands if q not in set(A) and all(d(q,p) in Sset for p in I)]
        surv=[]
        for q in cands:
            new={d(q,p) for p in A}-Sset
            if len(new)<=1 and all(t>mS for t in new): surv.append(q)
        pool={d(q,p) for q in surv for p in A}-Sset
        pool|={d(p,q) for p,q in itertools.combinations(surv,2)}-Sset
        pool={t for t in pool if t>mS}
        best=0
        for T in [None]+sorted(pool):
            pal=Sset|({T} if T else set())
            V=[q for q in surv if all(d(q,p) in pal for p in A)]
            best=max(best,maxclique(V,lambda p,q: d(p,q) in pal))
        results.append((len(S),len(I),len(cands),len(surv),11+best)); continue
    # five-distance cores: exact circle intersections (sympy) and circumcenters
    assert len(I)>=3
    x,y=sp.symbols('x y')
    def Qs(p): return p[0]**2+p[0]*p[1]+p[1]**2
    cands=set()
    for (A0,B0) in itertools.combinations(I,2):
        for a,b in itertools.product(S,repeat=2):
            sols=sp.solve([Qs((x-A0[0],y-A0[1]))-a,Qs((x-B0[0],y-B0[1]))-b],[x,y],dict=True)
            for s in sols:
                qx,qy=sp.nsimplify(s[x]),sp.nsimplify(s[y])
                if qx.is_real and qy.is_real: cands.add((sp.simplify(qx),sp.simplify(qy)))
    for t in itertools.combinations(I,3):
        if cross(*t)==0: continue
        sols=sp.solve([Qs((x-t[0][0],y-t[0][1]))-Qs((x-t[1][0],y-t[1][1])),Qs((x-t[0][0],y-t[0][1]))-Qs((x-t[2][0],y-t[2][1]))],[x,y],dict=True)
        for s in sols: cands.add((sp.simplify(s[x]),sp.simplify(s[y])))
    def ds(q,p): return sp.simplify(Qs((q[0]-p[0],q[1]-p[1])))
    cands=[q for q in cands if not any(q[0]==p[0] and q[1]==p[1] for p in A)]
    surv=[]
    for q in cands:
        newI={ds(q,p) for p in I}-Sset; newA={ds(q,p) for p in A}-Sset
        if len(newI)<=1 and len(newA)<=2: surv.append(q)
    assert all(v.is_rational for q in surv for v in q), 'irrational survivor'
    pool={ds(q,p) for q in surv for p in A}-Sset; pool|={ds(p,q) for p,q in itertools.combinations(surv,2)}-Sset
    best=0; bestpal=None
    for a,b in itertools.combinations(sorted(pool),2):
        if b<=mS: continue
        pal=Sset|{a,b}
        V=[q for q in surv if all(ds(q,p) in pal for p in A) and all(ds(q,p)!=b for p in I)]
        c=maxclique(V,lambda p,q: ds(p,q) in pal)
        if 11+c>best: best=11+c; bestpal=(a,b)
    for a in sorted(pool):
        pal=Sset|{a}; V=[q for q in surv if all(ds(q,p) in pal for p in A)]
        best=max(best,11+maxclique(V,lambda p,q: ds(p,q) in pal))
    results.append((len(S),len(I),len(cands),len(surv),best,bestpal))
print('core results (ndist, |I|, candidates, survivors, max |X|):')
for r in results: print('  ',r)
print('MAX over all lattice cores:',max(r[4] for r in results))
