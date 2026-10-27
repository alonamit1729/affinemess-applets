"""Independent enumeration of 11- and 12-point subsets of the triangular lattice (and of the square lattice)
with at most six distances, up to similarity.  Scope: minimum squared distance m over all lattice norms m <= MMAX,
and diameter^2 <= RATIO2 * m.  Method: fix one point at the origin and a point p1 realising the minimum distance
(one representative per symmetry orbit of norm m); DFS over the remaining candidates (norm in [m, RATIO2*m] from
every chosen point), keeping at most six distinct squared distances.  Output: similarity classes found, compared with
the one-point deletions of the eleven twelve-point cores of the paper (cores.json) when run on the triangular lattice."""
import sys, json, math, itertools, time
import numpy as np
lat=sys.argv[1] if len(sys.argv)>1 else 'tri'
MMAX=int(sys.argv[2]) if len(sys.argv)>2 else 16
RATIO2=int(sys.argv[3]) if len(sys.argv)>3 else 25
SIZES=tuple(int(x) for x in sys.argv[4].split(',')) if len(sys.argv)>4 else (11,12)
KMAX=max(SIZES); KMIN=min(SIZES)
if lat=='tri':
    Q=lambda p: p[0]*p[0]+p[0]*p[1]+p[1]*p[1]
    eu=lambda p: (p[0]+p[1]/2, p[1]*math.sqrt(3)/2)
    # symmetry group D6 on Eisenstein coordinates (a,b) ~ a + b*omega: units and conjugation
    def syms(p):
        a,b=p; out=set(); x=(a,b)
        for _ in range(6):
            out.add(x); out.add((x[0]+x[1],-x[1]))   # conjugation: a+b w -> a+b w^2 = (a-b) - b w ... use (a+b,-b)
            x=(-x[1],x[0]+x[1])                       # multiply by omega: (a+bw)w = aw + b w^2 = -b + (a-b) w
        return out
else:
    Q=lambda p: p[0]*p[0]+p[1]*p[1]
    eu=lambda p: (float(p[0]),float(p[1]))
    def syms(p):
        a,b=p; out=set()
        for x in [(a,b),(-b,a),(-a,-b),(b,-a)]:
            out.add(x); out.add((x[0],-x[1]))
        return out
def key(P):
    E=[eu(p) for p in P]
    D=np.array([[math.dist(a,b) for b in E] for a in E]); D/=D[D>0].min()
    return tuple(sorted(tuple(np.round(sorted(r),6)) for r in D))
R=int(math.isqrt(RATIO2*MMAX))+2
allpts=[(a,b) for a in range(-R,R+1) for b in range(-R,R+1) if 0<Q((a,b))<=RATIO2*MMAX]
norms=sorted({Q(p) for p in allpts if Q(p)<=MMAX})
found={k:{} for k in SIZES}
t0=time.time()
for m in norms:
    # representatives of norm-m points up to symmetry
    reps=[]; seen=set()
    for p in allpts:
        if Q(p)==m and p not in seen:
            reps.append(p); seen|=syms(p)
    cand=[p for p in allpts if m<=Q(p)<=RATIO2*m]
    for p1 in reps:
        C=[p for p in cand if p!=p1 and m<=Q((p[0]-p1[0],p[1]-p1[1]))<=RATIO2*m]
        n=len(C)
        # pairwise norms among candidates
        A=np.array(C);
        dx=A[:,0][:,None]-A[:,0][None,:]; dy=A[:,1][:,None]-A[:,1][None,:]
        if lat=='tri': NM=dx*dx+dx*dy+dy*dy
        else: NM=dx*dx+dy*dy
        ok=(NM>=m)&(NM<=RATIO2*m)
        q0=np.array([Q(p) for p in C]); q1=np.array([Q((p[0]-p1[0],p[1]-p1[1])) for p in C])
        nodes=0
        def dfs(chosen, pal, avail):
            global nodes
            nodes+=1
            k=len(chosen)+2
            if k in found:
                P=[(0,0),p1]+[C[i] for i in chosen]
                kk=key(P)
                if kk not in found[k]: found[k][kk]=P
            if k==KMAX: return
            if k+len(avail)<KMIN: return
            for t,i in enumerate(avail):
                newpal=pal|{q0[i],q1[i]}|{NM[i,j] for j in chosen}
                if len(newpal)>6: continue
                # restrict later candidates: index > i, compatible with i
                nav=[j for j in avail[t+1:] if ok[i,j]]
                dfs(chosen+[i], newpal, nav)
        dfs([], {m}, list(range(n)))
    print(f'm={m}: reps {len(reps)}, cand {len(cand)}, nodes {nodes}, sets so far {[len(found[k]) for k in SIZES]}, {time.time()-t0:.0f}s',flush=True)
json.dump({str(k):list(v.values()) for k,v in found.items()},open(f'latt_{lat}_{MMAX}_{RATIO2}_{KMAX}.json','w'))
if lat=='tri' and SIZES==(11,12):
    cores12=json.load(open('cores.json'))
    ck={}
    for L in cores12:
        k=key([tuple(p) for p in L]); ck.setdefault(k,[]).append(L)
    c11={}
    for L in cores12:
        for i in range(12):
            A=[tuple(p) for j,p in enumerate(L) if j!=i]
            if len({Q((a[0]-b[0],a[1]-b[1])) for a,b in itertools.combinations(A,2)})<=6: c11.setdefault(key(A),A)
    print('paper: 12-point core classes',len(ck),' 11-point deletion classes',len(c11))
    print('found 12-point classes',len(found[12]),' in paper:',sum(1 for k in found[12] if k in ck),' paper in found:',sum(1 for k in ck if k in found[12]))
    print('found 11-point classes',len(found[11]),' in paper:',sum(1 for k in found[11] if k in c11),' paper in found:',sum(1 for k in c11 if k in found[11]))
    extra=[P for k,P in found[11].items() if k not in c11]
    if extra: print('EXTRA 11-point lattice sets not among the paper deletions:',extra)
    extra=[P for k,P in found[12].items() if k not in ck]
    if extra: print('EXTRA 12-point lattice sets not among the paper cores:',extra)
    for P in found[11].values():
        S=sorted({Q((a[0]-b[0],a[1]-b[1])) for a,b in itertools.combinations(P,2)})
    print('distance-count histogram of found 11-sets:',
          {c:sum(1 for P in found[11].values() if len({Q((a[0]-b[0],a[1]-b[1])) for a,b in itertools.combinations(P,2)})==c) for c in (4,5,6)})
