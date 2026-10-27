"""All subsets of the Minkowski sum R12+R12 = {z^a+z^b : z^12=1} (73 points) with at most D distinct distances,
up to similarity.  Usage: python3 sweep.py D [nmin]   -> r12_subsets_D.json: {"n,d": [ {pts, dists} ... ]}.
DFS in increasing point index; first point restricted to one representative per D12-orbit; dedup by canonical
sorted-distance-matrix key (float, 6 digits, normalised by shortest distance)."""
import sys, math, cmath, itertools, json, time, collections
D=int(sys.argv[1]); NMIN=int(sys.argv[2]) if len(sys.argv)>2 else 3
Z=[cmath.exp(2j*math.pi*k/12) for k in range(12)]
pts=sorted({(round((Z[a]+Z[b]).real,9),round((Z[a]+Z[b]).imag,9)) for a in range(12) for b in range(12)})
N=len(pts); assert N==73
# orbit representatives of D12: group points by radius
rad=[round(math.hypot(*p),6) for p in pts]
reps=[]; seen=set()
for i,r in enumerate(rad):
    if r not in seen: seen.add(r); reps.append(i)
# squared distances -> bit index
d2={}
dist=[[0]*N for _ in range(N)]
for i in range(N):
    for j in range(N):
        if i==j: continue
        v=round((pts[i][0]-pts[j][0])**2+(pts[i][1]-pts[j][1])**2,7)
        if v not in d2: d2[v]=len(d2)
        dist[i][j]=1<<d2[v]
print('points',N,'distinct distances in the whole set',len(d2),flush=True)
vals=sorted(d2,key=d2.get)
out=collections.defaultdict(dict)   # (n,d) -> key -> record
cnt=collections.Counter(); nodes=[0]; t0=time.time()
def key(S):
    P=[pts[i] for i in S]
    Dm=[[math.dist(a,b) for b in P] for a in P]; m=min(x for r in Dm for x in r if x>0)
    return tuple(sorted(tuple(round(x/m,6) for x in sorted(r)) for r in Dm))
# key above is only a fingerprint (sorted rows); good enough for dedup at these sizes, confirmed later by exact compare
def rec(S,mask):
    nodes[0]+=1
    n=len(S); d=bin(mask).count('1')
    if n>=NMIN:
        k=key(S); cnt[(n,d)]+=1
        if k not in out[(n,d)]:
            out[(n,d)][k]={'pts':[pts[i] for i in S],'d2':[vals[b] for b in range(len(vals)) if mask>>b&1]}
    last=S[-1]
    for p in range(last+1,N):
        m=mask
        for q in S:
            m|=dist[p][q]
            if bin(m).count('1')>D: break
        else:
            S.append(p); rec(S,m); S.pop()
for r in reps:
    rec([r],0)
    print('rep',r,'radius',rad[r],'nodes',nodes[0],'%.0fs'%(time.time()-t0),flush=True)
res={}
for (n,d),dd in sorted(out.items()):
    res[f'{n},{d}']=list(dd.values())
    print(f'({n},{d}): raw {cnt[(n,d)]} distinct {len(dd)}')
json.dump(res,open(f'r12_subsets_{D}.json','w'))
print('nodes',nodes[0],'time %.0fs'%(time.time()-t0))
