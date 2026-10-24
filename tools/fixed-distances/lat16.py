import math, itertools, json, numpy as np, networkx as nx, sys, time
s3=math.sqrt(3)
R=float(sys.argv[1]) if len(sys.argv)>1 else 4.01
pts=[(a+b/2,b*s3/2) for a in range(-8,9) for b in range(-8,9)]
pts=[p for p in pts if math.hypot(*p)<=R]
def d2(p,q): return round((p[0]-q[0])**2+(p[1]-q[1])**2)
def key(P):
    D=np.array([[math.dist(a,b) for b in P] for a in P]); return tuple(sorted(tuple(np.round(sorted(r),5)) for r in D))
found={}
cands=[1,3,4,7,9,12,13,16]
t0=time.time()
for D7 in itertools.combinations(cands,7):
    G=nx.Graph(); G.add_nodes_from(range(len(pts)))
    for i,j in itertools.combinations(range(len(pts)),2):
        if d2(pts[i],pts[j]) in D7: G.add_edge(i,j)
    nmax=0
    for c in nx.find_cliques(G):
        nmax=max(nmax,len(c))
        if len(c)<16: continue
        for sub in itertools.combinations(c,16):
            P=[pts[i] for i in sub]
            if len({d2(a,b) for a,b in itertools.combinations(P,2)})==7:
                k=key(P)
                if k not in found: found[k]=P; print('found',D7,len(found),flush=True)
    print(D7,'max clique',nmax,'%.0fs'%(time.time()-t0),flush=True)
json.dump(list(found.values()),open('lattice_16_7.json','w'))
print('total',len(found))
