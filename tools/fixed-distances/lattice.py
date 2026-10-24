"""Find n-point subsets of the triangular / square lattice (and of regular polygons + centre) with exactly d distances,
by clique search on candidate point sets. Usage: python3 lattice.py n d"""
import itertools, math, sys, json
import networkx as nx
n=int(sys.argv[1]); d=int(sys.argv[2])
S3=math.sqrt(3)
def tri_pts(R):
    pts=[]
    for a in range(-R,R+1):
        for b in range(-R,R+1):
            x=a+b/2; y=b*S3/2
            if x*x+y*y<=R*R+1e-9: pts.append((x,y))
    return pts
def sq_pts(R):
    return [(a,b) for a in range(-R,R+1) for b in range(-R,R+1) if a*a+b*b<=R*R+1e-9]
def dkey(p,q): return round(math.dist(p,q)**2,6)
def search(pts,label):
    # distances available
    sq=sorted(set(dkey(p,q) for p,q in itertools.combinations(pts,2)))
    found=[]
    # the set must contain the smallest distance? not necessarily; but it contains some distance; WLOG scale so that... lattice has no scaling freedom.
    # choose the d-subset D of squared distances; graph: vertices pts, edge if dkey in D; find cliques of size n; keep those using exactly d distances.
    small=[s for s in sq if s<=sq[min(len(sq)-1, 14)]]  # restrict to small distances (diameter bounded)
    for D in itertools.combinations(small,d):
        Dset=set(D)
        G=nx.Graph(); G.add_nodes_from(range(len(pts)))
        for i,j in itertools.combinations(range(len(pts)),2):
            if dkey(pts[i],pts[j]) in Dset: G.add_edge(i,j)
        for C in nx.find_cliques(G):
            if len(C)>=n:
                for S in itertools.combinations(sorted(C),n):
                    P=[pts[i] for i in S]
                    used=set(dkey(p,q) for p,q in itertools.combinations(P,2))
                    if len(used)==d: found.append(P)
    # dedupe by congruence
    def dm(P): return sorted(round(math.dist(p,q),6) for p,q in itertools.combinations(P,2))
    def cong(P,Q):
        D1=[[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
        D2=[[math.dist(Q[i],Q[j]) for j in range(n)] for i in range(n)]
        if dm(P)!=dm(Q): return False
        for p in itertools.permutations(range(n)):
            if all(abs(D1[i][j]-D2[p[i]][p[j]])<1e-6 for i in range(n) for j in range(i+1,n)): return True
        return False
    uniq=[]
    for P in found:
        if not any(cong(P,Q) for Q in uniq): uniq.append(P)
    print(label,'found',len(uniq))
    for P in uniq: print('   dists',sorted(set(round(math.dist(p,q),4) for p,q in itertools.combinations(P,2))),[(round(x,3),round(y,3)) for x,y in P])
    return uniq
R=3 if n<=9 else 4
if n>=12: R=3
tri=search(tri_pts(R),'triangular')
sqr=search(sq_pts(3),'square')
json.dump({'tri':tri,'sq':sqr},open(f'lattice_{n}_{d}.json','w'))
