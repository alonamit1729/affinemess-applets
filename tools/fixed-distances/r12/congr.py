import math, numpy as np
def congruent(P,Q,tol=1e-6):
    """exact similarity test: a vertex bijection matching all pairwise distances (after scaling), by backtracking"""
    n=len(P)
    if len(Q)!=n: return False
    D1=np.array([[math.dist(a,b) for b in P] for a in P]); D2=np.array([[math.dist(a,b) for b in Q] for a in Q])
    m1=D1[D1>0].min(); m2=D2[D2>0].min(); D1=D1/m1; D2=D2/m2; T=tol*D2.max()
    s1=np.sort(D1[np.triu_indices(n,1)]); s2=np.sort(D2[np.triu_indices(n,1)])
    if np.abs(s1-s2).max()>T: return False
    r1=[np.sort(D1[i]) for i in range(n)]; r2=[np.sort(D2[i]) for i in range(n)]
    cand=[[j for j in range(n) if np.abs(r1[i]-r2[j]).max()<T] for i in range(n)]
    order=sorted(range(n),key=lambda i: len(cand[i]))
    perm=[-1]*n; used=[False]*n
    def bt(k):
        if k==n: return True
        i=order[k]
        for j in cand[i]:
            if used[j]: continue
            if all(abs(D1[i][order[t]]-D2[j][perm[order[t]]])<T for t in range(k)):
                perm[i]=j; used[j]=True
                if bt(k+1): return True
                used[j]=False; perm[i]=-1
        return False
    return bt(0)
