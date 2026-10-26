"""S_n-invariants of the cohomology of F(C,n) (Orlik-Solomon algebra of the braid arrangement) and of
M_{0,n+1} = F(C,n)/Aff(C); hence rational Betti numbers of the unlabelled shape spaces."""
import itertools, numpy as np, sys
from fractions import Fraction
n=int(sys.argv[1]) if len(sys.argv)>1 else 5
pairs=list(itertools.combinations(range(n),2)); N=len(pairs); pidx={p:i for i,p in enumerate(pairs)}
def perm_on_pairs(g):
    """g acts on generators e_{ij} -> e_{g(i)g(j)} (no sign: the hyperplane z_i=z_j maps to z_{g i}=z_{g j})"""
    return [pidx[tuple(sorted((g[a],g[b])))] for a,b in pairs]
def wedge_basis(k): return list(itertools.combinations(range(N),k))
def act_wedge(gp,k):
    B=wedge_basis(k); bi={b:i for i,b in enumerate(B)}; M=np.zeros((len(B),len(B)))
    for j,b in enumerate(B):
        img=[gp[x] for x in b]; order=sorted(range(k),key=lambda t: img[t]); sign=np.linalg.det(np.eye(k)[order]) if k>1 else 1.0
        M[bi[tuple(sorted(img))],j]=sign
    return M,B,bi
def relations(k):
    """R_k = span of r_T ^ (basis of Lambda^{k-2}), r_T = e_b e_c - e_a e_c + e_a e_b for triples T=(a,b,c) of pairs ij,jk,ik"""
    B=wedge_basis(k); bi={b:i for i,b in enumerate(B)}; rows=[]
    if k<2: return np.zeros((0,len(B)))
    for i,j,l in itertools.combinations(range(n),3):
        a,b,c=pidx[(i,j)],pidx[(j,l)],pidx[(i,l)]
        r={(b,c):1.0,(a,c):-1.0,(a,b):1.0}
        for extra in itertools.combinations(range(N),k-2):
            v=np.zeros(len(B))
            for (x,y),s in r.items():
                if x in extra or y in extra: continue
                idx=sorted((x,y)+extra); 
                # sign of sorting the sequence (x,y)+extra
                seq=list((x,y)+extra); perm=sorted(range(k),key=lambda t: seq[t]); sg=np.linalg.det(np.eye(k)[perm])
                v[bi[tuple(idx)]]+=s*sg
            if np.abs(v).max()>0: rows.append(v)
    if not rows: return np.zeros((0,len(B)))
    R=np.array(rows); 
    # orthonormal basis of the row space
    u,s,vt=np.linalg.svd(R,full_matrices=False); r=(s>1e-9).sum(); return vt[:r]
import math
perms=list(itertools.permutations(range(n)))
def cycle_type(g):
    seen=set(); ct=[]
    for i in range(n):
        if i in seen: continue
        c=0; j=i
        while j not in seen: seen.add(j); j=g[j]; c+=1
        ct.append(c)
    return tuple(sorted(ct))
reps={}
for g in perms: reps.setdefault(cycle_type(g),[]).append(g)
hF={}
for k in range(0,n):
    if k==0: hF[0]={ct:1.0 for ct in reps}; continue
    Rk=relations(k); P=Rk.T@Rk if len(Rk) else np.zeros((len(wedge_basis(k)),)*2)
    hF[k]={}
    for ct,gs in reps.items():
        M,B,bi=act_wedge(perm_on_pairs(gs[0]),k); hF[k][ct]=np.trace(M)-np.trace(P@M)
    print('H^%d(F(C,%d)): dim %d'%(k,n,round(hF[k][tuple([1]*n)])),flush=True)
# M_{0,n+1}: H^k(M) = H^k(F) - H^{k-1}(M) as virtual characters
hM={0:hF[0]}
for k in range(1,n-1): hM[k]={ct:hF[k][ct]-hM[k-1][ct] for ct in reps}
def inv(h): return sum(len(reps[ct])*h[ct] for ct in reps)/math.factorial(n)
print('Betti numbers:')
print('  F(C,%d)            :'%n,[round(hF[k][tuple([1]*n)]) for k in range(n)])
print('  F(C,%d)/S_%d        :'%(n,n),[round(inv(hF[k])) for k in range(n)])
print('  M_{0,%d}           :'%(n+1),[round(hM[k][tuple([1]*n)]) for k in range(n-1)])
print('  M_{0,%d}/S_%d       :'%(n+1,n),[round(inv(hM[k])) for k in range(n-1)])
print('  M_{0,%d}/(S_%d x Z2) :'%(n+1,n),[round(inv(hM[k])) if k%2==0 else 0 for k in range(n-1)], '(conjugation acts by (-1)^k)')
