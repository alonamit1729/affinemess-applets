import itertools, sys
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
n=int(sys.argv[1])
pairs=list(itertools.combinations(range(n),2))
# variables p_1..p_{n-2} (p_0=0, p_{n-1}=1). difference d_ij = p_j - p_i as affine form: (coeffs vector length n-2, const)
def diff(i,j):
    v=[Fr(0)]*(n-2); c=Fr(0)
    for idx,sgn in ((j,1),(i,-1)):
        if idx==0: continue
        if idx==n-1: c+=sgn
        else: v[idx-1]+=sgn
    return v,c
D={(i,j):diff(i,j) for i,j in pairs}
hyper=[]  # equations d_ij - d_kl = 0 as (vector, const)
for a,b in itertools.combinations(pairs,2):
    v=[x-y for x,y in zip(D[a][0],D[b][0])]; c=D[a][1]-D[b][1]
    if all(x==0 for x in v): continue
    hyper.append((tuple(v),c,a,b))
def rref(rows):
    rows=[list(r) for r in rows]; m=len(rows); ncol=n-1
    piv=[]; r=0
    for col in range(n-2):
        pr=None
        for i in range(r,m):
            if rows[i][col]!=0: pr=i;break
        if pr is None: continue
        rows[r],rows[pr]=rows[pr],rows[r]
        pv=rows[r][col]; rows[r]=[x/pv for x in rows[r]]
        for i in range(m):
            if i!=r and rows[i][col]!=0:
                f=rows[i][col]; rows[i]=[x-f*y for x,y in zip(rows[i],rows[r])]
        piv.append(col); r+=1
        if r==m: break
    for row in rows:
        if all(x==0 for x in row[:-1]) and row[-1]!=0: return None
    rows=[tuple(x) for x in rows[:r]]
    rows=[row for row in rows if any(x!=0 for x in row[:-1])]
    return tuple(sorted(rows))
def feasible(key):
    # maximize t s.t. p_{i+1}-p_i >= t, equations hold; variables p_1..p_{n-2}, t
    m=n-2
    A_eq=[]; b_eq=[]
    for row in key:
        A_eq.append([float(x) for x in row[:-1]]+[0.0]); b_eq.append(float(-row[-1]))
    A_ub=[]; b_ub=[]
    for i in range(n-1):
        rowv=[0.0]*(m+1)
        # p_{i+1}-p_i >= t  ->  -(p_{i+1}-p_i)+t <= 0
        c=0.0
        for idx,sgn in ((i+1,-1),(i,1)):
            if idx==0: continue
            if idx==n-1: c+=sgn
            else: rowv[idx-1]+=sgn
        rowv[m]=1.0; A_ub.append(rowv); b_ub.append(-c)
    res=linprog([0.0]*m+[-1.0],A_ub=A_ub,b_ub=b_ub,A_eq=A_eq if A_eq else None,b_eq=b_eq if b_eq else None,bounds=[(None,None)]*m+[(None,1)],method='highs')
    return res.status==0 and -res.fun>1e-9
def in_span(key,vec,c):
    # is equation vec.p + c = 0 implied by key?
    k2=rref(list(key)+[tuple(list(vec)+[c])])
    return k2 is not None and len(k2)==len(key)
def count_classes(key):
    # generic number of distinct differences on the flat
    cls=list(range(len(pairs)));
    def find(x):
        while cls[x]!=x: x=cls[x]
        return x
    for v,c,a,b in hyper:
        if in_span(key,v,c):
            ra,rb=find(pairs.index(a)),find(pairs.index(b))
            if ra!=rb: cls[ra]=rb
    # also pairs whose differences are identically equal (never happens for distinct pairs on a line except... skip)
    return len(set(find(i) for i in range(len(pairs))))
levels=[{(): None}]
flats={}  # key -> dim
flats[()]=n-2
frontier=[()]
for d in range(1,n-1):
    new=set()
    for key in frontier:
        for v,c,a,b in hyper:
            k2=rref(list(key)+[tuple(list(v)+[c])])
            if k2 is None or len(k2)!=len(key)+1: continue
            if k2 in new or k2 in flats: continue
            if feasible(k2): new.add(k2)
    for k2 in new: flats[k2]=n-2-d
    frontier=list(new)
    print('codim',d,'flats',len(new),file=sys.stderr)
info={key:(dim,count_classes(key)) for key,dim in flats.items()}
# maximal flats per k: flat F with count k is a component unless contained in a flat G with count k and bigger dim
def contains(G,F):  # G contains F  <=> every equation of G is implied by F
    return all(in_span(F,row[:-1],row[-1]) for row in G)
byk={}
for key,(dim,k) in info.items(): byk.setdefault(k,[]).append((dim,key))
print('n =',n,': expected dim for k is',n-2,'- (',len(pairs),'- k )')
for k in sorted(byk):
    comps=[]
    for dim,key in byk[k]:
        if any(dim2>dim and contains(key2,key) for dim2,key2 in byk[k]): continue
        comps.append((dim,key))
    # mod out reversal p_i -> 1 - p_{n-1-i}
    def reverse_key(key):
        rows=[]
        for row in key:
            v=list(row[:-1]); c=row[-1]
            # substitute p_i = 1 - q_{n-1-i}; for i in 1..n-2: p_i -> 1 - q_{n-1-i}
            nv=[Fr(0)]*(n-2); nc=c
            for i,coef in enumerate(v):
                idx=i+1; j=n-1-idx  # q_j index j in 1..n-2
                nc+=coef; nv[j-1]-=coef
            rows.append(tuple(nv+[nc]))
        return rref(rows)
    seen=set(); comps2=[]
    for dim,key in comps:
        rk=reverse_key(key)
        if key in seen or rk in seen: continue
        seen.add(key); comps2.append((dim,key))
    comps=comps2
    dims=sorted(set(d for d,_ in comps))
    expected=n-2-(len(pairs)-k)
    cnt={d:sum(1 for dd,_ in comps if dd==d) for d in dims}
    # sample points for 0-dim components
    ex=[]
    for dim,key in comps:
        if dim==0:
            # solve
            m=n-2; A=np.array([[float(x) for x in row[:-1]] for row in key]); b=np.array([float(-row[-1]) for row in key])
            p=np.linalg.lstsq(A,b,rcond=None)[0]; pts=[0.0]+list(p)+[1.0]
            # rescale to small integers if possible
            from fractions import Fraction
            fr=[Fraction(x).limit_denominator(1000) for x in pts]; import math
            L=1
            for f in fr: L=L*f.denominator//math.gcd(L,f.denominator)
            ex.append([int(f*L) for f in fr])
    print(f'  k={k}: expected dim {expected}; components by dim: {cnt}' + (f'; rigid examples: {ex}' if ex else ''))
    if n<=5:
        import sympy as sp
        for dim,key in comps:
            if dim>=1:
                syms=sp.symbols('p1:%d'%(n-1))
                eqs=[sum(sp.Rational(row[i])*syms[i] for i in range(n-2))+sp.Rational(row[-1]) for row in key]
                sol=sp.solve(eqs,syms,dict=True)[0] if eqs else {}
                free=[s_ for s_ in syms if s_ not in sol]
                expr=[0]+[sp.simplify(sol.get(s_,s_)) for s_ in syms]+[1]
                print('      family:',expr)
