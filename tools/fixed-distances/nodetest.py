"""Certified local test at the candidate points p of the stubborn d=5 types (interval arithmetic, 300 bits).
Setting: C = V(Js) is the nondegenerate locus of the type (all components of dimension 1), given by generators g_1..g_m
in the reduced ring (linear variables eliminated).  At p:
  * five generators G_1..G_5 whose gradients are certified independent define a smooth surface germ S through p;
    C lies in S.  For another generator g, Gamma = {g|_S = 0} is a curve germ in S containing C's germ.
  * if the gradient of g|_S at p is nonzero, Gamma is a smooth real arc; if it is zero and the intrinsic Hessian
    H = Hess g - sum(lambda_i Hess G_i) restricted to T_pS is indefinite, Gamma is an ordinary node with two real
    branches (Morse lemma).  In both cases every branch of Gamma is real, and C's germ (purely one-dimensional,
    contained in Gamma) consists of such branches: p is NOT isolated in the real locus.
  * The computation certifies det(H|T_pS) < 0 by interval arithmetic on msolve's isolating box of p, which covers
    both cases at once (if the gradient is nonzero the conclusion holds regardless of H).
A definite H (det > 0) is inconclusive here (smooth arc or isolated point) and is reported.  Usage: python3 nodetest.py idx"""
import sys, json, subprocess, re, itertools, time
from fractions import Fraction
import sympy as sp, mpmath
from mpmath import iv
mpmath.mp.prec=300; iv.prec=300
idx=int(sys.argv[1]); sys.argv=['x','5']; exec(open('certify5d.py').read().split('import os')[0])
src2=open('hidden2b.py').read(); i=src2.index('def reduce_linear'); j=src2.index('out=[]; t0'); exec(src2[i:j])
def sg(script,tmo=3600):
    fn=f'tmp/nt_{PID}.sing'; open(fn,'w').write(script)
    return subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=tmo).stdout
def boxes(out):
    s=out.replace('\n','').replace(' ','').rstrip(':')
    tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s); pos=0
    def build():
        nonlocal pos
        if tokens[pos]=='[':
            pos+=1; lst=[]
            while tokens[pos]!=']':
                lst.append(build())
                if tokens[pos]==',': pos+=1
            pos+=1; return lst
        tk=tokens[pos]; pos+=1
        if '^' in tk:
            num,den=tk.split('/'); base,ex=den.split('^'); return Fraction(int(num),int(base)**int(ex))
        return Fraction(tk)
    tree=build()
    return tree[1][1] if tree[0]==0 and isinstance(tree[1],list) and len(tree[1])>1 else []
VT=V+['t']; eqs=system(classes[idx]); rab=expand(ineq); t0=time.time()
# stage-1 system and its real solutions (isolating boxes)
ring='ring r=0,('+','.join(V)+'),dp;'
out=sg(ring+'\nideal i='+',\n'.join(eqs)+';\nideal M=minor(jacob(i),10); M=simplify(M,2+4+8);\nprint("MIN");\nint k; for(k=1;k<=size(M);k++){ print(string(M[k])); }\nquit;\n')
lines=[l.strip() for l in out.splitlines() if l.strip() and not l.startswith('//')]
minors=lines[lines.index('MIN')+1:]
fi,fo=f'tmp/nt_in_{PID}.txt',f'tmp/nt_out_{PID}.txt'
open(fi,'w').write(','.join(VT)+'\n0\n'+',\n'.join(eqs+minors+[rab])+'\n')
subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=3600)
B=boxes(open(fo).read()); print(f'type {idx}: {len(B)} real candidate points, %.0fs'%(time.time()-t0),flush=True)
# reduced saturated ideal (same as hidden_subsets.py)
prod='*'.join(facs)
out=sg(ring+'\nLIB "elim.lib";\nideal I='+',\n'.join(eqs)+';\nideal Js=sat(I,'+prod+'); list L=mstd(Js); print("PM"); print(L[2]);\nquit;\n')
lines=[l.strip() for l in out.splitlines() if l.strip() and not l.startswith('//')]
gens0=[l.rstrip(',') for l in lines[lines.index('PM')+1:]]
F,subs,free=reduce_linear(gens0); vs=[SYM[v] for v in free]
G=[sp.Poly(g,*vs).clear_denoms()[1].as_expr() for g in F]
grad=[[sp.diff(g,v) for v in vs] for g in G]; hess=[[[sp.diff(g,v,w) for w in vs] for v in vs] for g in G]
fG=[sp.lambdify(vs,g,'mpmath') for g in G]; fgrad=[[sp.lambdify(vs,d,'mpmath') for d in row] for row in grad]; fhess=[[[sp.lambdify(vs,h,'mpmath') for h in row] for row in H] for H in hess]
print(f'  {len(G)} generators in {free}, %.0fs'%(time.time()-t0),flush=True)
NS=len(vs)-2   # generators cutting out a smooth surface: codim(curve) - 1
vidx=[V.index(v) for v in free]
def I(lo,hi): return iv.mpf([str(lo),str(hi)])
def det(M):
    n=len(M)
    if n==1: return M[0][0]
    if n==2: return M[0][0]*M[1][1]-M[0][1]*M[1][0]
    return sum(((-1)**c)*M[0][c]*det([row[:c]+row[c+1:] for row in M[1:]]) for c in range(n))
def solve(A,b):
    """interval Gaussian elimination (square system)"""
    n=len(A); M=[row[:]+[b[i]] for i,row in enumerate(A)]
    for c in range(n):
        piv=max(range(c,n),key=lambda r: abs(iv.mpf(M[r][c]).mid)); M[c],M[piv]=M[piv],M[c]
        for r in range(n):
            if r!=c:
                f=M[r][c]/M[c][c]; M[r]=[M[r][k]-f*M[c][k] for k in range(n+1)]
    return [M[i][n]/M[i][i] for i in range(n)]
results=[]
for bx in B:
    P=[I(lo,hi) for lo,hi in bx]; x=[P[k] for k in vidx]
    pm=[float(xx.mid) for xx in x]
    Jm=[[float(sp.N(d.subs(dict(zip(vs,pm))))) for d in row] for row in grad]
    # choose 5 generators with independent gradients (numerically), certify the 5x5 minor by intervals
    import numpy as np
    A=np.array(Jm); best=None
    for rows in itertools.combinations(range(len(G)),NS):
        sub=A[list(rows)]; s=np.linalg.svd(sub,compute_uv=False)
        if s[-1]>1e-6 and (best is None or s[-1]>best[0]): best=(s[-1],rows)
    rows=list(best[1]); JI=[[iv.mpf(fgrad[r][c](*x)) for c in range(len(vs))] for r in rows]
    # certified rank 5: some 5x5 minor interval excludes 0
    cols5=None
    for cols in itertools.combinations(range(len(vs)),NS):
        dI=iv.mpf(det([[JI[r][c] for c in cols] for r in range(NS)]))
        if dI.a>0 or dI.b<0: cols5=list(cols); break
    rec={'point':points_from([ (lo+hi)/2 for lo,hi in bx]),'rank5_certified':cols5 is not None}
    if cols5 is None: rec['verdict']='could not certify a smooth surface through p'; results.append(rec); print('  ',rec['verdict'],flush=True); continue
    freec=[c for c in range(len(vs)) if c not in cols5]   # 2 free coordinates -> kernel basis
    K=[]
    for fc in freec:
        rhs=[-JI[r][fc] for r in range(NS)]; sol=solve([[JI[r][c] for c in cols5] for r in range(NS)],rhs)
        k=[iv.mpf(0)]*len(vs)
        for c,v in zip(cols5,sol): k[c]=v
        k[fc]=iv.mpf(1); K.append(k)
    verdict='inconclusive'; detail=[]
    for jg in range(len(G)):
        if jg in rows: continue
        gj=[iv.mpf(fgrad[jg][c](*x)) for c in range(len(vs))]
        # gradient of g|S on T_pS: components along K
        gt=[sum(gj[c]*k[c] for c in range(len(vs))) for k in K]
        if any(iv.mpf(v).a>0 or iv.mpf(v).b<0 for v in gt): verdict='gradient of g|S nonzero: smooth real arc through p — NOT isolated'; detail.append(jg); break
        # Lagrange multipliers (least squares): (J J^T) lambda = J grad g
        JJt=[[sum(JI[r][c]*JI[s][c] for c in range(len(vs))) for s in range(NS)] for r in range(NS)]
        lam=solve(JJt,[sum(JI[r][c]*gj[c] for c in range(len(vs))) for r in range(NS)])
        Hj=[[iv.mpf(fhess[jg][a][b](*x)) for b in range(len(vs))] for a in range(len(vs))]
        Hi=[[[iv.mpf(fhess[r][a][b](*x)) for b in range(len(vs))] for a in range(len(vs))] for r in rows]
        H=[[Hj[a][b]-sum(lam[r]*Hi[r][a][b] for r in range(NS)) for b in range(len(vs))] for a in range(len(vs))]
        HK=[[sum(K[u][a]*H[a][b]*K[w][b] for a in range(len(vs)) for b in range(len(vs))) for w in range(2)] for u in range(2)]
        dH=iv.mpf(HK[0][0]*HK[1][1]-HK[0][1]*HK[1][0])
        if dH.b<0: verdict=f'intrinsic Hessian indefinite (det in [{mpmath.nstr(dH.a,6)}, {mpmath.nstr(dH.b,6)}]): ordinary node with two real branches — NOT isolated'; detail.append(jg); break
        if dH.a>0: detail.append(('definite',jg,mpmath.nstr(dH.a,6)))
    rec['verdict']=verdict; rec['detail']=str(detail[:3])
    print('   p =',[(round(a,3),round(b,3)) for a,b in rec['point']],'->',verdict,flush=True)
    results.append(rec)
json.dump(results,open(f'nodetest_{idx}.json','w')); print('done %.0fs'%(time.time()-t0))
