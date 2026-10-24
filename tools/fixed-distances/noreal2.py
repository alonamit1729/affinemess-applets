"""Decide real points for the components where random slices found none.
(a) certificate of emptiness: a generator that is a positive-definite quadratic (minimum > 0) -> no real points;
(b) otherwise: numerical real point on the prime (least squares from random starts, inequations checked), then
    msolve on prime + dim hyperplanes through a rational approximation of that point (+ Rabinowitsch inequation):
    a certified real nondegenerate point.  Usage: python3 noreal2.py d"""
import sys, json, random
import numpy as np, sympy as sp
from fractions import Fraction
from scipy.optimize import least_squares
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
NR=json.load(open(f'noreal_{d}.json'))
syms=[SYM[v] for v in V]
def pd_certificate(gens):
    for g in gens:
        P=sp.Poly(sp.sympify(g.replace('^','**'),locals=SYM),*syms)
        if P.total_degree()!=2: continue
        # quadratic form: write as x^T A x + b x + c ; minimum = c - b A^-1 b /4 if A positive definite
        vs=[s for s in syms if P.degree(s)>0]
        A=sp.zeros(len(vs)); b=sp.zeros(len(vs),1); c=P.as_expr().subs({s:0 for s in vs})
        for i,u in enumerate(vs):
            A[i,i]=P.as_expr().coeff(u,2)
            b[i]=P.as_expr().coeff(u,1).subs({s:0 for s in vs})
            for j,w in enumerate(vs):
                if j>i: A[i,j]=A[j,i]=sp.Rational(1,2)*P.as_expr().coeff(u,1).coeff(w,1)
        try:
            if all(A[:k,:k].det()>0 for k in range(1,len(vs)+1)):
                mn=c-(b.T*A.inv()*b)[0]/4
                if mn>0: return g,mn
        except Exception: pass
    return None
res=[]
for rec in NR:
    gens=rec['gens']; k=rec['dim']
    cert=pd_certificate(gens)
    if cert:
        print(f'd={d} type {rec["idx"]} prime {rec["prime"]}: EMPTY — generator {cert[0]} has minimum {cert[1]} > 0',flush=True)
        res.append(dict(rec,gens=None,verdict='empty',certificate=cert[0])); continue
    F=[sp.sympify(g.replace('^','**'),locals=SYM) for g in gens]
    f=sp.lambdify([syms],F,'numpy'); Jf=sp.lambdify([syms],sp.Matrix(F).jacobian(syms),'numpy')
    vals=lambda x: [1.0]+list(x[7:7+d-1])
    def ok(x):
        v=vals(x); return min(abs(a-b) for a,b in itertools.combinations(v,2))>1e-3 and min(v)>1e-3
    found=None; random.seed(rec['idx'])
    for trial in range(300):
        x0=np.random.default_rng(trial).normal(0,1.5,len(V))
        r=least_squares(lambda x: np.array(f(x),dtype=float),x0,jac=lambda x: np.array(Jf(x),dtype=float),xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=400)
        if np.max(np.abs(r.fun))<1e-10 and ok(r.x): found=r.x; break
    if found is None:
        print(f'd={d} type {rec["idx"]} prime {rec["prime"]}: no numerical real point found in 300 starts — UNDECIDED',flush=True)
        res.append(dict(rec,gens=None,verdict='undecided')); continue
    # certified: hyperplanes through a rational approximation of the point with random integer normals
    certified=False
    for attempt in range(6):
        rng=random.Random(1000*rec['idx']+attempt)
        lins=[]
        for _ in range(k):
            a=[rng.randint(-5,5) for _ in XV]
            if all(c==0 for c in a): a[0]=1
            q=[Fraction(float(found[i])).limit_denominator(10**6) for i in range(7)]
            const=sum(ai*qi for ai,qi in zip(a,q))
            lins.append('+'.join(f'({ai})*{v}' for ai,v in zip(a,XV))+f'-({const})')
        st,md=parse_real(run_msolve(gens+[expand(l) for l in lins]+[expand(ineq)],V+['t']))
        if st=='zerodim' and md:
            # pick the solution closest to the numerical point
            best=min(md,key=lambda m: sum((float(m[i])-found[i])**2 for i in range(7)))
            certified=True; break
    if certified:
        print(f'd={d} type {rec["idx"]} prime {rec["prime"]} dim {k}: REAL — certified nondegenerate real point near',[round(float(found[i]),3) for i in range(7)],flush=True)
        res.append(dict(rec,gens=None,verdict='real',point=points_from(best)))
    else:
        print(f'd={d} type {rec["idx"]} prime {rec["prime"]}: numerical point found but certification failed ({st}) — UNDECIDED',flush=True)
        res.append(dict(rec,gens=None,verdict='undecided_numeric',numpoint=points_from([Fraction(float(v)) for v in found])))
json.dump(res,open(f'noreal2_{d}.json','w'))
