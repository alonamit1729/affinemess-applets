"""d=6 type 242, primes 1 and 3: Singular crashes computing the singular locus even in the reduced ring.
Here: minimal generating set from Singular (mstd), Jacobian minors with sympy, and msolve decides whether the
singular locus is empty / finite (then real points -> sphere tests) / positive-dimensional."""
import sys, json, subprocess, itertools
from fractions import Fraction
import sympy as sp
d=6; sys.argv=['x','6']; exec(open('certify5d.py').read().split('import os')[0])
src2=open('hidden2b.py').read(); i=src2.index('def reduce_linear'); j=src2.index('out=[]; t0'); exec(src2[i:j])
def run_msolve(eqs,vars_,tmo=1800):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
VT=V+['t']; idx=242; eqs=system(classes[idx]); rab=expand(ineq)
primes=decompose(eqs); out=[]
for pi in (1,3):
    p=primes[pi]; F,subs,free=reduce_linear(p['gens']); vs=[SYM[v] for v in free]
    gens=[intstr(g,vs) for g in F]
    ring='ring r=0,('+','.join(free)+'),dp;'
    fn=f'tmp/h242_{PID}.sing'; open(fn,'w').write(ring+'\nideal P='+',\n'.join(gens)+';\nlist L=mstd(P); print("DIM"); print(dim(L[1])); ideal Pm=L[2]; print("PM"); print(Pm); ideal M=minor(jacob(Pm),'+str(len(free))+'-dim(L[1])); M=simplify(M,2+4+8); print("MIN"); print(M);\nquit;\n')
    r=subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=1800)
    lines=[l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
    dimP=int(lines[lines.index('DIM')+1]); Pm=[l.rstrip(',') for l in lines[lines.index('PM')+1:lines.index('MIN')]]; minors=[l.rstrip(',') for l in lines[lines.index('MIN')+1:]]
    full={}
    for v in V:
        if v in free: continue
        e=subs[SYM[v]]
        for _ in range(len(subs)): e=sp.expand(e.subs(subs))
        full[v]=e
    extra=[intstr(SYM[v]-full[v],[SYM[w] for w in V]) for v in full]
    print(f'prime {pi}: dim {dimP} in {len(free)} free vars, {len(Pm)} minimal gens, {len(minors)} minors',flush=True)
    st,md=parse_real(run_msolve(Pm+sorted(minors)+extra+[rab],VT))
    rec={'prime':pi,'dim':dimP,'free':free,'nminors':len(minors),'singular_locus':st,'real_singular':len(md),'isolated':0,'undecided':0}
    for mid in md:
        iso=True
        for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
            cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
            sph=intstr(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),[SYM[v] for v in V])
            try: st2,md2=parse_real(run_msolve(eqs+[rab,sph],VT))
            except subprocess.TimeoutExpired: iso=None; break
            if st2!='zerodim' or md2: iso=False; break
        if iso is True: rec['isolated']+=1
        if iso is None: rec['undecided']+=1
    print(rec,flush=True); out.append(rec)
json.dump(out,open('hidden_242.json','w')); print('done')
