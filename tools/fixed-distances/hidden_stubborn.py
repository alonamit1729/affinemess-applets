"""d=5 types whose primary decomposition does not finish: work with the saturated ideal Js = I : f^inf (f = product of
the inequation factors), which is computed quickly.  Singular: minimal generating set Pm of Js (mstd) and the
10x10 minors of its Jacobian (dim Js = 1 in 11 variables); msolve on Js + minors + Rabinowitsch inequation:
  empty   -> every nondegenerate point of V(Js) is a smooth point of a curve: no isolated real point;
  zerodim -> finitely many candidates: real ones are sphere-tested;
  posdim  -> inconclusive (Js not reduced along a component).  Usage: python3 hidden_stubborn.py idx [idx ...]"""
import sys, json, subprocess, time
from fractions import Fraction
import sympy as sp
idxs=[int(a) for a in sys.argv[1:]]; sys.argv=['x','5']; exec(open('certify5d.py').read().split('import os')[0])
def run_msolve(eqs,vars_,tmo=1800):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
VT=V+['t']; out=[]
for idx in idxs:
    t0=time.time(); eqs=system(classes[idx]); rab=expand(ineq); prod='*'.join(facs)
    ring='ring r=0,('+','.join(V)+'),dp;'
    fn=f'tmp/hs_{PID}.sing'
    open(fn,'w').write(ring+'\nLIB "elim.lib";\nideal I='+',\n'.join(eqs)+';\npoly f='+prod+';\nideal Js=sat(I,f); list L=mstd(Js); print("DIM"); print(dim(L[1])); ideal Pm=L[2]; print("NGEN"); print(size(Pm));\nideal M=minor(jacob(Pm),'+str(len(V))+'-dim(L[1])); M=simplify(M,2+4+8); print("NMIN"); print(size(M)); print("GENS"); print(Pm); print("MINORS"); print(M);\nquit;\n')
    try: r=subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=3600)
    except subprocess.TimeoutExpired: print(idx,'Singular timeout',flush=True); out.append({'idx':idx,'status':'singular timeout'}); continue
    lines=[l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
    try:
        dimJ=int(lines[lines.index('DIM')+1]); ngen=int(lines[lines.index('NGEN')+1]); nmin=int(lines[lines.index('NMIN')+1])
        gens=[l.rstrip(',') for l in lines[lines.index('GENS')+1:lines.index('MINORS')]]; minors=[l.rstrip(',') for l in lines[lines.index('MINORS')+1:]]
    except ValueError:
        print(idx,'Singular output incomplete:',lines[-3:],flush=True); out.append({'idx':idx,'status':'singular error'}); continue
    print(idx,'dim',dimJ,'minimal gens',ngen,'minors',nmin,'%.0fs'%(time.time()-t0),flush=True)
    try: st,md=parse_real(run_msolve(gens+minors+[rab],VT))
    except subprocess.TimeoutExpired: st,md='msolve timeout',[]
    rec={'idx':idx,'dim':dimJ,'ngens':ngen,'nminors':nmin,'singular_locus':st,'real_singular':len(md),'isolated':0,'undecided':0}
    for mid in md:
        iso=True
        for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
            cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
            sph=intstr(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),[SYM[v] for v in V])
            try: st2,md2=parse_real(run_msolve(eqs+[rab,sph],VT))
            except subprocess.TimeoutExpired: iso=None; break
            if st2!='zerodim' or md2: iso=False; break
        if iso is True: rec['isolated']+=1; rec.setdefault('isolated_points',[]).append(points_from(mid))
        if iso is None: rec['undecided']+=1
    print(rec,'%.0fs'%(time.time()-t0),flush=True); out.append(rec); json.dump(out,open('hidden_stubborn_5.json','w'))
print('done')
