"""H5: no hidden rigid configurations for d >= 5.
A real nondegenerate configuration of a type that is isolated in the real locus must be a point where the Jacobian of
the ten distance equations has rank < 10 (at a full-rank point the real locus is a manifold of dimension nvars-10 = d-4 >= 1).
So for each type: msolve on {equations} + {all 10x10 minors of the Jacobian} + Rabinowitsch inequation.
  empty      -> no rank-deficient nondegenerate point at all: no hidden rigid configuration;
  zerodim    -> finitely many such points; for each real one, sphere tests (radii 1e-2, 1e-3, 1e-6) decide isolation;
  posdim     -> rank drops along a component (reported; handled separately).
Usage: python3 hidden.py d"""
import sys, json, time
import sympy as sp
from fractions import Fraction
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
syms=[SYM[v] for v in V]
def run_msolve(eqs,vars_,tmo=300):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
out=[]; t0=time.time()
import os
done={x['idx'] for x in json.load(open(f'hidden_{d}.json'))} if os.path.exists(f'hidden_{d}.json') else set()
if done: out=json.load(open(f'hidden_{d}.json'))
for ci,col in enumerate(classes):
    if ci in done: continue
    eqs=system(col)
    ring='ring r=0,('+','.join(V)+'),dp;'
    script=ring+'\nideal i='+',\n'.join(eqs)+';\nideal M=minor(jacob(i),10);\nM=simplify(M,2+4+8);\nprint("MIN");\nint k; for(k=1;k<=size(M);k++){ print(string(M[k])); }\nquit;\n'
    open(f'tmp/hd_{PID}.sing','w').write(script)
    rr=subprocess.run([SG,'-q',f'tmp/hd_{PID}.sing'],capture_output=True,text=True,timeout=900)
    lines=[l.strip() for l in rr.stdout.splitlines() if l.strip() and not l.startswith('//')]
    minors=lines[lines.index('MIN')+1:]
    try:
        o=run_msolve(eqs+minors+[expand(ineq)],V+['t']); st,md=parse_real(o)
    except Exception as ex:
        st,md='error:'+str(ex)[:60],[]
    rec={'idx':ci,'col':list(col),'status':st,'n_real_singular':len(md),'isolated':[]}
    for mid in md:
        rec['isolated'].append({'point':points_from(mid),'isolated':'deferred to hidden2 (singular points of the reduced components)'}); continue
        iso=True
        for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
            cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]; sph=str(sp.Poly(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),*[SYM[v] for v in V]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
            try: st2,md2=parse_real(run_msolve(eqs+[expand(ineq),sph],V+['t']))
            except subprocess.TimeoutExpired: iso=None; break
            if st2!='zerodim' or md2: iso=False; break
        rec['isolated'].append({'point':points_from(mid),'isolated':iso})
    out.append(rec)
    flag='' if st=='empty' else ('  <-- '+st)
    print(ci,'/',len(classes),st,'real singular',len(md),'%.0fs'%(time.time()-t0),flag,flush=True)
    json.dump(out,open(f'hidden_{d}.json','w'))
print('done')
