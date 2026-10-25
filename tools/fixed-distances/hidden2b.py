"""Retry of the singular-locus computation for primes on which Singular's slocus crashed: eliminate the variables that
occur linearly with constant coefficient (an isomorphism of the component), compute slocus in the smaller ring, and lift
the real singular points back for the sphere test in V.  Usage: python3 hidden2b.py d"""
import sys, json, time, subprocess
from fractions import Fraction
import sympy as sp
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
def run_msolve(eqs,vars_,tmo=600):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
VT=V+['t']
H2=json.load(open(f'hidden2_{d}.json'))
jobs=[(r['idx'],q['prime']) for r in H2 for q in r['primes'] if q.get('slocus') in ('error','timeout')]
print(d,'primes to retry:',jobs,flush=True)
def reduce_linear(gens):
    F=[sp.sympify(g.replace('^','**'),locals=SYM) for g in gens]; subs={}
    changed=True
    while changed:
        changed=False
        for i,f in enumerate(F):
            if f==0: continue
            P=sp.Poly(f,*[SYM[v] for v in V if SYM[v] not in subs])
            for v in P.gens:
                if P.degree(v)==1:
                    c=sp.Poly(f,v).all_coeffs()[0]
                    if c.is_number and c!=0:
                        sol=sp.solve(f,v)[0]; subs[v]=sol
                        F=[sp.expand(g.subs(v,sol)) for g in F]; changed=True; break
            if changed: break
    F=[g for g in F if g!=0]
    free=[v for v in V if SYM[v] not in subs]
    return F,subs,free
out=[]; t0=time.time()
C={x['idx']:x for x in json.load(open(f'cert5{d}.json'))}
cache={}
for idx,pi in jobs:
    if idx not in cache: cache[idx]=decompose(system(classes[idx]))
    p=cache[idx][pi]; F,subs,free=reduce_linear(p['gens'])
    gens=[str(sp.Poly(g,*[SYM[v] for v in free]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','') for g in F]   # Singular parses a^2/4 as a^(2/4)
    ring='ring r=0,('+','.join(free)+'),dp;'
    script=ring+'\nLIB "sing.lib";\nideal P='+',\n'.join(gens)+';\nprint("DIMP"); print(dim(std(P)));\nideal S=std(slocus(P));\nprint("SL"); print(dim(S)); print(S);\nquit;\n'
    open(f'tmp/hd2b_{PID}.sing','w').write(script)
    rec={'idx':idx,'prime':pi,'dim':p['dim'],'eliminated':len(subs),'free':free}
    try:
        rr=subprocess.run([SG,'-q',f'tmp/hd2b_{PID}.sing'],capture_output=True,text=True,timeout=1800)
        lines=[l.strip() for l in rr.stdout.splitlines() if l.strip() and not l.startswith('//')]
        if 'SL' not in lines: rec['slocus']='error'; rec['singular_out']=rr.stdout[-300:]; out.append(rec); print(idx,pi,'ERROR again',flush=True); continue
        k=lines.index('SL'); sdim=int(lines[k+1]); sgens=[l.rstrip(',') for l in lines[k+2:]]; rec['dim_reduced']=int(lines[lines.index('DIMP')+1])
    except subprocess.TimeoutExpired:
        rec['slocus']='timeout'; out.append(rec); print(idx,pi,'TIMEOUT',flush=True); continue
    rec['slocus_dim']=sdim
    if sdim==0:
        # real singular points in the reduced variables; lift and sphere-test in V (with the inequation)
        # express every eliminated variable as a polynomial in the free ones: add equations var - expr
        # fully resolve the chained substitutions, then clear denominators for msolve
        full={}
        for v in V:
            if v in free: continue
            e=subs[SYM[v]]
            for _ in range(len(subs)): e=sp.expand(e.subs(subs))
            full[v]=e
        extra=[str(sp.Poly(SYM[v]-full[v],*[SYM[w] for w in V]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','') for v in full]
        assert all(set(map(str,full[v].free_symbols))<=set(free) for v in full)
        st,md=parse_real(run_msolve(sgens+extra+[expand(ineq)],VT)); rec['real_singular']=len(md); rec['isolated']=0
        for mid in md:
            iso=True
            for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
                cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
                sph=str(sp.Poly(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),*[SYM[v] for v in V]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
                try: st2,md2=parse_real(run_msolve(system(classes[idx])+[expand(ineq),sph],VT))
                except subprocess.TimeoutExpired: iso=None; break
                if st2!='zerodim' or md2: iso=False; break
            if iso is True: rec['isolated']+=1; rec.setdefault('isolated_points',[]).append(points_from(mid))
            if iso is None: rec['undecided']=rec.get('undecided',0)+1
    out.append(rec); print(idx,pi,'eliminated',len(subs),'vars; slocus dim',sdim,'real singular',rec.get('real_singular'),'isolated',rec.get('isolated'),'%.0fs'%(time.time()-t0),flush=True)
    json.dump(out,open(f'hidden2b_{d}.json','w'))
print('done')
