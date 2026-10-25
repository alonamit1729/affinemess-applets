"""Finish the d=6 follow-up: the four primes where Singular's slocus crashes even in the reduced ring (singular locus
computed by hand as P + minors of size codim of the Jacobian), and the two primes whose singular locus is a curve
(recursive candidates as in hidden3.py, in the reduced ring).  Lifted candidates are sphere-tested in V."""
import sys, json, time, subprocess
from fractions import Fraction
import sympy as sp
d=int(sys.argv[1]); sys.argv=['x',str(d)]; exec(open('certify5d.py').read().split('import os')[0])
src=open('hidden2b.py').read(); i=src.index('def reduce_linear'); j=src.index('out=[]; t0'); exec(src[i:j])
def run_msolve(eqs,vars_,tmo=900):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
VT=V+['t']
def sg(script,tmo=1800):
    open(f'tmp/hd2c_{PID}.sing','w').write(script)
    r=subprocess.run([SG,'-q',f'tmp/hd2c_{PID}.sing'],capture_output=True,text=True,timeout=tmo)
    return [l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
def singular_locus(gens,free):
    """P + minors of size codim of the Jacobian (what slocus computes), returned as (dim, generators)"""
    ring='ring r=0,('+','.join(free)+'),dp;'
    lines=sg(ring+'\nideal P='+',\n'.join(gens)+';\nlist L=mstd(P); ideal Pm=L[2]; int c='+str(len(free))+'-dim(L[1]);\nideal S=std(Pm+minor(jacob(Pm),c));\nprint("SL"); print(dim(S)); print(S);\nquit;\n')
    k=lines.index('SL'); return int(lines[k+1]),[l.rstrip(',') for l in lines[k+2:]]
def primes_of(gens,free):
    ring='ring r=0,('+','.join(free)+'),dp;'
    lines=sg(ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(gens)+';\nlist L=minAssGTZ(i);\nint k; for(k=1;k<=size(L);k++){ print("PRIME"); print(dim(std(L[k]))); print(L[k]); }\nquit;\n')
    out=[]; cur=None
    for l in lines:
        if l=='PRIME': cur={'dim':None,'gens':[]}; out.append(cur); continue
        if cur is None: continue
        if cur['dim'] is None: cur['dim']=int(l); continue
        cur['gens'].append(l.rstrip(','))
    return out
def candidates(gens,free,depth,log):
    pts=[]
    for q in primes_of(gens,free):
        if q['dim']==0: pts.append(q['gens']); log.append(f'{"  "*depth}0-dim prime')
        else:
            sdim,sg_=singular_locus(q['gens'],free); log.append(f'{"  "*depth}{q["dim"]}-dim prime: singular locus dim {sdim}')
            if sdim==0: pts.append(sg_)
            elif sdim>0 and depth<3: pts+=candidates(sg_,free,depth+1,log)
            elif sdim>0: log.append('recursion limit')
    return pts
import os
H2b=json.load(open(f'hidden2b_{d}.json')) if os.path.exists(f'hidden2b_{d}.json') else []
jobs=[(r['idx'],r['prime'],r.get('slocus'),r.get('slocus_dim')) for r in H2b if r.get('slocus')=='error' or (r.get('slocus_dim') or 0)>0]
print('jobs',jobs,flush=True)
out=[]; cache={}; t0=time.time()
for idx,pi,status,sdim0 in jobs:
    if idx not in cache: cache[idx]=decompose(system(classes[idx]))
    p=cache[idx][pi]; F,subs,free=reduce_linear(p['gens']); vs=[SYM[v] for v in free]
    gens=[intstr(g,vs) for g in F]
    full={}
    for v in V:
        if v in free: continue
        e=subs[SYM[v]]
        for _ in range(len(subs)): e=sp.expand(e.subs(subs))
        full[v]=e
    extra=[intstr(SYM[v]-full[v],[SYM[w] for w in V]) for v in full]
    rec={'idx':idx,'prime':pi,'log':[],'candidate_systems':0,'real_candidates':0,'isolated':0,'undecided':0}
    sdim,sgens=singular_locus(gens,free); rec['slocus_dim']=sdim; rec['log'].append(f'singular locus dim {sdim}')
    systems=[]
    if sdim==0: systems=[sgens]
    elif sdim>0: systems=candidates(sgens,free,1,rec['log'])
    for S in systems:
        rec['candidate_systems']+=1
        st,md=parse_real(run_msolve(S+extra+[expand(ineq)],VT)); rec['real_candidates']+=len(md)
        for mid in md:
            iso=True
            for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
                cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
                sph=intstr(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),[SYM[v] for v in V])
                try: st2,md2=parse_real(run_msolve(system(classes[idx])+[expand(ineq),sph],VT))
                except subprocess.TimeoutExpired: iso=None; break
                if st2!='zerodim' or md2: iso=False; break
            if iso is True: rec['isolated']+=1; rec.setdefault('isolated_points',[]).append(points_from(mid))
            if iso is None: rec['undecided']+=1
    out.append(rec); print(idx,pi,'slocus dim',sdim,'candidate systems',rec['candidate_systems'],'real candidates',rec['real_candidates'],'isolated',rec['isolated'],'undecided',rec['undecided'],'%.0fs'%(time.time()-t0),flush=True)
    for l in rec['log']: print('   ',l,flush=True)
    json.dump(out,open(f'hidden2c_{d}.json','w'))
print('done')
