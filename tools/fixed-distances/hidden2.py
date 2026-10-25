"""Per-component follow-up for the types whose rank-deficient locus has positive dimension.
Decompose the type (minimal primes); for each non-degenerate prime with real points: Singular slocus(prime)
 (singular locus of the reduced component); dim -1: smooth, no isolated real points possible; dim 0: real singular
 points (msolve, inequation imposed) -> sphere tests in V; dim >= 1: reported.  Usage: python3 hidden2.py d [timeout]"""
import sys, json, time, subprocess
from fractions import Fraction
d=int(sys.argv[1]); TMO=int(sys.argv[2]) if len(sys.argv)>2 else 1200
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
VT=V+['t']
H=json.load(open(f'hidden_{d}.json')); todo=[r for r in H if r['status']!='empty']
C={r['idx']:r for r in json.load(open(f'cert5{d}.json'))}
import os
out=json.load(open(f'hidden2_{d}.json')) if os.path.exists(f'hidden2_{d}.json') else []
done={r['idx'] for r in out}
print(d,'types to follow up:',len(todo),flush=True)
t0=time.time()
for r in todo:
    if r['idx'] in done: continue
    col=classes[r['idx']]; eqs=system(col); rab=expand(ineq)
    primes=decompose(eqs)
    rec={'idx':r['idx'],'col':list(col),'primes':[]}
    if primes is None:
        rec['note']='decomposition did not finish'; out.append(rec); json.dump(out,open(f'hidden2_{d}.json','w')); print(r['idx'],'no decomposition',flush=True); continue
    cp=C[r['idx']]['primes']
    for pi,p in enumerate(primes):
        if p['degenerate']: continue
        realflag=cp[pi].get('real') if pi<len(cp) else None
        if realflag is False: continue
        ring='ring r=0,('+','.join(V)+'),dp;'
        script=ring+'\nLIB "sing.lib";\nideal P='+',\n'.join(p['gens'])+';\nideal S=std(slocus(P));\nprint("SL"); print(dim(S)); print(S);\nquit;\n'
        open(f'tmp/hd2_{PID}.sing','w').write(script)
        info={'prime':pi,'dim':p['dim']}
        try:
            rr=subprocess.run([SG,'-q',f'tmp/hd2_{PID}.sing'],capture_output=True,text=True,timeout=TMO)
            lines=[l.strip() for l in rr.stdout.splitlines() if l.strip() and not l.startswith('//')]
            if 'SL' not in lines:
                print('   Singular output without SL for type',r['idx'],'prime',pi,':',rr.stdout[-600:],rr.stderr[-300:],flush=True)
                info['slocus']='error'; rec['primes'].append(info); continue
            k=lines.index('SL'); sdim=int(lines[k+1]); sgens=[l.rstrip(',') for l in lines[k+2:]]
        except subprocess.TimeoutExpired:
            info['slocus']='timeout'; rec['primes'].append(info); continue
        info['slocus_dim']=sdim
        if sdim==0:
            st,md=parse_real(run_msolve(sgens+[rab],VT)); info['real_singular']=len(md); info['isolated']=0
            for mid in md:
                iso=True
                for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
                    cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]; sph=str(sp.Poly(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),*[SYM[v] for v in V]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
                    try: st2,md2=parse_real(run_msolve(eqs+[rab,sph],VT))
                    except subprocess.TimeoutExpired: iso=None; break
                    if st2!='zerodim' or md2: iso=False; break
                if iso is True: info['isolated']+=1; info.setdefault('isolated_points',[]).append(points_from(mid))
                if iso is None: info['undecided']=info.get('undecided',0)+1
        rec['primes'].append(info)
    out.append(rec)
    print(r['idx'],[(q['dim'],q.get('slocus_dim',q.get('slocus')),q.get('real_singular'),q.get('isolated')) for q in rec['primes']],'%.0fs'%(time.time()-t0),flush=True)
    json.dump(out,open(f'hidden2_{d}.json','w'))
print('done')
