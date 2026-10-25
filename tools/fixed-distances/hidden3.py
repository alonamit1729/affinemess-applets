"""Recursive follow-up for components whose singular locus has positive dimension.
A real point p isolated in V_R lying on a prime P must be an isolated real point of Sing(P)_R, hence (recursively) a
real point of a zero-dimensional piece of the chain P ⊃ Sing(P) ⊃ Sing(Sing(P)) ⊃ ... or a singular point of a
positive-dimensional piece.  We collect all such candidate real points and sphere-test them in V.
Usage: python3 hidden3.py d"""
import sys, json, time, subprocess
from fractions import Fraction
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
def run_msolve(eqs,vars_,tmo=600):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
VT=V+['t']; ring='ring r=0,('+','.join(V)+'),dp;'
def sg(script,tmo=1800):
    open(f'tmp/hd3_{PID}.sing','w').write(script)
    r=subprocess.run([SG,'-q',script and f'tmp/hd3_{PID}.sing'],capture_output=True,text=True,timeout=tmo)
    return [l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
def primes_of(gens):
    lines=sg(ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(gens)+';\nlist L=minAssGTZ(i);\nint k; for(k=1;k<=size(L);k++){ print("PRIME"); print(dim(std(L[k]))); print(L[k]); }\nquit;\n')
    out=[]; cur=None
    for l in lines:
        if l=='PRIME': cur={'dim':None,'gens':[]}; out.append(cur); continue
        if cur is None: continue
        if cur['dim'] is None: cur['dim']=int(l); continue
        cur['gens'].append(l.rstrip(','))
    return out
def slocus_of(gens):
    lines=sg(ring+'\nLIB "sing.lib";\nideal P='+',\n'.join(gens)+';\nideal S=std(slocus(P));\nprint("SL"); print(dim(S)); print(S);\nquit;\n')
    k=lines.index('SL'); return int(lines[k+1]),[l.rstrip(',') for l in lines[k+2:]]
def candidates(gens,depth,log):
    """real points that could be isolated in the real locus of V(gens): zero-dim pieces and singular points, recursively"""
    pts=[]
    for q in primes_of(gens):
        if q['dim']==0:
            st,md=parse_real(run_msolve(q['gens']+[expand(ineq)],VT)); pts+=md; log.append(f'{"  "*depth}0-dim prime: {len(md)} real points')
        else:
            sdim,sg_=slocus_of(q['gens']); log.append(f'{"  "*depth}{q["dim"]}-dim prime: singular locus dim {sdim}')
            if sdim==0:
                st,md=parse_real(run_msolve(sg_+[expand(ineq)],VT)); pts+=md; log.append(f'{"  "*depth}  {len(md)} real singular points')
            elif sdim>0 and depth<3: pts+=candidates(sg_,depth+1,log)
            elif sdim>0: log.append(f'{"  "*depth}  recursion limit')
    return pts
H2=json.load(open(f'hidden2_{d}.json'))
todo=[r for r in H2 if any(q.get('slocus_dim',0)>0 or q.get('slocus')=='timeout' for q in r['primes'])]
print(d,'types with positive-dimensional singular loci:',[r['idx'] for r in todo],flush=True)
out=[]; t0=time.time()
for r in todo:
    col=classes[r['idx']]; eqs=system(col); rab=expand(ineq)
    primes=decompose(eqs); C={x['idx']:x for x in json.load(open(f'cert5{d}.json'))}; cp=C[r['idx']]['primes']
    rec={'idx':r['idx'],'log':[],'candidates':0,'isolated':0,'undecided':0}
    for pi,p in enumerate(primes):
        if p['degenerate'] or (pi<len(cp) and cp[pi].get('real') is False): continue
        sdim,sgens=slocus_of(p['gens'])
        if sdim<=0: continue
        rec['log'].append(f'prime {pi} (dim {p["dim"]}): singular locus of dim {sdim}')
        pts=candidates(sgens,1,rec['log'])
        for mid in pts:
            rec['candidates']+=1; iso=True
            for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
                cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
                sph=str(sp.Poly(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),*[SYM[v] for v in V]).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
                try: st2,md2=parse_real(run_msolve(eqs+[rab,sph],VT))
                except subprocess.TimeoutExpired: iso=None; break
                if st2!='zerodim' or md2: iso=False; break
            if iso is True: rec['isolated']+=1; rec.setdefault('isolated_points',[]).append(points_from(mid))
            if iso is None: rec['undecided']+=1
    out.append(rec); print(r['idx'],'candidates',rec['candidates'],'isolated',rec['isolated'],'undecided',rec['undecided'],'%.0fs'%(time.time()-t0),flush=True)
    for l in rec['log']: print('   ',l)
    json.dump(out,open(f'hidden3_{d}.json','w'))
print('done')
