"""Stage 2 (unified, robust): for every type whose rank-deficient locus is non-empty, decide whether any real
nondegenerate point is isolated in the real locus.
  decompose the type (Singular, 300 s / 600 s); for each non-degenerate prime with real points:
    eliminate variables occurring linearly with constant coefficient (isomorphism), then in the reduced ring:
      singular locus = Pm + (codim)-minors of jacob(Pm), Pm a minimal generating set (mstd);
      dim -1 -> smooth: no isolated real point on this component;
      dim  0 -> real singular points (msolve, lifted, inequation imposed) -> sphere tests in V (isolation certificate);
      dim >0 -> recurse: primes of the singular locus; 0-dim ones -> points; positive-dim ones -> their singular loci.
  Everything that fails (timeouts, Singular crashes) is recorded as unresolved, never silently dropped.
Usage: python3 hidden_final.py d"""
import sys, json, time, subprocess, os
from fractions import Fraction
import sympy as sp
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
src2=open('hidden2b.py').read(); i=src2.index('def reduce_linear'); j=src2.index('out=[]; t0'); exec(src2[i:j])
def run_msolve(eqs,vars_,tmo=600):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
VT=V+['t']
def sg(script,tmo):
    fn=f'tmp/hf_{PID}.sing'; open(fn,'w').write(script)
    r=subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=tmo)
    return [l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
def decompose2(eqs):
    ring='ring r=0,('+','.join(V)+'),dp;'
    for algo,tmo in (('minAssGTZ',300),('minAssChar',600)):
        script=ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(eqs)+';\nlist L='+algo+'(i);\nint k; ideal G; int bad;\nfor(k=1;k<=size(L);k++){ G=std(L[k]); print("PRIME"); print(dim(G)); bad=0; '+' '.join(f'if(reduce({f},G)==0){{bad=1;}}' for f in facs)+' print(bad); print(L[k]); }\nquit;\n'
        try: lines=sg(script,tmo); break
        except subprocess.TimeoutExpired: lines=None
    if lines is None: return None
    primes=[]; cur=None
    for line in lines:
        if line=='PRIME': cur={'dim':None,'degenerate':None,'gens':[]}; primes.append(cur); continue
        if cur is None: continue
        if cur['dim'] is None and re.fullmatch(r'-?\d+',line): cur['dim']=int(line); continue
        if cur['degenerate'] is None and re.fullmatch(r'[01]',line): cur['degenerate']=bool(int(line)); continue
        cur['gens'].append(line.rstrip(','))
    return primes
def singular_locus(gens,free,tmo=900):
    ring='ring r=0,('+','.join(free)+'),dp;'
    lines=sg(ring+'\nideal P='+',\n'.join(gens)+';\nlist L=mstd(P); ideal Pm=L[2]; int c='+str(len(free))+'-dim(L[1]);\nideal S=std(Pm+minor(jacob(Pm),c));\nprint("SL"); print(dim(S)); print(S);\nquit;\n',tmo)
    if 'SL' not in lines: raise RuntimeError('Singular gave no result: '+' | '.join(lines[-3:]))
    k=lines.index('SL'); return int(lines[k+1]),[l.rstrip(',') for l in lines[k+2:]]
def primes_of(gens,free,tmo=900):
    ring='ring r=0,('+','.join(free)+'),dp;'
    lines=sg(ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(gens)+';\nlist L=minAssGTZ(i);\nint k; for(k=1;k<=size(L);k++){ print("PRIME"); print(dim(std(L[k]))); print(L[k]); }\nquit;\n',tmo)
    out=[]; cur=None
    for l in lines:
        if l=='PRIME': cur={'dim':None,'gens':[]}; out.append(cur); continue
        if cur is None: continue
        if cur['dim'] is None: cur['dim']=int(l); continue
        cur['gens'].append(l.rstrip(','))
    if not out: raise RuntimeError('no primes returned')
    return out
def candidates(gens,free,depth,log):
    pts=[]
    for q in primes_of(gens,free):
        if q['dim']==0: pts.append(q['gens']); log.append(f'{"  "*depth}0-dim prime of the singular locus')
        else:
            sdim,sg_=singular_locus(q['gens'],free); log.append(f'{"  "*depth}{q["dim"]}-dim prime: singular locus dim {sdim}')
            if sdim==0: pts.append(sg_)
            elif sdim>0 and depth<3: pts+=candidates(sg_,free,depth+1,log)
            elif sdim>0: raise RuntimeError('recursion limit')
    return pts
H=json.load(open(f'hidden_{d}.json')); todo=[r for r in H if r['status']!='empty']
C={r['idx']:r for r in json.load(open(f'cert5{d}.json'))}
fn_out=f'hidden_final_{d}.json'
out=json.load(open(fn_out)) if os.path.exists(fn_out) else []
done={r['idx'] for r in out}
print(d,'types:',len(todo),'already done:',len(done),flush=True)
t0=time.time()
for r in todo:
    if r['idx'] in done: continue
    idx=r['idx']; col=classes[idx]; eqs=system(col); rab=expand(ineq)
    rec={'idx':idx,'col':list(col),'primes':[],'unresolved':[],'isolated':0,'undecided':0,'log':[]}
    primes=decompose2(eqs)
    if primes is None:
        rec['unresolved'].append('decomposition did not finish'); out.append(rec); json.dump(out,open(fn_out,'w')); print(idx,'UNRESOLVED: no decomposition','%.0fs'%(time.time()-t0),flush=True); continue
    cp=C[idx]['primes']
    for pi,p in enumerate(primes):
        if p['degenerate']: continue
        if pi<len(cp) and cp[pi].get('real') is False: continue
        info={'prime':pi,'dim':p['dim']}
        try:
            F,subs,free=reduce_linear(p['gens']); vs=[SYM[v] for v in free]
            gens=[intstr(g,vs) for g in F]; info['eliminated']=len(subs)
            full={}
            for v in V:
                if v in free: continue
                e=subs[SYM[v]]
                for _ in range(len(subs)): e=sp.expand(e.subs(subs))
                full[v]=e
            extra=[intstr(SYM[v]-full[v],[SYM[w] for w in V]) for v in full]
            if not gens: sdim,sgens=-1,[]; info['note']='zero ideal in the free variables: a smooth affine component'
            else: sdim,sgens=singular_locus(gens,free)
            info['slocus_dim']=sdim
            systems=[]
            if sdim==0: systems=[sgens]
            elif sdim>0: systems=candidates(sgens,free,1,rec['log'])
            info['real_candidates']=0; info['isolated']=0; info['undecided']=0
            for SYS in systems:
                st,md=parse_real(run_msolve(SYS+extra+[rab],VT)); info['real_candidates']+=len(md)
                for mid in md:
                    iso=True
                    for eps in (Fraction(1,100),Fraction(1,1000),Fraction(1,10**6)):
                        cen=[Fraction(mid[k]).limit_denominator(10**12) for k in range(len(V))]
                        sph=intstr(sp.sympify('+'.join(f'({v}-({cen[k]}))**2' for k,v in enumerate(V))+f'-({eps*eps})',locals=SYM),[SYM[v] for v in V])
                        try: st2,md2=parse_real(run_msolve(eqs+[rab,sph],VT))
                        except subprocess.TimeoutExpired: iso=None; break
                        if st2!='zerodim' or md2: iso=False; break
                    if iso is True: info['isolated']+=1; info.setdefault('isolated_points',[]).append(points_from(mid))
                    if iso is None: info['undecided']+=1
            rec['isolated']+=info['isolated']; rec['undecided']+=info['undecided']
        except Exception as ex:
            info['error']=str(ex)[:200]; rec['unresolved'].append(f'prime {pi}: {str(ex)[:120]}')
        rec['primes'].append(info)
    out.append(rec); json.dump(out,open(fn_out,'w'))
    print(idx,[(q['dim'],q.get('slocus_dim',q.get('error','?')[:30]),q.get('real_candidates'),q.get('isolated')) for q in rec['primes']],'unresolved' if rec['unresolved'] else '','%.0fs'%(time.time()-t0),flush=True)
print('done')
