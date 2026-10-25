"""Smoothness certificate via generator subsets.  Sing(P) is contained in {rank Jac(S) < codim} for every subset S of
the generators with |S| = codim, so adding the codim-minors of several random subsets to P gives a superset of the
singular locus; if msolve finds it empty, P is smooth and its real points are never isolated.
Usage: python3 hidden_subsets.py d idx prime|sat [nsubsets]"""
import sys, json, subprocess, random, time
import sympy as sp
d=int(sys.argv[1]); idx=int(sys.argv[2]); which=sys.argv[3]; nsub=int(sys.argv[4]) if len(sys.argv)>4 else 30
sys.argv=['x',str(d)]; exec(open('certify5d.py').read().split('import os')[0])
src2=open('hidden2b.py').read(); i=src2.index('def reduce_linear'); j=src2.index('out=[]; t0'); exec(src2[i:j])
def run_msolve(eqs,vars_,tmo=3600):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    return open(fo).read()
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
def sg(script,tmo=1800):
    fn=f'tmp/hsub_{PID}.sing'; open(fn,'w').write(script)
    r=subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=tmo)
    return [l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
VT=V+['t']; eqs=system(classes[idx]); rab=expand(ineq); t0=time.time()
if which=='sat':
    prod='*'.join(facs); ring='ring r=0,('+','.join(V)+'),dp;'
    lines=sg(ring+'\nLIB "elim.lib";\nideal I='+',\n'.join(eqs)+';\nideal Js=sat(I,'+prod+'); list L=mstd(Js); print("DIM"); print(dim(L[1])); print("PM"); print(L[2]);\nquit;\n')
    dimP=int(lines[lines.index('DIM')+1]); gens0=[l.rstrip(',') for l in lines[lines.index('PM')+1:]]
    F,subs,free=reduce_linear(gens0); vs=[SYM[v] for v in free]; gens=[intstr(g,vs) for g in F]
    full={}
    for v in V:
        if v in free: continue
        e=subs[SYM[v]]
        for _ in range(len(subs)): e=sp.expand(e.subs(subs))
        full[v]=e
    extra=[intstr(SYM[v]-full[v],[SYM[w] for w in V]) for v in full]
    print(f'  saturated ideal: {len(gens0)} minimal generators; after eliminating {len(subs)} linear variables: {len(gens)} generators in {len(free)} variables',flush=True)
else:
    primes=decompose(eqs); p=primes[int(which)]
    F,subs,free=reduce_linear(p['gens']); vs=[SYM[v] for v in free]; gens=[intstr(g,vs) for g in F]; dimP=p['dim']
    full={}
    for v in V:
        if v in free: continue
        e=subs[SYM[v]]
        for _ in range(len(subs)): e=sp.expand(e.subs(subs))
        full[v]=e
    extra=[intstr(SYM[v]-full[v],[SYM[w] for w in V]) for v in full]
codim=len(free)-dimP
print(f'type {idx} ({which}): dim {dimP}, {len(gens)} generators in {len(free)} variables, codim {codim}',flush=True)
random.seed(idx); minors=set()
ring='ring r=0,('+','.join(free)+'),dp;'
for k in range(nsub):
    S=random.sample(range(len(gens)),codim)
    lines=sg(ring+'\nideal S='+',\n'.join(gens[s] for s in S)+';\nideal M=minor(jacob(S),'+str(codim)+'); M=simplify(M,2+4+8); print("MIN"); print(M);\nquit;\n')
    if 'MIN' not in lines: print('  subset',k,'Singular failed',flush=True); continue
    minors.update(l.rstrip(',') for l in lines[lines.index('MIN')+1:])
print(f'  {len(minors)} distinct minors from {nsub} subsets, %.0fs'%(time.time()-t0),flush=True)
# Groebner basis of P + minors in the reduced ring (Singular), then msolve on the (small) basis with the inequations
try:
    lines=sg(ring+'\nideal P='+',\n'.join(gens)+';\nideal M='+',\n'.join(sorted(minors))+';\nideal G=std(P+M); print("GDIM"); print(dim(G)); print("GB"); print(G);\nquit;\n',tmo=3600)
    gdim=int(lines[lines.index('GDIM')+1]); gb=[l.rstrip(',') for l in lines[lines.index('GB')+1:]]
    print(f'  Groebner basis of P + minors: dim {gdim}, {len(gb)} elements, %.0fs'%(time.time()-t0),flush=True)
    if gdim<0: st,md='empty',[]
    else:
        try: st,md=parse_real(run_msolve(gb+extra+[rab],VT))
        except subprocess.TimeoutExpired: st,md='msolve timeout',[]
except subprocess.TimeoutExpired: st,md='singular timeout',[]
print(f'  superset of the singular locus (with inequations): {st}, real points: {len(md)}, %.0fs'%(time.time()-t0),flush=True)
json.dump({'idx':idx,'which':which,'dim':dimP,'ngens':len(gens),'nminors':len(minors),'status':st,'real_points':[points_from(m) for m in md]},open(f'hidden_subsets_{d}_{idx}_{which}.json','w'))
