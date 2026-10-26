"""Degeneration poset of the colouring types of 5-point d-distance sets (d = 4..7): for every realizable
positive-dimensional type T (saturated ideal from sat5{d}.json, x1 = 1) and every coincidence equation
(two distance values equal, or one distance value zero), decompose the resulting ideal over Q into minimal primes
(Singular), keep the non-degenerate real components, and identify each component's colouring type and dimension from
a sample point (random slice -> msolve).  Output: degen_{d}.json: list of {type, eq, components:[{dim, col_type, d_new,
coincident, sample}]}.  Usage: python3 degen.py d"""
import sys, json, os, subprocess, re, time, random, itertools, math
sys.argv=[sys.argv[0]]+sys.argv[1:]
src=open('certify5d.py').read().split('import os')[0]
exec(src)   # classes, system, facs, V, XV, SV, SG, MS, PID, expand, run_msolve, parse_real, decompose, points_from
import numpy as np
d=int(sys.argv[1])
SAT=json.load(open(f'sat5{d}.json'))
Sval=['1']+SV
coinc=[(f'{Sval[a]}-{Sval[b]}',['eq',a,b]) for a,b in itertools.combinations(range(d),2)]+[(f'{Sval[k]}',['zero',k]) for k in range(1,d)]
# canonical colouring type of a sample configuration (fewer distances): partition of the 10 edges by distance value
E=list(itertools.combinations(range(5),2))
def coltype(P,tol=1e-6):
    ds=[math.dist(P[i],P[j]) for i,j in E]; m=max(ds)
    vals=[]; col=[]
    for x in ds:
        for k,v in enumerate(vals):
            if abs(x-v)<tol*m: col.append(k); break
        else: vals.append(x); col.append(len(vals)-1)
    return canon(col), len(vals)
def decompose_sat(eqs,tmo=150):
    ring='ring r=0,('+','.join(V)+'),dp;'
    script=ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(eqs)+';\nlist L=minAssGTZ(i);\nint k; ideal G;\nfor(k=1;k<=size(L);k++){ G=std(L[k]); print("PRIME"); print(dim(G)); print(L[k]); }\nquit;\n'
    open(f'tmp/sg_in_{PID}.sing','w').write(script)
    try: r=subprocess.run([SG,'-q',f'tmp/sg_in_{PID}.sing'],capture_output=True,text=True,timeout=tmo)
    except subprocess.TimeoutExpired: return None
    primes=[]; cur=None
    for line in r.stdout.splitlines():
        line=line.strip()
        if line=='PRIME': cur={'dim':None,'gens':[]}; primes.append(cur); continue
        if cur is None or line.startswith('//') or not line: continue
        if cur['dim'] is None and re.fullmatch(r'-?\d+',line): cur['dim']=int(line); continue
        g=line.split('=',1)[1].strip() if '=' in line else line
        cur['gens'].append(g.rstrip(','))
    return primes
outfn=f'degen_{d}.json'
out=json.load(open(outfn)) if os.path.exists(outfn) else []
done={(o['type'],tuple(o['eq'])) for o in out}
t0=time.time()
for k,v in sorted(SAT.items(),key=lambda kv:int(kv[0])):
    idx=int(k)
    if v.get('timeout'): continue
    gens=v['gens']
    for eq,lab in coinc:
        if (idx,tuple(lab)) in done: continue
        rec={'type':idx,'dim':v['dim'],'eq':lab,'components':[]}
        primes=decompose_sat(gens+[expand(eq)])
        if primes is None: rec['timeout']=True; out.append(rec); json.dump(out,open(outfn,'w')); print(idx,lab,'TIMEOUT',flush=True); continue
        random.seed(idx*1000+hash(tuple(lab))%1000)
        for p in primes:
            comp={'dim':p['dim'],'real':False}
            # real sample point: random slices (dim of them) -> msolve
            for attempt in range(6):
                lins=['+'.join(f'({random.randint(-5,5)})*{x}' for x in XV)+f'-({random.randint(-3,3)})' for _ in range(p['dim'])]
                st,md=parse_real(run_msolve(p['gens']+[expand(l) for l in lins],V))
                if st=='zerodim' and md:
                    P=points_from(md[0]); comp['real']=True; comp['sample']=[list(map(float,q)) for q in P]
                    comp['coincident']=any(math.dist(P[i],P[j])<1e-7 for i,j in E)
                    if not comp['coincident']:
                        ct,nv=coltype(P); comp['col_type']=list(ct); comp['d_new']=nv
                    break
                if st=='empty' and p['dim']==0: break
            rec['components'].append(comp)
        out.append(rec); json.dump(out,open(outfn,'w'))
        print(idx,lab,[(c['dim'],c.get('d_new'),c.get('coincident')) for c in rec['components'] if c['real']],'%.0fs'%(time.time()-t0),flush=True)
print('done')
