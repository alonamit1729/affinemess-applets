"""Sampled degeneration data (fallback for the exact decomposition): random slices of (saturated ideal + coincidence equation) -> real points -> their colouring types.  Degeneration poset of the colouring types of 5-point d-distance sets (d = 4..7): for every realizable
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
outfn=f'degen_s_{d}.json'
out=json.load(open(outfn)) if os.path.exists(outfn) else []
done={(o['type'],tuple(o['eq'])) for o in out}
t0=time.time()
def run_msolve_t(eqs,vars_,tmo=90):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    try: subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=tmo)
    except subprocess.TimeoutExpired: return None
    return open(fo).read()
for k,v in sorted(SAT.items(),key=lambda kv:int(kv[0])):
    idx=int(k)
    if v.get('timeout'): continue
    gens=v['gens']; dimT=v['dim']
    for eq,lab in coinc:
        if (idx,tuple(lab)) in done: continue
        rec={'type':idx,'dim':dimT,'eq':lab,'samples':[],'timeouts':0}
        random.seed(idx*1000+hash(tuple(lab))%1000)
        seen=set()
        for attempt in range(5):
            lins=['+'.join(f'({random.randint(-5,5)})*{x}' for x in XV)+f'-({random.randint(-3,3)})' for _ in range(dimT-1)]
            o=run_msolve_t(gens+[expand(eq)]+[expand(l) for l in lins],V)
            if o is None: rec['timeouts']+=1; continue
            st,md=parse_real(o)
            if st=='posdim': rec['posdim']=True; break
            for mid in md:
                P=points_from(mid)
                if any(math.dist(P[i],P[j])<1e-7 for i,j in E): key=('degenerate',)
                else:
                    ct,nv=coltype(P); key=(nv,)+tuple(ct)
                if key not in seen:
                    seen.add(key); rec['samples'].append({'coincident':key[0]=='degenerate','d_new':None if key[0]=='degenerate' else key[0],'col_type':None if key[0]=='degenerate' else list(key[1:]),'sample':[list(map(float,q)) for q in P]})
        out.append(rec); json.dump(out,open(outfn,'w'))
        print(idx,lab,sorted({(s_['d_new']) for s_ in rec['samples']},key=str),'timeouts',rec['timeouts'],'%.0fs'%(time.time()-t0),flush=True)
print('done')
