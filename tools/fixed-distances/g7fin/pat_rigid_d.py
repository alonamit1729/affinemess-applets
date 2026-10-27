"""Copy of g7fin/pat_rigid.py reading g7fin/patterns_{n}d{k}.json and writing g7fin/pat_rigid_{n}d{k}.json.
Same method: squared-distance equations with p0=(0,0), p1=(x1,0), colour of edge (0,1) has value 1, free values s_c,
Rabinowitsch variable t for 'all values distinct and nonzero'; msolve Groebner basis mod p screens (zero-dimensional iff a pure
power of every variable occurs among the leading monomials), Singular dim(std(I)) over Q certifies.  Difference from
pat_rigid.py: Singular is run in every case (also when the mod-p screen says positive-dimensional) and the dimension over Q
is recorded ('dimQ'), plus the mod-p screen result and timings.  Usage: python3 g7fin/pat_rigid_d.py n k"""
import sys, json, itertools, subprocess, re, time, os
from os import getpid
import sympy as sp
n=int(sys.argv[1]); K=int(sys.argv[2]); PID=str(getpid())
MS='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
SG='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/Singular'
pats=json.load(open(f'g7fin/patterns_{n}d{K}.json'))
out=[]; t0=time.time()
fn=f'g7fin/pat_rigid_{n}d{K}.json'
done=json.load(open(fn)) if os.path.exists(fn) else []
for pi,pat in enumerate(pats):
    if pi<len(done): continue
    t1=time.time()
    k=pat['k']; col={tuple(map(int,key.split(','))):c for key,c in pat['col'].items()}
    # relabel colours so that the colour of edge (0,1) is 0
    c01=col[(0,1)]; relab={c01:0}
    for c in range(k):
        if c!=c01: relab[c]=len(relab)
    col={e:relab[c] for e,c in col.items()}
    XV=['x1']+[v for i in range(2,n) for v in (f'x{i}',f'y{i}')]; SV=[f's{c}' for c in range(1,k)]; V=XV+SV+['t']
    SYM={v:sp.Symbol(v) for v in V}
    coords={0:(0,0),1:(SYM['x1'],0)}; coords.update({i:(SYM[f'x{i}'],SYM[f'y{i}']) for i in range(2,n)})
    S=[1]+[SYM[v] for v in SV]
    eqs=[]
    for (i,j),c in col.items():
        (a,b),(cc,d)=coords[i],coords[j]
        eqs.append(sp.expand((a-cc)**2+(b-d)**2-S[c]))
    facs=[S[a]-S[b] for a,b in itertools.combinations(range(k),2)]+[S[c] for c in range(1,k)]
    ineq=SYM['t']*sp.Mul(*facs)-1
    eqs.append(sp.expand(ineq))
    seqs=[str(e).replace('**','^').replace(' ','') for e in eqs]
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(V)+'\n0\n'+',\n'.join(seqs)+'\n')
    try:
        subprocess.run([MS,'-f',fi,'-o',fo,'-g','2'],capture_output=True,text=True,timeout=1200)
        gb=open(fo).read(); body=gb[gb.index('['):] if '[' in gb else ''
        if body.strip() in ('[1]','[1,]') or body.strip().startswith('[1]'): modp='empty_modp'
        else:
            pure=set()
            for term in body.strip('[]\n ').split(','):
                term=term.strip().lstrip('-')
                if not term: continue
                lm=term.split('+')[0].split('-')[0].strip()
                fs=[f for f in lm.split('*') if not f.strip().lstrip('-').isdigit()]
                if len(fs)!=1: continue
                v=fs[0].split('^')[0].strip()
                if v in V: pure.add(v)
            modp='zerodim_modp' if all(v in pure for v in V) else 'posdim_modp'
    except subprocess.TimeoutExpired: modp='timeout_modp'
    tms=time.time()-t1
    # Singular over Q, always
    script='ring r=0,('+','.join(V)+'),dp;\nideal i='+',\n'.join(seqs)+';\nprint(dim(std(i))); quit;\n'
    open(f'tmp/sg_in_{PID}.sing','w').write(script)
    dimq=None
    try:
        r=subprocess.run([SG,'-q',f'tmp/sg_in_{PID}.sing'],capture_output=True,text=True,timeout=1800)
        dimq=int(r.stdout.strip().splitlines()[-1])
        st='zerodim' if dimq==0 else ('empty' if dimq<0 else 'posdim')
    except Exception as ex: st=modp+' (Q timeout)'
    done.append({'name':pat['name'],'k':k,'status':st,'modp':modp,'dimQ':dimq,'nvars':len(V),
                 't_msolve':round(tms,1),'t_singular':round(time.time()-t1-tms,1)})
    json.dump(done,open(fn,'w'))
    print(pi,pat['name'],k,st,modp,dimq,'%.0fs'%(time.time()-t0),flush=True)
import collections; print(collections.Counter(r['status'] for r in done))
