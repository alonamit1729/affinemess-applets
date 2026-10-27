"""Is every known n-point <=6-distance pattern rigid?  For each pattern class in g7fin/patterns_{n}new.json: msolve on the
squared-distance equations (scale fixed by the first colour, distinct nonzero values via a Rabinowitsch variable):
'empty' (no realization with distinct values — should not happen for a pattern read off a real set), 'zerodim' (finitely
many realizations: rigid) or 'posdim'.  Usage: python3 g7fin/pat_rigid.py n"""
import sys, json, itertools, subprocess, re, time
from os import getpid
import sympy as sp
from fractions import Fraction
n=int(sys.argv[1]); PID=str(getpid())
MS='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
pats=json.load(open(f'g7fin/patterns_{n}new.json'))
def run_msolve(eqs,vars_,tmo=1800):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    try: subprocess.run([MS,'-f',fi,'-o',fo,'-p','64'],capture_output=True,text=True,timeout=tmo)
    except subprocess.TimeoutExpired: return None
    return open(fo).read()
def status(out):
    s=out.replace('\n','').replace(' ','')
    if s.startswith('[-1]'): return 'empty'
    if s.startswith('[1,'): return 'posdim'
    return 'zerodim'
out=[]; t0=time.time()
fn=f'g7fin/pat_rigid_{n}new.json'
import os
done=json.load(open(fn)) if os.path.exists(fn) else []
for pi,pat in enumerate(pats):
    if pi<len(done): continue
    k=pat['k']; col={tuple(map(int,key.split(','))):c for key,c in pat['col'].items()}
    # relabel colours so that the colour of edge (0,1) is 0
    c01=col[(0,1)]; relab={c01:0}; 
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
    subprocess.run([MS,'-f',fi,'-o',fo,'-g','2'],capture_output=True,text=True,timeout=600)
    gb=open(fo).read(); body=gb[gb.index('['):] if '[' in gb else ''
    if body.strip() in ('[1]','[1,]') or body.strip().startswith('[1]'): st='empty'
    else:
        pure=set()
        for term in body.strip('[]\n ').split(','):
            term=term.strip().lstrip('-')
            if not term: continue
            lm=term.split('+')[0].split('-')[0].strip()
            facs=[f for f in lm.split('*') if not f.strip().lstrip('-').isdigit()]
            if len(facs)!=1: continue
            v=facs[0].split('^')[0].strip()
            if v in V: pure.add(v)
        st='zerodim_modp' if all(v in pure for v in V) else 'posdim_modp'
    if st=='zerodim_modp':
        SG='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/Singular'
        script='ring r=0,('+','.join(V)+'),dp;\nideal i='+',\n'.join(seqs)+';\nprint(dim(std(i))); quit;\n'
        open(f'tmp/sg_in_{PID}.sing','w').write(script)
        try:
            r=subprocess.run([SG,'-q',f'tmp/sg_in_{PID}.sing'],capture_output=True,text=True,timeout=900)
            dimq=int(r.stdout.strip().splitlines()[-1]); st='zerodim' if dimq==0 else ('empty' if dimq<0 else 'posdim')
        except Exception as ex: st='zerodim_modp (Q timeout)'
    done.append({'name':pat['name'],'k':k,'status':st}); json.dump(done,open(fn,'w'))
    print(pi,pat['name'],k,st,'%.0fs'%(time.time()-t0),flush=True)
import collections; print(collections.Counter(r['status'] for r in done))
