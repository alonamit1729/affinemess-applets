"""Stratification of the one-parameter families of 5-point d-distance sets (d = 4 from the primes in cert54.json).
For each real non-degenerate 1-dimensional prime P (with x1 = 1, i.e. one representative of the mirror pair):
  * exact boundary: real roots of P + (one coincidence equation), one equation per pair of distance values and per
    distance value = 0 (coincident points); identified up to similarity with the rigid 5-point sets of fewer distances
    or recorded as degenerate (coincident points);
  * numerical tracing of the real curve from each boundary root and from random interior points: arcs between boundary
    roots, closed loops, and ends at infinity (the base edge becomes infinitely short relative to the rest);
  * everything is then mapped to shape space (unlabelled, up to similarity) to get the incidence graph.
Usage: python3 strat.py d     (d = 4 reads cert54.json; writes strat_d.json)"""
import json, sys, math, itertools, subprocess, re, time, random, collections
from os import getpid
import numpy as np, sympy as sp
PID=str(getpid())
MS='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
d=int(sys.argv[1]); n=5
E=list(itertools.combinations(range(n),2))
XV=['x1','x2','y2','x3','y3','x4','y4']; SV=[f's{k}' for k in range(1,d)]; V=XV+SV
SYM={v:sp.Symbol(v) for v in V}
Sval=['1']+SV
# ---------- exact boundary roots ----------
def run_msolve(eqs,vars_):
    fi,fo=f'tmp/ms_in_{PID}.txt',f'tmp/ms_out_{PID}.txt'
    open(fi,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f',fi,'-o',fo,'-p','200'],capture_output=True,text=True,timeout=900)
    return open(fo).read()
from fractions import Fraction
def parse_real(out):
    s=out.replace('\n','').replace(' ','')
    if s.startswith('[-1]'): return 'empty',[]
    if s.startswith('[1,'): return 'posdim',[]
    tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s.rstrip(':'))
    if not tokens or tokens[0]!='[': return 'error',[]
    def build(i):
        if tokens[i]=='[':
            lst=[]; i+=1
            while tokens[i]!=']':
                v,i=build(i); lst.append(v)
                if tokens[i]==',': i+=1
            return lst,i+1
        tk=tokens[i]
        if '^' in tk:
            num,den=tk.split('/'); base,ex=den.split('^'); return Fraction(int(num),int(base)**int(ex)),i+1
        return Fraction(tk),i+1
    tree,_=build(0)
    sols=tree[1][1] if isinstance(tree[1],list) and len(tree[1])>1 else []
    return 'zerodim',[[float((lo+hi)/2) for lo,hi in sol] for sol in sols]
def expand(expr): return str(sp.expand(sp.sympify(expr,locals=SYM))).replace('**','^').replace(' ','')
coinc=[(f'{Sval[a]}-{Sval[b]}',('eq',a,b)) for a,b in itertools.combinations(range(d),2)]+[(f'{Sval[k]}',('zero',k)) for k in range(1,d)]
# ---------- shapes ----------
def pts_of(v): return np.array([[0,0],[v[0],0],[v[1],v[2]],[v[3],v[4]],[v[5],v[6]]],float)
def dmat(P): return np.array([[math.dist(a,b) for b in P] for a in P])
PERMS=list(itertools.permutations(range(n)))
def canon(P,tol=1e-6):
    """canonical form of a point set up to similarity (incl. reflections) and relabelling: lexicographically minimal
    normalised distance matrix over all relabellings; rounding makes it hashable."""
    D=dmat(P); m=D[D>1e-9].min() if (D>1e-9).any() else 1.0; D=D/m
    best=None
    for p in PERMS:
        t=tuple(round(D[p[i],p[j]],6) for i in range(n) for j in range(i+1,n))
        if best is None or t<best: best=t
    return best
def nd(P):
    D=dmat(P); vals=sorted(D[np.triu_indices(n,1)]); m=vals[-1]
    k=1
    for a,b in zip(vals,vals[1:]):
        if b-a>1e-7*m: k+=1
    return k
def coincident(P):
    return any(math.dist(P[i],P[j])<1e-7*max(1.0,dmat(P).max()) for i,j in E)
# ---------- tracing ----------
class Curve:
    def __init__(self,gens):
        self.syms=[SYM[v] for v in V]
        self.F=[]
        for g in list(gens)+['x1-1']:
            f=sp.sympify(g.replace('^','**'),locals=SYM); m=max(abs(c) for c in sp.Poly(f,*self.syms).coeffs())
            self.F.append(sp.expand(f/m))
        self.f=sp.lambdify(self.syms,self.F,'numpy')
        J=sp.Matrix(self.F).jacobian(self.syms); self.j=sp.lambdify(self.syms,J,'numpy')
    def res(self,v): return np.array(self.f(*v),float)
    def jac(self,v): return np.array(self.j(*v),float)
    def project(self,v,it=30):
        for _ in range(it):
            r=self.res(v)
            if np.abs(r).max()<1e-13: return v,True
            J=self.jac(v); dv=np.linalg.lstsq(J,-r,rcond=1e-7)[0]; v=v+dv
            if np.abs(dv).max()>1e3: return v,False
        return v,np.abs(self.res(v)).max()<1e-9
    def tangents(self,v):
        J=self.jac(v); u,s,vt=np.linalg.svd(J)
        null=[vt[i] for i in range(len(V)) if (s[i] if i<len(s) else 0)<1e-7]
        return null
def svals(v): return [1.0]+list(v[7:7+d-1])
def margin(v):
    s=svals(v); return min(min(abs(a-b) for a,b in itertools.combinations(s,2)), min(s[1:]))
def boundary_hit(v):
    s=svals(v); sc=max(1.0,max(abs(x) for x in v))
    for (a,b) in itertools.combinations(range(d),2):
        if abs(s[a]-s[b])<1e-6*sc: return ('eq',a,b)
    for k in range(1,d):
        if s[k]<1e-6*sc: return ('zero',k)
    if max(abs(x) for x in v[:7])>1e4: return ('inf',)
    return None
def trace(C,v0,dir0,roots,start=None,home=None,maxsteps=300000):
    """follow the curve from v0 in direction dir0 until it arrives at one of the exact boundary roots (proximity),
    returns to home (loop), or escapes to infinity.  Returns (endpoint, path samples, hit)."""
    v=v0.copy(); t=dir0/np.linalg.norm(dir0); path=[v.copy()]; steps=0; left_start=False
    RT=np.array(roots) if len(roots) else np.zeros((0,len(V)))
    while steps<maxsteps:
        steps+=1
        T=C.tangents(v)
        if not T: return v,path,('singular',)
        tt=max(T,key=lambda u: abs(u@t))
        if tt@t<0: tt=-tt
        if len(T)>1 and abs(tt@t)<0.5: return v,path,('singular',)
        t=tt
        sc=max(1.0,max(abs(x) for x in v))
        dr=np.abs(RT-v).max(axis=1) if len(RT) else np.array([np.inf])
        h=1e-3*sc if dr.min()>1e-2*sc else 2e-4*sc
        w,ok=C.project(v+h*t)
        if ok and (w-v)@t<=0: ok=False
        if not ok:
            for mult in (0.25,2,4,8,16,32):
                w,ok=C.project(v+mult*h*t)
                if ok and (w-v)@t>0.3*mult*h: break
                ok=False
            if not ok: return v,path,('lost',)
        v=w; path.append(v.copy())
        sc=max(1.0,max(abs(x) for x in v))
        if sc>300: return v,path,('inf',)
        dr=np.abs(RT-v).max(axis=1) if len(RT) else np.array([np.inf])
        if start is not None and not left_start:
            if dr[start]>3e-3*sc: left_start=True
            else: continue
        k=int(np.argmin(dr))
        if dr[k]<6e-4*sc: return RT[k],path,('root',k)
        if home is not None and steps>50 and np.abs(v-home).max()<6e-4*sc: return v,path,('loop',)
    return v,path,('maxsteps',)
def ratio_path(path):
    A=np.array(path)[:,:7]; P=np.stack([np.zeros((len(A),2)),np.stack([A[:,0],np.zeros(len(A))],1),A[:,1:3],A[:,3:5],A[:,5:7]],1)
    D=np.linalg.norm(P[:,:,None,:]-P[:,None,:,:],axis=3); iu=np.triu_indices(5,1); Dv=D[:,iu[0],iu[1]]
    return Dv.max(axis=1)/Dv.min(axis=1)
def ratio(v):
    D=dmat(pts_of(v)); vals=D[np.triu_indices(n,1)]; return vals.max()/vals.min()
# ---------- main ----------
R=json.load(open(f'cert5{d}.json'))
import os
outfn=f'strat_{d}.json'
out=json.load(open(outfn)) if os.path.exists(outfn) else []
done={(o['type'],o['prime']) for o in out}
t0=time.time()
def fp_path(path):
    """permutation- and similarity-invariant fingerprint of every point of a path: sorted distances / diameter"""
    A=np.array(path)[:,:7]; P=np.stack([np.zeros((len(A),2)),np.stack([A[:,0],np.zeros(len(A))],1),A[:,1:3],A[:,3:5],A[:,5:7]],1)
    D=np.linalg.norm(P[:,:,None,:]-P[:,None,:,:],axis=3); iu=np.triu_indices(5,1); Dv=np.sort(D[:,iu[0],iu[1]],axis=1)
    return Dv/Dv[:,-1:]
def truncate_fp(fp):
    """cut the path where it becomes nearly degenerate (min distance < 2% of the diameter), so that the two traversals
    of an arc, and its symmetric copies, are compared over the same shape range"""
    bad=np.where(fp[:,0]<0.05)[0]
    if len(bad)==0: return fp
    if bad[0]==0:   # starts degenerate: cut the initial part
        good=np.where(fp[:,0]>=0.05)[0]
        if len(good)==0: return fp
        fp=fp[good[0]:]; bad=np.where(fp[:,0]<0.05)[0]
        if len(bad)==0: return fp
    return fp[:bad[0]] if bad[0]>2 else fp
def same_arc(fpA,fpB,tol=4e-3):
    """every one of 40 samples of A lies on B and vice versa"""
    for X,Y in ((fpA,fpB),(fpB,fpA)):
        idx=np.linspace(0,len(X)-1,40).astype(int)
        for k in idx:
            if np.abs(Y-X[k]).max(axis=1).min()>tol: return False
    return True
def process(idx,col,pi,gens):
        rec={'type':idx,'col':col,'prime':pi,'roots':[],'arcs':[],'loops':0}
        for eq,lab in coinc:
            st,md=parse_real(run_msolve(gens+['x1-1',expand(eq)],V))
            if st=='posdim': rec.setdefault('posdim_eqs',[]).append(lab); continue
            for v in md:
                P=pts_of(v)
                if not any(np.abs(np.array(q['v'])-np.array(v)).max()<1e-8 for q in rec['roots']):
                    rec['roots'].append({'v':v,'lab':lab,'ndist':nd(P),'coincident':coincident(P),'canon':canon(P)})
        C=Curve(gens)
        geo=[]   # geometric arcs: {'ends':[end,end], 'fp': fingerprint path (subsampled), 'mult', 'ratio_max', 'canon_mid', 'len'}
        def endrep(i):
            q=rec['roots'][i]
            return {'end':'degenerate','root':i,'canon':q['canon']} if q['coincident'] else {'root':i,'canon':q['canon']}
        def same_end(a,b):
            if a.get('end')!=b.get('end'): return False
            if 'canon' not in a or 'canon' not in b: return True
            return max(abs(x-y) for x,y in zip(a['canon'],b['canon']))<2e-5
        def add_arc(e1,e2,path,hit):
            fp=truncate_fp(fp_path(path))
            for g in geo:
                if ((same_end(g['ends'][0],e1) and same_end(g['ends'][1],e2)) or (same_end(g['ends'][0],e2) and same_end(g['ends'][1],e1))) and same_arc(fp,g['fpd']):
                    g['mult']+=1; return
            rs=ratio_path(path); kmax=int(np.argmax(rs))
            sub=fp[np.linspace(0,len(fp)-1,min(len(fp),400)).astype(int)]
            geo.append({'ends':[e1,e2],'fpd':fp,'fp':sub.tolist(),'mult':1,'ratio_max':float(rs[kmax]),'ratio_min':float(rs.min()),
                        'canon_rmax':canon(pts_of(path[kmax])),'canon_mid':canon(pts_of(path[len(path)//2])),'len':len(path),'hit':list(hit)})
        for i,q in enumerate(rec['roots']):
            v0=np.array(q['v']); T=C.tangents(v0)
            dirs=[]
            for u in T: dirs+= [u,-u]
            if len(T)>1:
                for _ in range(6):
                    c=np.random.randn(len(T)); dirs.append(sum(ci*u for ci,u in zip(c,T)))
            starts=[]
            for u in dirs:
                w=None
                for off in (1e-4,1e-3,1e-2,3e-2):
                    w1,ok=C.project(v0+off*u/np.linalg.norm(u))
                    if ok and boundary_hit(w1) is None: w=w1; break
                if w is None: continue
                if any(np.abs(s_-w).max()<1e-6 for s_ in starts): continue
                starts.append(w)
                vend,path,hit=trace(C,w,u,[q['v'] for q in rec['roots']],start=i)
                e1=endrep(i)
                if hit[0]=='root': e2=endrep(hit[1])
                elif hit[0]=='inf': e2={'end':'inf'}
                elif hit[0] in ('lost','singular','maxsteps') and max(abs(x) for x in vend)>100: e2={'end':'inf'}; hit=('inf',)
                else: e2={'end':hit[0]}
                add_arc(e1,e2,[v0]+path,hit)
        random.seed(idx*100+pi)
        for _ in range(3):
            lin='+'.join(f'({random.randint(-5,5)})*{v}' for v in XV)+f'-({random.randint(-3,3)})'
            st,md=parse_real(run_msolve(gens+['x1-1',expand(lin)],V))
            for v in md:
                w=np.array(v)
                if margin(w)<=1e-9: continue
                T=C.tangents(w)
                if len(T)!=1: continue
                vend,path,hit=trace(C,w,T[0],[q['v'] for q in rec['roots']],home=w)
                if hit[0]=='loop': add_arc({'end':'loop'},{'end':'loop'},path,hit); rec['loops']+=1
                elif hit[0]=='root': pass
                else: rec.setdefault('slice_odd',[]).append(list(hit))
        rec['arcs']=[{k:v for k,v in g.items() if k!='fpd'} for g in geo]
        out.append(rec); json.dump(out,open(outfn,'w'))
        kinds=collections.Counter(tuple(sorted(e.get('end','vertex') for e in g['ends'])) for g in geo)
        print(f"type {idx}: roots {len(rec['roots'])} geometric arcs {len(geo)} {dict(kinds)} mults {collections.Counter(g['mult'] for g in geo)} {time.time()-t0:.0f}s",flush=True)
if os.path.exists(f'sat5{d}.json'):
    SAT=json.load(open(f'sat5{d}.json'))
    cols={r['idx']:r['col'] for r in R}
    for k,v in sorted(SAT.items(),key=lambda kv:int(kv[0])):
        idx=int(k)
        if v.get('timeout') or v['dim']!=1 or (idx,-1) in done: continue
        process(idx,cols[idx],-1,v['gens'])
else:
  for r in R:
    if r.get('status')!='posdim': continue
    for pi,p in enumerate(r['primes']):
        if p['dim']!=1 or p['degenerate'] or not p.get('real_slices',p.get('real')) or (r['idx'],pi) in done: continue
        st,_=parse_real(run_msolve(p['gens']+['x1+1'],V))
        if st!='empty': continue
        process(r['idx'],r['col'],pi,p['gens'])
print('done',len(out),'curves')
