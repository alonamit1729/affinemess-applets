"""Strategy-B test: are configurations near the collision ends of one-parameter families extendable by a point that
keeps the distance palette?  For each (5,d) family type (saturated ideal, x1 = 1): trace a short way from every
degenerate root (two points coinciding) into the family and count one-point extensions at several shape distances
from the end; compare with extension counts at random interior points and at the rigid sets of the cell."""
import sys, json, math, itertools, random
d=int(sys.argv[1]); sys.argv=['strat.py',str(d)]
src=open('strat.py').read().split("if os.path.exists(f'sat5{d}.json'):")[0]
exec(src)
def palette(P,tol=1e-7):
    vals=[]
    for i,j in itertools.combinations(range(len(P)),2):
        x=math.dist(P[i],P[j])
        if not any(abs(x-v)<tol*max(1,x) for v in vals): vals.append(x)
    return vals
def extensions(P,tol=1e-6):
    """points q with all distances to P in the palette (candidates: circle-circle intersections)"""
    vals=palette(P); mx=max(vals); found=[]
    def ok(q):
        return all(any(abs(math.dist(q,p)-v)<tol*mx for v in vals) for p in P) and all(math.dist(q,p)>tol*mx for p in P)
    for (i,j) in itertools.combinations(range(len(P)),2):
        a,b=P[i],P[j]; D=math.dist(a,b)
        for r1 in vals:
            for r2 in vals:
                if D>r1+r2+tol or D<abs(r1-r2)-tol: continue
                x=(r1*r1-r2*r2+D*D)/(2*D); h2=r1*r1-x*x
                if h2<-tol: continue
                h=math.sqrt(max(h2,0)); ux,uy=(b[0]-a[0])/D,(b[1]-a[1])/D
                for s in (1,-1):
                    q=(a[0]+x*ux-s*h*uy,a[1]+x*uy+s*h*ux)
                    if ok(q) and not any(math.dist(q,f)<1e-6*mx for f in found): found.append(q)
    return found
SAT=json.load(open(f'sat5{d}.json')); R=json.load(open(f'cert5{d}.json'))
random.seed(1)
res=[]; t0=time.time()
for k,v in sorted(SAT.items(),key=lambda kv:int(kv[0])):
    idx=int(k)
    if v['dim']!=1: continue
    gens=v['gens']; C=Curve(gens)
    roots=[]
    for eq,lab in coinc:
        st,md=parse_real(run_msolve(gens+['x1-1',expand(eq)],V))
        for x in md:
            if not any(np.abs(np.array(q)-np.array(x)).max()<1e-8 for q in roots): roots.append(x)
    degs=[r for r in roots if coincident(pts_of(r))]
    for r in degs:
        v0=np.array(r); T=C.tangents(v0)
        for u in [t for tt in T for t in (tt,-tt)]:
            w=None
            for off in (1e-3,1e-2,3e-2):
                w1,ok=C.project(v0+off*u/np.linalg.norm(u))
                if ok and boundary_hit(w1) is None: w=w1; break
            if w is None: continue
            vend,path,hit=trace(C,w,u,roots,start=roots.index(r),maxsteps=4000)
            fp=fp_path(path)[:,0]   # min distance / diameter along the path
            rec={'type':idx,'end':'collision','samples':[]}
            for target in (0.02,0.05,0.1,0.2):
                ks=np.where(fp>=target)[0]
                if len(ks)==0: continue
                P=[tuple(p) for p in pts_of(path[ks[0]])]
                rec['samples'].append({'fp0':float(fp[ks[0]]),'ext':len(extensions(P))})
            if rec['samples']: res.append(rec)
    # interior points: random slices
    for _ in range(2):
        lin='+'.join(f'({random.randint(-5,5)})*{x}' for x in XV)+f'-({random.randint(-3,3)})'
        st,md=parse_real(run_msolve(gens+['x1-1',expand(lin)],V))
        for x in md[:2]:
            P=[tuple(p) for p in pts_of(np.array(x))]
            if coincident(P) or margin(np.array(x))<1e-6: continue
            res.append({'type':idx,'end':'interior','samples':[{'fp0':float(fp_path([x])[0,0]),'ext':len(extensions(P))}]})
    print(idx,'done %.0fs'%(time.time()-t0),flush=True)
# rigid sets of the cell
rig=[]
for r in R:
    for P in r['rigid'][:3]:
        rig.append(len(extensions([tuple(p) for p in P])))
json.dump({'families':res,'rigid_ext':rig},open(f'g7fin/extend_ends_{d}.json','w'))
import collections
def summ(lst): return f'n={len(lst)} extendable={sum(1 for x in lst if x>0)} mean={np.mean(lst):.2f}' if lst else 'none'
for target in (0.02,0.05,0.1,0.2):
    L=[s_['ext'] for rec in res if rec['end']=='collision' for s_ in rec['samples'] if abs(s_['fp0']-target)<0.5*target+1e-9]
    print(f'near collision ends, min/diam ~{target}:',summ(L))
print('interior points:',summ([s_['ext'] for rec in res if rec['end']=='interior' for s_ in rec['samples']]))
print('rigid sets:',summ(rig))
