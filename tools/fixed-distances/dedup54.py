"""H2 + H3 for (5,4).
For every zero-dimensional system (223 zero-dimensional types + the 12 zero-dimensional non-degenerate primes of
positive-dimensional types):
  * msolve -p 200: isolating boxes of the real solutions (certified);
  * msolve -P 2: rational parametrization; Sturm count (sympy, squarefree part) of the real roots of the eliminating
    polynomial, compared with the number of boxes (independent real-root count);
  * exact congruence classes: the group G of colour-permuting automorphisms of the colouring, times the four frame
    symmetries, acts on the roots; the image of a root box under g in G is computed by interval arithmetic
    (mpmath.iv, 300 bits) and matched to the unique isolating box it meets.  Orbits = congruence classes, exactly."""
import subprocess, json, itertools, re, sys, time
from fractions import Fraction
import sympy as sp
import mpmath
mpmath.mp.prec=300
src=open('certify54.py').read().split('results=[]')[0]; exec(src)
C=json.load(open('cert54.json'))
XVAR=['x1','x2','y2','x3','y3','x4','y4']
def boxes(out):
    s=out.replace('\n','').replace(' ','')
    tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s.rstrip(':'))
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
    if tree[0]==-1: return []
    return tree[1][1] if isinstance(tree[1],list) and len(tree[1])>1 else []
def rur_count(eqs):
    open('ms_in.txt','w').write(','.join(VARS)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f','ms_in.txt','-o','ms_rur.txt','-P','2','-t','2'],capture_output=True,text=True,timeout=600)
    s=open('ms_rur.txt').read().replace('\n','')
    m=re.search(r'\[1,\s*\[\[(\d+),\s*\[([-\d,\s]+)\]\]',s)
    deg=int(m.group(1)); cf=[int(x) for x in m.group(2).split(',')]
    assert len(cf)==deg+1
    x=sp.Symbol('x'); P=sp.Poly(list(reversed(cf)),x)
    Q=sp.sqf_part(P)
    return deg, int(Q.count_roots())
def group(col):
    G=[]
    for sg in itertools.permutations(range(5)):
        tau={}; ok=True
        for k,(i,j) in enumerate(E):
            a=col[k]; b=col[E.index(tuple(sorted((sg[i],sg[j]))))]
            if tau.setdefault(a,b)!=b: ok=False; break
        if ok and len(set(tau.values()))==4: G.append((sg,tau))
    return G
iv=mpmath.iv
def image(box,sg,tau,rot,refl):
    """box: list of (lo,hi) for VARS[:10]; returns intervals for the 7 coordinates of the framed image root."""
    B=[iv.mpf([str(lo),str(hi)]) for lo,hi in box[:10]]
    P=[(iv.mpf(0),iv.mpf(0)),(B[0],iv.mpf(0)),(B[1],B[2]),(B[3],B[4]),(B[5],B[6])]
    svals=[iv.mpf(1),B[7],B[8],B[9]]
    lam=1/svals[tau[0]]; sc=iv.sqrt(lam)
    Q=[P[sg[i]] for i in range(5)]
    Q=[((q[0]-Q[0][0])*sc,(q[1]-Q[0][1])*sc) for q in Q]
    a,b=Q[1]; r=iv.sqrt(a*a+b*b); c=a/r; s=b/r
    if rot: c,s=-c,-s
    R=[(x*c+y*s, -x*s+y*c) for x,y in Q]
    if refl: R=[(x,-y) for x,y in R]
    return [R[1][0],R[2][0],R[2][1],R[3][0],R[3][1],R[4][0],R[4][1]]
def meets(I,box):
    # compare exactly: interval endpoints (mpf, dyadic) converted to Fractions
    def fr(x): return Fraction(mpmath.nstr(x,120,strip_zeros=False)) if False else Fraction(int(mpmath.floor(x*2**290)),2**290)
    return all(Fraction(int(mpmath.floor(I[k].a*2**290)),2**290)<=box[k][1] and Fraction(int(mpmath.ceil(I[k].b*2**290)),2**290)>=box[k][0] for k in range(7))
t0=time.time(); total_orbits=0; report=[]; mism=[]
for r in C:
    systems=[]
    if r['status']=='zerodim': systems.append(('type',system(classes[r['idx']])))
    elif r['status']=='posdim':
        for pi,p in enumerate(r['primes']):
            if p['dim']==0 and not p['degenerate']: systems.append((f'prime{pi}',p['gens']+[expand(ineq)]))
    if not systems: continue
    col=classes[r['idx']]; G=group(col)
    bx=[]; sysinfo=[]
    for tag,eqs in systems:
        out=run_msolve(eqs); b=[[(lo,hi) for lo,hi in sol] for sol in boxes(out)]
        deg,sturm=rur_count(eqs)
        sysinfo.append({'system':tag,'complex_degree':deg,'real_sturm':sturm,'real_boxes':len(b)})
        if sturm!=len(b): mism.append((r['idx'],tag,'sturm',sturm,'boxes',len(b)))
        bx+=b
    # all roots of the type pooled (the group permutes the zero-dimensional primes among themselves)
    parent=list(range(len(bx)))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    amb=0
    for ai,box in enumerate(bx):
        for sg,tau in G:
            for rot in (0,1):
                for refl in (0,1):
                    I=image(box,sg,tau,rot,refl)
                    hits=[bi for bi,b2 in enumerate(bx) if meets(I,b2)]
                    if len(hits)!=1: amb+=1; continue
                    parent[find(ai)]=find(hits[0])
    roots={}
    for a in range(len(bx)): roots.setdefault(find(a),[]).append(a)
    norb_type=len(roots); nreal_type=len(bx)
    reps=[points_from([(lo+hi)/2 for lo,hi in bx[min(v)][:7]]) for v in roots.values()]
    if amb: mism.append((r['idx'],'ambiguous images',amb))
    report.append({'idx':r['idx'],'systems':sysinfo,'real_roots':nreal_type,'orbits':norb_type,'G':len(G),'ambiguous':amb,'reps':reps})
    if norb_type!=len(r['rigid']): mism.append((r['idx'],'orbits',norb_type,'numerical classes (unnormalised)',len(r['rigid'])))
    total_orbits+=norb_type
    print(r['idx'],r['status'],'real',nreal_type,'orbits',norb_type,'numerical',len(r['rigid']),'|G|',len(G),'%.0fs'%(time.time()-t0),flush=True)
print('TOTAL congruence classes (exact orbits):',total_orbits)
print('mismatches:',mism)
json.dump({'total_orbits':total_orbits,'mismatches':[list(map(str,m)) for m in mism],'report':report},open('dedup54.json','w'))
json.dump([{'idx':x['idx'],'rigid':x['reps']} for x in report],open('cert54_rigid_exact.json','w'))
