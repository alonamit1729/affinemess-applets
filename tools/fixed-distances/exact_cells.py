"""Exact lists for the cells classified by the complete pattern enumeration + real root isolation (g7count/):
exact_cells.json = {"n,d": [{"pts":..., "label":...}]} with (10,5), (10,6), (11,6), (12,6) from sets10_by_d.json /
sets1113_by_d.json and (16,7) from the two Erdős–Fishburn sets.  Labels: regular-polygon subsets (with/without centre),
triangular/square lattice pieces, otherwise the distance spectrum."""
import json, math, itertools, cmath
S10=json.load(open('g7count/sets10_by_d.json')); S11=json.load(open('g7count/sets1113_by_d.json'))
def circ(p,q,r):
    (ax,ay),(bx,by),(cx,cy)=p,q,r
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-12: return None
    return (((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d,((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d)
def polygon_label(P):
    n=len(P); best=(0,None)
    for i,j,l in itertools.combinations(range(n),3):
        c=circ(P[i],P[j],P[l])
        if c is None: continue
        R=math.dist(c,P[i]); on=[t for t in range(n) if abs(math.dist(c,P[t])-R)<1e-7]
        if len(on)>best[0]: best=(len(on),(c,R,on))
    n_on,(c,R,on)=best
    if n_on<n-1: return None
    angs=sorted(math.degrees(math.atan2(P[t][1]-c[1],P[t][0]-c[0]))%360 for t in on)
    gaps=[(angs[(i+1)%n_on]-angs[i])%360 for i in range(n_on)]
    g=min(gaps); N=360/g
    if abs(N-round(N))>1e-6: return None
    N=int(round(N))
    if any(abs(x/g-round(x/g))>1e-6 for x in gaps): return None
    off=[t for t in range(n) if t not in on]
    centre=len(off)==1 and math.dist(P[off[0]],c)<1e-7
    if off and not centre: return None
    missing=N-n_on
    # gap structure: runs of missing vertices between present ones
    runs=sorted([int(round(x/g))-1 for x in gaps if int(round(x/g))>1],reverse=True)
    core=f'regular {N}-gon'+(' with its centre' if centre else '')
    if missing==0: return core,N
    if missing==1: return core+' minus a vertex',N
    desc=' + '.join(f'{r} adjacent' if r>1 else '1' for r in runs)
    return core+f' minus {missing} vertices (missing: {desc})',N
def lattice_label(P):
    n=len(P)
    diffs=[(P[j][0]-P[i][0],P[j][1]-P[i][1]) for i in range(n) for j in range(n) if i!=j]
    v=min(diffs,key=lambda d: math.hypot(*d)); lv=math.hypot(*v)
    for kind,ang in (('triangular',60),('square',90)):
        w=(v[0]*math.cos(math.radians(ang))-v[1]*math.sin(math.radians(ang)),v[0]*math.sin(math.radians(ang))+v[1]*math.cos(math.radians(ang)))
        det=v[0]*w[1]-v[1]*w[0]; ok=True
        for d in diffs:
            a=(d[0]*w[1]-d[1]*w[0])/det; b=(v[0]*d[1]-v[1]*d[0])/det
            if abs(a-round(a))>1e-6 or abs(b-round(b))>1e-6: ok=False;break
        if ok: return f'piece of the {kind} lattice'
    return None
def spectrum(P):
    ds=sorted(set(round(math.dist(a,b),6) for a,b in itertools.combinations(P,2))); m=ds[0]
    return ', '.join('%.4g'%(x/m) for x in ds)
out={}
for cell,S in [('10,5',S10['5']),('10,6',S10['6']),('11,6',S11['11,6']),('12,6',S11['12,6'])]:
    L=[]
    for r in S:
        P=[tuple(p) for p in r['pts']]
        pl=polygon_label(P); N=None
        if pl: lab,N=pl
        else:
            lab=lattice_label(P) or ('distances 1 : '+spectrum(P))
        L.append({'pts':P,'label':lab,'N':N,'pat':r['pat']})
    out[cell]=L
    import collections; print(cell,len(L),collections.Counter(x['label'].split(' minus')[0].split(' (')[0] for x in L).most_common(8))
json.dump(out,open('exact_cells.json','w'))
