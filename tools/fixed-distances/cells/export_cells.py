"""Build cells/exact_cells2.json for build_data2.py from cells/rigid_{n}_{d}.json and cells/fam_{n}_{d}.json:
{"n,d": {"rigid":[{pts,label}], "fam":[{pts,dim,col,label}], "stats":{...}}}.  Labels: polygon subsets, lattice pieces,
otherwise the distance spectrum; plus "maximal" / "subset of k of the (n+1)-point sets" when the (n+1,d) cell is exact.
Family labels: the catalog's structural label when a catalog family has the same colouring type, else generic.
Usage: python3 cells/export_cells.py n d [n d ...]"""
import json, sys, os, math, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__)); SCR=os.path.dirname(HERE)
src=open(os.path.join(SCR,'exact_cells.py')).read()
exec(src[:src.index('out={}')].replace("S10=json.load(open('g7count/sets10_by_d.json')); S11=json.load(open('g7count/sets1113_by_d.json'))",""))
sys.path.insert(0,os.path.join(SCR,'r12')); from congr import congruent
A=json.load(open(os.path.join(SCR,'applet_data.json')))
OUT=os.path.join(HERE,'exact_cells2.json')
out=json.load(open(OUT)) if os.path.exists(OUT) else {}
def dvals(P): return len(set(round(math.dist(a,b),7) for a,b in itertools.combinations(P,2)))
import cmath
sys.path.insert(0,os.path.join(SCR,'g7fin','gen')); import patlib
Z=[cmath.exp(2j*math.pi*k/12) for k in range(12)]
R12SUM=sorted({(round((Z[a]+Z[b]).real,9),round((Z[a]+Z[b]).imag,9)) for a in range(12) for b in range(12)})
PHI2=(1+5**.5)**2/4
P10=[(math.cos(2*math.pi*k/5),math.sin(2*math.pi*k/5)) for k in range(5)]+[(math.cos(2*math.pi*k/5+math.pi/5)/PHI2,math.sin(2*math.pi*k/5+math.pi/5)/PHI2) for k in range(5)]
_MCACHE={}
def fits_in(P,M):
    """is P similar to a subset of the point set M?  Anchor the shortest edge of P on a pair of M; a pair is tried only
    if every distance of P, rescaled, occurs among the distances of M (precomputed, rounded)."""
    key=id(M)
    if key not in _MCACHE:
        Mc=[complex(*m) for m in M]
        dset={round(abs(a-b),5) for a,b in itertools.combinations(Mc,2)}
        pairs=[(U,W,abs(W-U)) for U in Mc for W in Mc if U!=W]
        _MCACHE[key]=(Mc,dset,pairs)
    Mc,dset,pairs=_MCACHE[key]
    Pc=[complex(*p) for p in P]
    i,j=min(itertools.combinations(range(len(Pc)),2),key=lambda e: abs(Pc[e[0]]-Pc[e[1]]))
    dmin=abs(Pc[i]-Pc[j]); pd=sorted({abs(a-b)/dmin for a,b in itertools.combinations(Pc,2)})
    for U,W,L in pairs:
        if any(round(x*L,5) not in dset for x in pd): continue
        s=(W-U)/(Pc[j]-Pc[i]); t=U-s*Pc[i]
        for conj in (False,True):
            Q=[s*(p.conjugate() if conj else p)+t for p in Pc]
            if all(min(abs(q-m) for m in Mc)<1e-6 for q in Q): return True
    return False
def special_label(P):
    if fits_in(P,R12SUM): return 'piece of the twelvefold set R₁₂ + R₁₂'
    if fits_in(P,P10): return 'subset of the regular pentagon with the inner pentagon of its pentagram'
    return None
def canon_col(col,n): return patlib.canon(n,tuple(patlib.rgf(list(col))))
def coltype(P,n,tol=1e-6):
    E=list(itertools.combinations(range(n),2)); ds=[math.dist(P[i],P[j]) for i,j in E]; mx=max(ds); vals=[]; col=[]
    for x in ds:
        for k,v in enumerate(vals):
            if abs(x-v)<tol*mx: col.append(k); break
        else: vals.append(x); col.append(len(vals)-1)
    perms=[[E.index(tuple(sorted((p[a],p[b])))) for a,b in E] for p in itertools.permutations(range(n))]
    best=None
    for pe in perms:
        c=[col[pe[i]] for i in range(len(E))]; mp={}; o=[]
        for x in c:
            if x not in mp: mp[x]=len(mp)
            o.append(mp[x])
        t=tuple(o)
        if best is None or t<best: best=t
    return best
args=[int(x) for x in sys.argv[1:]]
for n,d in zip(args[::2],args[1::2]):
    rig=json.load(open(os.path.join(HERE,f'rigid_{n}_{d}.json'))); fam=json.load(open(os.path.join(HERE,f'fam_{n}_{d}.json')))
    sup_key=f'{n+1},{d}'
    sup=None
    if os.path.exists(os.path.join(HERE,f'rigid_{n+1}_{d}.json')):
        sup=[r['pts'] for r in json.load(open(os.path.join(HERE,f'rigid_{n+1}_{d}.json')))]
    elif sup_key in A['cats'] and A['cats'][sup_key].get('complete'):
        sup=[[ [float(x) for x in c] for c in it['coords']] for it in A['cats'][sup_key]['items'] if it['kind']=='rigid']
    def fpkey(P):
        ds=sorted(math.dist(a,b) for a,b in itertools.combinations(P,2)); m=ds[0]
        return tuple(round(x/m,4) for x in ds)
    supindex={}
    if sup is not None:
        for k,Q in enumerate(sup):
            for i in range(len(Q)):
                Qi=[q for j,q in enumerate(Q) if j!=i]
                if dvals(Qi)==d: supindex.setdefault(fpkey(Qi),[]).append((k,Qi))
    R=[]; cnt=collections.Counter(); nmax=0
    for r in rig:
        P=[tuple(p) for p in r['pts']]
        if r.get('rank')==1: lab='collinear: '+spectrum(P)
        else:
            pl=polygon_label(P); lab=(pl[0] if pl else None) or lattice_label(P) or special_label(P) or ('distances 1 : '+spectrum(P))
        if sup is not None:
            ks=[]
            for k,Qi in supindex.get(fpkey(P),[]):
                if k not in ks and congruent(P,Qi): ks.append(k)
            if ks: lab+='; a subset of '+('%d of the %d-point %d-distance sets'%(len(ks),n+1,d) if len(ks)>1 else 'one %d-point %d-distance set'%(n+1,d))
            else: lab+='; maximal (in no %d-point %d-distance set)'%(n+1,d); nmax+=1
        cnt[lab.split(';')[0].split(' (')[0]]+=1
        R.append({'pts':P,'label':lab})
    # families: match catalog labels by colouring type
    key=f'{n},{d}'
    catfam={}
    if key in A['cats']:
        for it in A['cats'][key]['items']:
            if it['kind']=='family':
                P=[[float(x) for x in c] for c in it['coords']]
                catfam[canon_col(coltype(P,n),n)]=it['label']
    F=[]
    for f in fam:
        ct=canon_col(f['col'],n); sizes=sorted([f['col'].count(k) for k in range(d)],reverse=True)
        lab=catfam.get(ct) or f'{f["dim"]}-parameter family; distance classes of sizes '+'·'.join(map(str,sizes))
        if lab.startswith('piece of the triangular lattice'): lab='one-parameter family through a piece of the triangular lattice (the lattice piece is one point of the curve); distance classes of sizes '+'·'.join(map(str,sizes))
        if f.get('components',1)>1: lab+=f' ({f["components"]} components)'
        if f.get('ranks')==[1]: lab+='; collinear'
        F.append({'pts':f['point'],'dim':f['dim'],'col':f['col'],'label':lab})
    C=json.load(open(os.path.join(HERE,f'cert_{n}_{d}.json')))
    famidx={f['idx'] for f in fam}
    ncomplex=sum(1 for r in C.values() if not r['rigid'] and r['idx'] not in famidx)
    out[key]={'rigid':R,'fam':F,'stats':{'nmax':nmax if sup is not None else None,'kinds':dict(cnt.most_common()),'catfam_matched':sum(1 for f in fam if canon_col(f['col'],n) in catfam),'npat':len(C),'ncomplex':ncomplex}}
    print(key,'rigid',len(R),'fam',len(F),'maximal',nmax if sup is not None else '?','kinds',cnt.most_common(5),'catalog families matched',out[key]['stats']['catfam_matched'],'/',len(catfam))
json.dump(out,open(OUT,'w'))
