"""Assemble applet_data.json with per-dimension counts/bounds."""
import collections, json, math, itertools, os, glob
from collections import Counter
import build_data as B
from build_data import item, poly, SRC, same, dvals, pd
import examples as EX
import symfam as SF
from extend import extensions
S3=math.sqrt(3)
def local_dim(P):
    P=[tuple(map(float,p)) for p in P]
    col,vals=SF.colouring(P,1e-7*max(SF.dists(P)))
    G,dim=SF.generic_point_and_dim(P,col,len(vals))
    return dim,G,SF.type_hash(P,col)

import numpy as np
from scipy.optimize import least_squares
_rng=np.random.default_rng(5)
def beautify(G,dim):
    """pick a representative of the family through G with a less extreme distance ratio"""
    n=len(G); E=list(itertools.combinations(range(n),2))
    col,vals=SF.colouring([tuple(p) for p in G],1e-7*max(SF.dists(G))); d=len(vals); NV=2*n-3
    def unpack(z):
        R=[(0.0,0.0),(z[0],0.0)]
        for i in range(2,n): R.append((z[2*i-3],z[2*i-2]))
        return R
    def res(z):
        R=unpack(z); v=z[NV:]
        return [math.dist(R[i],R[j])-v[col[k]] for k,(i,j) in enumerate(E)]+[v[0]-1.0]
    best=G; bestr=min(pd(G))/max(pd(G))
    for _ in range(10):
        z0=np.concatenate([_rng.uniform(-1.5,1.5,NV),_rng.uniform(0.3,2.5,d)])
        r=least_squares(res,z0,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=300)
        if max(abs(r.fun))>1e-11: continue
        Q=unpack(r.x); ds=pd(Q)
        if min(ds)<1e-3*max(ds): continue
        if len(dvals(Q))!=d: continue
        dq,_,_=local_dim(Q)
        if dq!=dim: continue
        rat=min(ds)/max(ds)
        if rat>bestr: best,bestr=Q,rat
    return settle([tuple(p) for p in best],dim)
def settle(G,dim):
    """walk along the family until the shortest/longest distance ratio and the smallest gap between distance
    values are comfortably inside the bounds the applet's flexer uses (same score as in the page's JS)"""
    n=len(G); E=list(itertools.combinations(range(n),2))
    col,vals=SF.colouring([tuple(p) for p in G],1e-7*max(SF.dists(G))); d=len(vals); NV=2*n-3
    ang=math.atan2(G[1][1]-G[0][1],G[1][0]-G[0][0]); co,si=math.cos(-ang),math.sin(-ang)
    Q=[((x-G[0][0])*co-(y-G[0][1])*si,(x-G[0][0])*si+(y-G[0][1])*co) for x,y in G]
    z=[Q[1][0]]
    for i in range(2,n): z+=[Q[i][0],Q[i][1]]
    v=[0.0]*d
    for k,(i,j) in enumerate(E): v[col[k]]=math.dist(Q[i],Q[j])
    z=np.array(z+v); v0=v[0]
    def unpack(z):
        R=[(0.0,0.0),(z[0],0.0)]
        for i in range(2,n): R.append((z[2*i-3],z[2*i-2]))
        return R
    def res(z):
        R=unpack(z); vv=z[NV:]
        return [math.dist(R[i],R[j])-vv[col[k]] for k,(i,j) in enumerate(E)]+[vv[0]-v0]
    def project(z):
        r=least_squares(res,z,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=60); return r.x if max(abs(r.fun))<1e-10 else None
    def score(z):
        R=unpack(z); ds=[math.dist(R[i],R[j]) for i,j in E]; mx=max(ds); vv=sorted(z[NV:])
        sep=min(vv[k+1]-vv[k] for k in range(d-1))/mx if d>1 else 1
        return min(min(ds)/mx/0.1, sep/0.012)
    def tangent(z):
        f0=np.array(res(z)); h=1e-7; J=np.zeros((len(f0),len(z)))
        for i in range(len(z)):
            zz=z.copy(); zz[i]+=h; J[:,i]=(np.array(res(zz))-f0)/h
        U,S,Vt=np.linalg.svd(J); return Vt[-dim:]
    sc=score(z); tries=0
    while sc<1.4 and tries<400:
        tries+=1; best=None; bs=sc
        for u in tangent(z):
            for sg in (1,-1):
                zt=project(z+sg*0.03*u)
                if zt is None: continue
                st=score(zt)
                if st>bs: bs,best=st,zt
        if best is None: break
        z=best; sc=bs
    R=unpack(z)
    if len(dvals(R))!=d: return G
    return [tuple(p) for p in R]

import ast, re
def struct_label(st):
    m=re.match(r"\((.*?)\)(.*)",st); tup=ast.literal_eval('('+m.group(1)+')'); rest=m.group(2).strip()
    kind,mm,c,a,b,g,p=tup
    def gon(k,mm): return ('%d regular %d-gon%s'%(k,mm,'' if k==1 else 's')) if mm>2 else ('%d antipodal pair%s'%(k,'' if k==1 else 's'))
    if kind=='C':
        parts=([ 'the centre'] if c else [])+[gon(a,mm)+' (free radii and rotations)']
        lab='%d-fold rotational symmetry: '%mm+' + '.join(parts)
    elif kind=='D':
        parts=(['the centre'] if c else [])
        if a: parts.append(gon(a,mm)+' on the mirror axes')
        if b: parts.append(gon(b,mm)+' on the other axes')
        if g: parts.append('%d generic orbit%s of %d points'%(g,'' if g==1 else 's',2*mm))
        lab='dihedral symmetry D%d: '%mm+' + '.join(parts)
    else:
        lab='mirror symmetry: %d point%s on the axis + %d mirrored pair%s'%(c,'' if c==1 else 's',a,'' if a==1 else 's')
    k=re.search(r"\+(\d+) coincidence",rest)
    if k: lab+=', with %s extra coincidence%s imposed'%(k.group(1),'' if k.group(1)=='1' else 's')
    return lab

def listing_note(nrigid,nfam,complete):
    """describe what the catalog list contains, matching its actual content"""
    if nrigid and nfam:
        t='Rigid configurations first, then one representative of each family (a family = a realizable pattern of equal pairs; its dimension is measured in the full configuration space).'
    elif nfam:
        t='One representative of each family (a family = a realizable pattern of equal pairs; its dimension is measured in the full configuration space). Hover or tap a card to flex it.'
    elif nrigid:
        t='All configurations shown are rigid.'
    else:
        t=''
    if not complete: t+=' The list is a lower bound.'
    return t
class Pool:
    """examples for one cell, deduped: rigid by similarity, families by type"""
    def __init__(s): s.rigid=[]; s.fam=[]; s.exact={}
    def add(s,P,label,src,exact=None):
        ds=pd(P)
        if min(ds)<1e-3*max(ds): return   # near-degenerate escape of a numerical solver
        dim,G,h=local_dim(P)
        if dim==0:
            if any(same(P,Q) for Q,_,_ in s.rigid): return
            if exact: s.exact[len(s.rigid)]=exact
            s.rigid.append((P,label,src))
        else:
            if any(hh==h and dd==dim for _,_,_,hh,dd in s.fam): return
            G=beautify(G,dim)
            s.fam.append((G,label,src,h,dim))
    def counts(s):
        c={}
        if s.rigid: c[0]=len(s.rigid)
        for *_,dim in s.fam: c[dim]=c.get(dim,0)+1
        return c
_item=item
def item(P,label,kind='rigid',dim=0,**kw):
    if kind=='family' and dim>=1:
        try: P=settle([tuple(map(float,p)) for p in P],dim)
        except Exception as e: print('settle failed',label[:40],e)
    return _item(P,label,kind,dim,**kw)
# ---------------- stratification: incidence of the one-parameter families with the rigid sets ----------------
def add_stratification(cells,cats):
    import itertools as _it
    E5=list(_it.combinations(range(5),2))
    perms5=[[E5.index(tuple(sorted((p[a],p[b])))) for a,b in E5] for p in _it.permutations(range(5))]
    def ccol(rgs):
        best=None
        for pe in perms5:
            c=[rgs[pe[i]] for i in range(10)]; mp={}; out=[]
            for x in c:
                if x not in mp: mp[x]=len(mp)
                out.append(mp[x])
            t=tuple(out)
            if best is None or t<best: best=t
        return best
    for d in (4,5):
        fn=f'strat_graph_{d}.json'
        if not os.path.exists(fn) or f'5,{d}' not in cats: continue
        G=json.load(open(fn)); verts={v['id']:v for v in G['vertices']}
        fam_items={}
        for k,it in enumerate(cats[f'5,{d}']['items']):
            if it['kind']=='family' and it['dim']==1: fam_items[ccol(it['cls'])]=it
        tcol={int(t):ccol(c) for t,c in G['types'].items()}
        per_type={}; per_vertex={}
        for e in G['edges']:
            per_type.setdefault(e['type'],[]).append(e)
            for end in e['ends']:
                if end[0]=='v': per_vertex.setdefault(end[1],[]).append(e)
        nvv=sum(1 for e in G['edges'] if all(x[0]=='v' for x in e['ends'])); ninf=sum(1 for e in G['edges'] if any(x[1]=='inf' for x in e['ends']))
        ndeg=sum(1 for e in G['edges'] if any(x[1]=='degenerate' for x in e['ends']))
        # catalog numbers of the (5,d-1) families, for junction points
        lowfam={}
        if f'5,{d-1}' in cats:
            lowcols={r['idx']:tuple(r['col']) for r in json.load(open(f'cert5{d-1}.json'))} if os.path.exists(f'cert5{d-1}.json') else {}
            lowbycol={ccol(c):t for t,c in lowcols.items()}
            for k,it in enumerate(cats[f'5,{d-1}']['items']):
                if it['kind']=='family': 
                    t=lowbycol.get(ccol(it['cls']))
                    if t is not None: lowfam[t]=k+1
        def vname(i):
            v=verts[i]
            if v['item'] is None:
                import re as _re
                m=_re.search(r'type (\d+)',v['label'] or '')
                t=int(m.group(1)) if m else None
                return f"interior point of {d-1}-distance family #{lowfam[t]}" if t in lowfam else f"interior point of a {d-1}-distance family"
            if v['cell']==d-1 and (v['label'] or '').startswith('rigid'): return f"{d-1}-distance set #{v['item']+1}"
            if v['cell']==d-1: return (v['label'] or 'rigid set').split(';')[0]
            return f"{v['label']} ({v['cell']} distances)"
        matched=0
        for t,edges in per_type.items():
            it=fam_items.get(tcol[t])
            if it is None: continue
            matched+=1
            parts=[]
            for e in edges:
                a,b=e['ends']
                def nm(x): return vname(x[1]) if x[0]=='v' else ('two points merge' if x[1]=='degenerate' else 'one distance → 0 (infinity)')
                parts.append(f'{nm(a)} ↔ {nm(b)}' if a[0]=='v' and b[0]=='v' else f'{nm(a) if a[0]=="v" else nm(b)} → {nm(b) if a[0]=="v" else nm(a)}')
            it['strat']=f'{len(parts)} arc{"s" if len(parts)!=1 else ""} of this family:\n'+'\n'.join('• '+p for p in parts)
        for i,edges in per_vertex.items():
            v=verts[i]
            if v['cell'] is None: continue
            key=f"5,{v['cell']}"
            if key in cats and v['item'] is not None and v['item']<len(cats[key]['items']):
                it=cats[key]['items'][v['item']]
                it['strat']=(it.get('strat','')+(' · ' if it.get('strat') else '')+f"endpoint of {len(edges)} arcs of {d}-distance families")
        # connectivity and first Betti number of the graph (rigid sets + arcs between them)
        par={i:i for i in verts}
        def find(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        for e in G['edges']:
            a,b=e['ends']
            if a[0]=='v' and b[0]=='v': par[find(a[1])]=find(b[1])
        comps=len({find(i) for i in verts}); b1=nvv-len(verts)+comps
        cells[f'5,{d}']['strat']={'arcs_vv':nvv,'arcs_degenerate':ndeg,'arcs_inf':ninf,'vertices':len(verts),'components':comps,'b1':b1}
        nrig=sum(1 for it in cats[f'5,{d}']['items'] if it['kind']=='rigid')
        junc=sum(1 for v in verts.values() if v['cell']==d-1 and v['item'] is None); nrv=len(verts)-junc
        byc=collections.Counter(v['cell'] for v in verts.values() if v['item'] is not None)
        lower=('the 34 three-distance sets and the regular pentagon' if d==4 else
               ', '.join(f"{byc[c]} with {c} distances" for c in sorted(byc,reverse=True)) + (f'; a further {junc} endpoints are interior points of one-parameter families of {d-1} distances (the curve meets the family there, not at a rigid set)' if junc else ''))
        cells[f'5,{d}']['note']+=(f' Structure of the solution set: the one-parameter families are curves whose ends are rigid configurations with fewer distances (two distance values merge there), degenerate configurations (two points merge) or infinity (one distance becomes negligible). Together with their endpoints — {nrv} rigid sets: {lower} — they form a graph with {nvv} arcs between rigid sets ({"connected" if comps==1 else f"{comps} components"}, first Betti number {b1}), {ndeg} arcs ending in a degenerate configuration and {ninf} escaping to infinity; every rigid set with fewer distances is the limit of some family.'+(f' The {nrig} rigid {d}-distance sets are isolated points of the solution set.' if nrig else '')+' Each card below lists the arcs through it.')
        print(f'stratification (5,{d}): {matched} families annotated, arcs vv {nvv}, deg {ndeg}, inf {ninf}')
    # degeneration posets: boundary strata of the higher-dimensional families (exact for d = 5; sampled for d = 6, 7)
    for d in (5,6,7):
        fn=f'degen_poset_{d}.json'
        if not os.path.exists(fn) or f'5,{d}' not in cats: continue
        PO=json.load(open(fn)); cols={r['idx']:tuple(r['col']) for r in json.load(open(f'cert5{d}.json'))}
        bycol={ccol(c):int(t) for t,c in cols.items() if t in {int(k) for k in PO}}
        ann=0; srcs=set()
        for it in cats[f'5,{d}']['items']:
            if it['kind']!='family' or (d==5 and it['dim']<2): continue
            t=bycol.get(ccol(it['cls']))
            if t is None: continue
            v=PO[str(t)]; srcs.add(v['src'])
            cnt={}
            for x in v['to']:
                if x=='degenerate': key='degenerate configurations (two points merge)'
                else:
                    dn,idx,st=x; key=f'{"rigid " if st=="rigid" else ""}{dn}-distance {"sets" if st=="rigid" else "families"}'
                cnt[key]=cnt.get(key,0)+1
            parts=[f'{n} {k}' if not k.startswith('degenerate') else k for k,n in sorted(cnt.items(),key=lambda kv:(-len(kv[0]),kv[0]))]
            it['strat']=('boundary strata' + (' (sampled)' if v['src']=='sampled' else ' (exact)') + ': ' + '; '.join(parts)) if parts else ''
            ann+=1
        if d>=6 or any(it.get('strat','').startswith('boundary') for it in cats[f'5,{d}']['items']):
            cells[f'5,{d}']['note']+=(' Boundary structure: for each family the strata of fewer distances in its closure were determined from the saturated ideal with one coincidence equation added — '
                +('by exact decomposition into irreducible components (Singular), with the type of each component read off a certified real point' if srcs=={'exact'} else
                  'exactly for some families and, where the decomposition did not finish, by random slices through the boundary (certified real points, so every listed stratum is genuinely in the closure, but a stratum can be missed)')
                +'; the cards list them.')
        print(f'degeneration poset (5,{d}): {ann} families annotated ({",".join(sorted(srcs))})')

def main():
    cells={}; cats={}
    def cell(n,d,dims,status,marks=(),note='',finite=None,extra_line=None):
        cells[f'{n},{d}']={'n':n,'d':d,'status':status,'dims':{str(k):v for k,v in sorted(dims.items())},'marks':list(marks),'note':note,'finite':finite,
                           'krull':(2*n-4-math.comb(n,2)+d) if isinstance(n,int) else None}
        if extra_line: cells[f'{n},{d}']['extra_line']=extra_line
    def cat(n,d,items,complete,title,note='',sources=()):
        cats[f'{n},{d}']={'items':items,'complete':complete,'title':title,'note':note,'sources':[SRC[k] for k in sources]}
    X=lambda lo,hi=None: {'lo':lo,'hi':hi}
    E=lambda v: {'lo':v,'hi':v}
    # ---------- n = 3, 4 ----------
    cell(3,1,{0:E(1)},'exact',note='The equilateral triangle.'); cat(3,1,[item([(0,0),(1,0),(.5,S3/2)],'equilateral triangle')],True,'The equilateral triangle')
    cell(3,2,{1:E(1)},'exact',note='Isosceles triangles: one parameter (the apex angle).'); cat(3,2,[item([(0,0),(1,0),(.5,1.1)],'isosceles triangle (apex angle free)','family',1)],True,'Isosceles triangles')
    cell(3,3,{2:E(1)},'exact',note='Scalene triangles: two parameters.'); cat(3,3,[item([(0,0),(1,0),(.3,.8)],'scalene triangle','family',2)],True,'Scalene triangles')
    cell(4,2,{0:E(6)},'exact',marks='†',note='The six two-distance quadruples (Einhorn–Schoenberg 1966); re-verified by exact computation (5 colouring types, 6 rigid configurations up to similarity).'); cat(4,2,[item(P,l) for l,P in B.two],True,'The six two-distance quadruples',sources=('ES',))
    cell(4,3,{1:E(9)},'exact',marks='‡',note='Nine one-parameter families, one per realizable 3-colouring of the edges of K₄ (all nine types are realizable; certified by exact computation); no rigid configurations (every component has dimension ≥ 2·4−4−6+3 = 1). The six two-distance quadruples are where families meet.')
    cat(4,3,[item(P,l,'family',1) for l,P in B.fam43],True,'Nine one-parameter families','One representative per family; the family is traced by moving the free parameter.',sources=('us',))
    S4=json.load(open('strata_4_4_1.json'))['realizable']
    items=[]
    for k,(l,P) in enumerate(B.fam44):
        if P is None: P=sorted(S4[k]['sample'].items())[0][1]
        items.append(item(P,l,'family',2))
    cell(4,4,{2:E(7)},'exact',marks='‡',note='Seven two-parameter families (one per realizable 4-colouring of K₄). All 7 colouring types with four colours are realizable and each has dimension exactly 2 (certified by exact computation).')
    cat(4,4,items,True,'Seven two-parameter families',sources=('us',))
    S5=json.load(open('strata_4_5_1.json'))['realizable']
    cell(4,5,{3:E(2)},'exact',marks='‡',note='One coincidence among the six distances: the equal pair is adjacent or opposite (certified by exact computation).')
    cat(4,5,[item(sorted(S5[k]['sample'].items())[0][1],l,'family',3) for k,(l,_) in enumerate(B.fam45)],True,'Two three-parameter families',sources=('us',))
    cell(4,6,{4:E(1)},'exact',marks='‡',note='Generic quadruples.'); cat(4,6,[item(P,l,'family',4) for l,P in B.fam46],True,'Generic quadruples')
    # ---------- n = 5 ----------
    cell(5,2,{0:E(1)},'exact',marks='†',note='The regular pentagon (Einhorn–Schoenberg 1966); re-verified by exact computation (17 colouring types, one rigid configuration).'); cat(5,2,[item(poly(5,range(5)),'regular pentagon',Nhint=5)],True,'The regular pentagon',sources=('ES',))
    FP=json.load(open('five_page.json')); FIVE=[x['pts'] for x in json.load(open('five_all.json'))]; SIX=json.load(open('six.json'))
    HEPT=poly(7,range(7)); HEXC=[(0,0)]+poly(6,range(6))
    def precise(pts,pool):
        for P in pool:
            if same(pts,P,1e-4): return P
        raise SystemExit('no match')
    items=[]
    for g in FP:
        if g['title'] in ('Six points','Seven points'): continue
        for it in g['items']: items.append(item(precise(it['pts'],FIVE),f"{g['title']}: {it['label']}",note=it['dist']))
    cell(5,3,{0:E(34)},'exact',marks='†',note='Finite: 34 rigid configurations (Shinohara 2004, Theorem 1), re-verified by exact computation (124 colouring types of K₅ with three colours; 35 certified real roots, 34 configurations up to similarity) and by two numerical computations (extension of the four-point families, and exhaustive colouring enumeration). Each of our 34 sets was matched to one of Shinohara’s figures 501–534 by its pair of distance ratios (α, β); all 34 correspond (his table gives 2 sin 4π/9 for the β of fig. 524, where the construction in the text yields 2 sin 7π/18).')
    cat(5,3,items,True,'The 34 five-point three-distance sets','Grouped by the figure they come from.',sources=('Sh04','us'))
    def numeric_cell(d,N):
        R=list(N.values()); ncls=json.load(open(f'strata_5_{d}_1.json'))['classes']
        rig=[]; fams={}; items=[]
        for r in R:
            sizes=sorted([r['col'].count(k) for k in range(d)],reverse=True)
            for P in r['rigid']:
                rig.append(P); items.append(item(P,'rigid; distance classes of sizes '+'·'.join(map(str,sizes)),'rigid',0))
            for k,P in sorted(r['sample'].items()):
                fams[int(k)]=fams.get(int(k),0)+1; items.append(item(P,f'{int(k)}-parameter family; distance classes of sizes '+'·'.join(map(str,sizes)),'family',int(k)))
        dims={}
        if rig: dims[0]=E(len(rig))
        for k,v in fams.items(): dims[k]=E(v)
        e=2*5-4-10+d
        note=f'{len(R)} of the {ncls} colourings of K₅ with {d} colours (up to symmetry) are realizable (numerical realizability test). Every component has dimension ≥ {e}' + (', so there are no rigid configurations' if e>=1 else '') + '. Exact certification of this cell is in progress.'
        cell(5,d,dims,'exact',marks='*',note=note)
        cat(5,d,items,True,f'Five points, {d} distances',listing_note(len(rig),sum(fams.values()),True)+' Types are described by the sizes of the distance classes.',sources=('us',))
    import os
    NCLS={4:371,5:443,6:249,7:76}
    for d in (4,5,6,7):
        # certified classification (msolve + Singular) when complete; numerical samples only supply representatives
        cert=json.load(open(f'cert5{d}.json')) if os.path.exists(f'cert5{d}.json') else []
        e=2*5-4-10+d
        stub={}
        for r in cert:
            if r.get('unresolved'):
                fn=f'stubborn_{d}_{r["idx"]}.json'
                if os.path.exists(fn):
                    sb=json.load(open(fn))
                    if sb['dim_saturated']==e and sb['slice_hits']>0 and sb['ranks'] and all(rk==7+d-1-e for rk,_ in sb['ranks']) : stub[r['idx']]=sb
        certified=len(cert)==NCLS[d] and all((not r.get('unresolved')) or r['idx'] in stub for r in cert)
        N={tuple(r['col']):r for r in json.load(open(f'strata_5_{d}_1.json'))['realizable']}
        if d==4:   # exact similarity classes of the rigid roots (dedup54.py: group orbits with interval arithmetic)
            EXACT={x['idx']:x['rigid'] for x in json.load(open('cert54_rigid_exact.json'))}
            for r in cert: r['rigid']=EXACT.get(r['idx'],[])
        if not certified:
            numeric_cell(d,N); continue
        ncls=len(cert)
        rig=[]; fams={}; items=[]; nreal=0
        for r in cert:
            col=tuple(r['col']); sizes=sorted([col.count(k) for k in range(d)],reverse=True)
            realdims=set()
            if r['idx'] in stub:
                realdims.add(e); r=dict(r,primes=[{'real':True,'dim':e,'point':stub[r['idx']]['points'][0]}])
            for p in r['primes']:
                if p.get('real') and p['dim']>=1: realdims.add(p['dim'])
            if d==4:
                for p in r['primes']:
                    if p['dim']==1 and not p['degenerate'] and p.get('real_slices',0)>0: realdims.add(1)
            if r['rigid'] or realdims: nreal+=1
            for P in r['rigid']:
                rig.append(P); items.append(item(P,'rigid; distance classes of sizes '+'·'.join(map(str,sizes)),'rigid',0))
            for k in sorted(realdims):
                fams[k]=fams.get(k,0)+1
                rep=None
                nr=N.get(col)
                if nr and str(k) in nr['sample']: rep=nr['sample'][str(k)]
                if rep is None:
                    for p in r['primes']:
                        if p.get('real') and p['dim']==k and 'point' in p: rep=p['point']; break
                if rep is None: print('no representative for type',d,r['idx'],'dim',k); continue
                items.append(item(rep,f'{k}-parameter family; distance classes of sizes '+'·'.join(map(str,sizes)),'family',k))
        dims={}
        if rig: dims[0]=E(len(rig))
        for k,v in fams.items(): dims[k]=E(v)
        e=2*5-4-10+d
        note=(f'Certified by exact computation: for each of the {ncls} colourings of K₅ with {d} colours (up to symmetry) the polynomial system was solved exactly '
              f'(msolve for zero-dimensional systems with certified real-root isolation; Singular minimal primes for positive-dimensional ones, each component tested for real points). '
              f'{nreal} types are realizable. Every component has dimension ≥ {e}' + (', so there are no rigid configurations' if e>=1 else '') + '.')
        if stub: note+=f' For {len(stub)} type{"s" if len(stub)>1 else ""} the primary decomposition did not finish; there the certificate is weaker: the saturated ideal (inequations enforced) has dimension exactly {e}, so every component has dimension {e}, real points were certified by random slices, and the Jacobian has full rank at them numerically (a smooth point of an {e}-dimensional component).'
        hidden={5:'A search for hidden rigid configurations (real points of a positive-dimensional type that are isolated in the real locus, which the dimension count cannot see) found none: for 384 of the 443 types the Jacobian of the ten equations has full rank on the whole nondegenerate variety; for the other 59 the real rank-deficient points were shown to be smooth points of their components, or the components’ singular loci have no real points. For four types (380, 381, 426, 427), whose decomposition does not finish, the 20 candidate points were settled by a certified local test (interval arithmetic on the isolating boxes): at each, a smooth surface cut out by generators with independent gradients contains the curve, and the intrinsic Hessian of a further generator restricted to it is indefinite, so by the Morse lemma the point is a node of two real branches — not isolated.',
                6:'A search for hidden rigid configurations (real-isolated points of a positive-dimensional type) found none: full Jacobian rank everywhere for 216 of the 249 types; for the other 33 every component is smooth at the real rank-deficient points or has a singular locus without real points.',
                7:'A search for hidden rigid configurations (real-isolated points of a positive-dimensional type) found none: full Jacobian rank everywhere for 70 of the 76 types; for the other 6 every component is smooth at the real rank-deficient points or has a singular locus without real points.'}
        if d in hidden: note+=' '+hidden[d]
        if d==4: note+=' Rigid configurations are counted up to similarity: the real roots of each type were grouped into orbits of the colour-permuting automorphisms of the colouring (with rescaling) by interval arithmetic on the certified isolating boxes, and the real-root counts were confirmed by Sturm sequences on msolve’s rational parametrizations. (A first version of this cell reported 481: its dedup compared distances without rescaling and counted thirteen configurations twice.) The earlier random search had the right total but had misclassified two rigid types as families, and had called the equally spaced collinear quintuple rigid when it in fact lies on a one-parameter family (the lattice piece 0, −u, v, u+v, −u−v with the angle between u and v free).'
        cell(5,d,dims,'exact',marks='‡',note=note)
        cat(5,d,items,True,f'Five points, {d} distances',listing_note(len(rig),sum(fams.values()),True)+' Types are described by the sizes of the distance classes.',sources=('cert',))
    # ---------- (6,3), (7,3) ----------
    six=[g for g in FP if g['title']=='Six points'][0]['items']; seven=[g for g in FP if g['title']=='Seven points'][0]['items']
    cell(6,3,{0:E(9)},'exact',note='Nine rigid configurations: Shinohara (2004, Theorem 2) lists the six maximal six-point sets (his figs. 601–606); the other three are the six-point subsets of the two seven-point sets (the regular hexagon, the heptagon minus a vertex, the hexagon-with-centre minus a vertex). Our enumeration by extension of the 34 five-point sets gives the same nine.')
    cat(6,3,[item(precise(it['pts'],SIX),it['label'],note=it['dist']) for it in six],True,'The nine six-point three-distance sets',sources=('Sh04','us'))
    cell(7,3,{0:E(2)},'exact',note='The regular heptagon and the regular hexagon with its centre; g(3) = 7.')
    cat(7,3,[item(precise(it['pts'],[HEPT,HEXC]),it['label'],note=it['dist'],Nhint=7 if 'hept' in it['label'] else 6) for it in seven],True,'The two seven-point three-distance sets',sources=('EF','Sh04'))
    # ---------- exact maximal cells ----------
    L94=json.load(open('lattice_9_4.json'))['tri']
    def uniq(L):
        out=[]
        for P in L:
            if not any(same(P,Q) for Q in out): out.append(P)
        return out
    L94=uniq(L94)
    def lat_label(P):
        m=min(pd(P)); deg=max(sum(1 for q in P if abs(math.dist(p,q)-m)<1e-6) for p in P)
        return 'piece of the triangular lattice: a regular hexagon, its centre, and two adjacent points of the next ring' if deg>=6 else 'piece of the triangular lattice: rows of 1, 3, 3 and 2 points'
    FOUR=[(1,0),(0,1),(-1,0),(0,-1),(0,-S3),((S3-1)/2,-(S3-1)/2),(-(S3-1)/2,-(S3-1)/2),((S3+1)/2,-(S3+1)/2),(-(S3+1)/2,-(S3+1)/2)]
    it4=item(FOUR,'three equilateral triangles with a common centre and parallel sides (Erdős–Fishburn’s description; 3-fold symmetry about (0, −1/√3)). Equivalently: square ABCD (side √2); on sides AD and CD the outward equilateral apexes H, I and the inward apexes F, G of isosceles triangles with apex angle 150°; E on the axis beyond D with DE = DF = √3 − 1. Distances √3−1, √2, 2, √3+1 (the 12-gon chords c₁, c₂, c₃, c₅).')
    it4['exact']=[['1','0'],['0','1'],['−1','0'],['0','−1'],['0','−√3'],['(√3−1)/2','−(√3−1)/2'],['−(√3−1)/2','−(√3−1)/2'],['(√3+1)/2','−(√3+1)/2'],['−(√3+1)/2','−(√3+1)/2']]
    it4['frame']='square ABCD inscribed in the unit circle, D at the bottom'
    cell(9,4,{0:E(4)},'exact',note='g(4) = 9. Erdős–Fishburn (1996, Theorem 1): exactly four configurations — the regular nonagon, two pieces of the triangular lattice, and three concentric equilateral triangles with parallel sides (the only known maximum set for any k that is neither a regular polygon, with or without centre, nor a lattice piece). Finite, so no families.',finite=True)
    cat(9,4,[item(poly(9,range(9)),'regular nonagon',Nhint=9)]+[item(P,lat_label(P)) for P in L94]+[it4],True,'The four nine-point four-distance sets',sources=('EF',))
    # (8,4): Shinohara 2008, Theorem 1.2(a)
    NINE=[poly(9,range(9))]+L94+[FOUR]
    h=.5+S3/2; FIGE=[(.5,.5),(-.5,.5),(-.5,-.5),(.5,-.5),(0,h),(0,-h),(h,0),(-h,0)]
    eight=[(poly(8,range(8)),'regular octagon',8),(poly(7,range(7))+[(0,0)],'regular heptagon with its centre',7),(FIGE,'square with the four outward equilateral apexes (Shinohara’s fig. 1(e)); distances 1, √2, √(2+√3), 1+√3',None)]
    nine_names=['the regular nonagon','the lattice piece (hexagon, centre, two adjacent outer points)','the lattice piece (rows 1, 3, 3, 2)','the three concentric triangles']
    nine_names=[nine_names[0]]+[('the lattice piece (hexagon, centre, two adjacent outer points)' if 'hexagon' in lat_label(P) else 'the lattice piece (rows 1, 3, 3, 2)') for P in L94]+[nine_names[3]]
    for P,name in zip(NINE,nine_names):
        for i in range(9):
            Q=[p for j,p in enumerate(P) if j!=i]
            if len(dvals(Q))!=4 or any(same(Q,R) for R,_,_ in eight): continue
            eight.append((Q,f'{name} minus one point',9 if name.startswith('the regular') else None))
    assert len(eight)==15, len(eight)
    cell(8,4,{0:E(15)},'exact',note='Shinohara (2008, Theorem 1.2): every eight-point four-distance set is the regular octagon, the regular heptagon with its centre, the square with four outward equilateral apexes, or an eight-point subset of one of the four nine-point four-distance sets. Deleting one point from each nine-point set and keeping the subsets with exactly four distances gives twelve more up to congruence (one from the nonagon, eight from the two lattice pieces, three from the concentric triangles): 15 in all. Our earlier extension search had found 14 of them.',finite=True)
    cat(8,4,[item(P,lab,Nhint=N) if N else item(P,lab) for P,lab,N in eight],True,'The fifteen eight-point four-distance sets',sources=('Sh08','EF'))
    L125=uniq(json.load(open('lattice_12_5.json'))['tri'])
    cell(12,5,{0:E(1)},'exact',note='g(5) = 12; the unique maximum five-distance set is a piece of the triangular lattice (Shinohara 2008).',finite=True)
    cat(12,5,[item(P,'piece of the triangular lattice: rows of 3, 4, 3 and 2 points') for P in L125],True,'The unique twelve-point five-distance set',sources=('Sh08',))
    # (7,4): Lan–Wei 2013, Theorem 8 — the 42 sets reconstructed from Fig. 4 and the proof (seven4.py)
    S7=json.load(open('seven4.json')); assert len(S7)==42
    def exact_sub(P_exact,P_num,Q):
        keep=[min(range(len(P_num)),key=lambda i: math.dist(P_num[i],q)) for q in Q]
        assert len(set(keep))==len(Q) and all(math.dist(P_num[i],q)<1e-9 for i,q in zip(keep,Q))
        return [P_exact[i] for i in keep]
    FOUR_EX=[['1','0'],['0','1'],['−1','0'],['0','−1'],['0','−√3'],['(√3−1)/2','−(√3−1)/2'],['−(√3−1)/2','−(√3−1)/2'],['(√3+1)/2','−(√3+1)/2'],['−(√3+1)/2','−(√3+1)/2']]
    FIGE_EX=[['1/2','1/2'],['−1/2','1/2'],['−1/2','−1/2'],['1/2','−1/2'],['0','(1+√3)/2'],['0','−(1+√3)/2'],['(1+√3)/2','0'],['−(1+√3)/2','0']]
    PHI2=(1+5**.5)**2/4
    P10=poly(5,range(5))+[(math.cos(2*math.pi*k/5+math.pi/5)/PHI2,math.sin(2*math.pi*k/5+math.pi/5)/PHI2) for k in range(5)]
    P10_EX=[[f'cos(2π·{k}/5)' if k else '1',f'sin(2π·{k}/5)' if k else '0'] for k in range(5)]+[[f'(3−√5)/2·cos({2*k+1}π/5)',f'(3−√5)/2·sin({2*k+1}π/5)'] for k in range(5)]
    LW741_EX=[['0','0'],['1','0'],['1/2','−√3/2'],['−√3/2','1/2'],['1/2','1+√3/2'],['−(1+√3)/2','(1+√3)/2'],['−1−√3/2','1/2']]
    LW742_EX=[['0','0'],['1','0'],['1/2','−√3/2'],['0','1'],['−√3/2','1/2'],['(1−√3)/2','−(1+√3)/2'],['−(1+√3)/2','(1−√3)/2']]
    groups={'poly':'regular polygon: ','ef9':'','lat':'','sh8':'','dp10':'','sqlat':'','lw741':'','lw742':''}
    counters={}
    items74=[]
    for s7 in S7:
        P=[tuple(p) for p in s7['pts']]; src=s7['src']; lab=s7['label']; ex=None; N=None
        if src=='poly':
            N=9 if 'nonagon' in lab else (8 if 'octagon' in lab else 7)
            if 'nonagon' in lab: lab='regular nonagon minus two vertices (Lan–Wei figs. 701–703, 708: four versions, by the gap between the missing vertices)'
            elif 'octagon' in lab: lab='regular octagon minus a vertex (Lan–Wei fig. 707)'
            else: lab='regular heptagon with its centre, minus a vertex (Lan–Wei fig. 704)'
        elif src=='ef9': ex=exact_sub(FOUR_EX,FOUR,P); lab='seven points of the three concentric equilateral triangles (the nine-point four-distance set; Lan–Wei figs. 705, 709–716)'
        elif src=='sh8': ex=exact_sub(FIGE_EX,FIGE,P); lab='seven points of the square with four outward equilateral apexes (Lan–Wei figs. 717, 718)'
        elif src=='dp10': ex=exact_sub(P10_EX,P10,P); lab='seven points of a regular pentagon and the inner pentagon of its pentagram (the ten-point five-distance set; Lan–Wei figs. 719, 720; maximal)'
        elif src=='lat': lab='piece of the triangular lattice (Lan–Wei figs. 706, 721–739)'
        elif src=='sqlat': lab='piece of the square lattice: distances 1, √2, 2, √5 (Lan–Wei fig. 740; maximal)'
        elif src=='lw741': ex=LW741_EX; lab='distances 1, √(2+√3), 1+√3, √(5+2√3) (Lan–Wei fig. 741; maximal)'
        elif src=='lw742': ex=LW742_EX; lab='distances 1, √2, √(2+√3), √(4+√3) (Lan–Wei fig. 742; maximal)'
        it=item(P,lab,Nhint=N) if N else item(P,lab)
        if ex: it['exact']=ex
        items74.append(it)
    cell(7,4,{0:E(42)},'exact',note='Lan and Wei (2013, Theorem 8): exactly 42 seven-point four-distance sets, all drawn in their Fig. 4. The list here was reconstructed from their descriptions (four versions of the nonagon minus two vertices, the octagon minus a vertex, the heptagon-with-centre minus a vertex, nine subsets of the nine-point set of three concentric triangles, twenty pieces of the triangular lattice, two subsets of the square with four apexes, two subsets of the ten-point five-distance set, a square-lattice piece, and two sets given by their distance patterns), each verified to have exactly four distances; the count 42 is reproduced. Five are maximal (not contained in a larger four-distance set). Our earlier extension search had found 15 of them.',finite=True)
    cat(7,4,items74,True,'The 42 seven-point four-distance sets',sources=('LW13',))
    # (11,5): Wei 2011, Theorem 13
    L115=uniq(json.load(open('lattice_11_5.json')))
    assert len(L115)==3
    def lab115(P):
        T=L125[0]
        for i in range(12):
            if same([p for j,p in enumerate(T) if j!=i],P):
                k=sum(1 for j in range(12) if j!=i and abs(math.dist(T[i],T[j])-1)<1e-6)
                return f'the twelve-point five-distance set minus a point that had {k} nearest neighbours'
        raise SystemExit('11-set not a subset of the 12-set')
    cell(11,5,{0:E(4)},'exact',note='Wei (2011, Theorem 13): exactly four eleven-point five-distance sets — the regular 11-gon and three pieces of the triangular lattice (all three are subsets of the unique twelve-point five-distance set). The lattice pieces were re-derived here by a clique search in the lattice with distance set {1, √3, 2, √7, 3}.',finite=True)
    cat(11,5,[item(poly(11,range(11)),'regular 11-gon',Nhint=11)]+[item(P,lab115(P)) for P in L115],True,'The four eleven-point five-distance sets',sources=('Wei11',))
    L136=uniq(json.load(open('lattice_13_6.json'))['tri'])
    cell(13,6,{0:E(3)},'exact',note='g(6) = 13; exactly three configurations (Wei 2012, confirmed by Szöllősi–Östergård 2020).',finite=True)
    cat(13,6,[item(poly(13,range(13)),'regular 13-gon',Nhint=13),item(poly(12,range(12))+[(0,0)],'regular 12-gon with its centre',Nhint=12)]+[item(P,'piece of the triangular lattice: a regular hexagon, its centre and the six points of the next ring') for P in L136],True,'The three thirteen-point six-distance sets',sources=('Wei12','SO20'))
    # ---------- bounds cells ----------
    gmax={1:3,2:5,3:7,4:9,5:12,6:13,7:16}; fmin={3:1,4:2,5:2,6:3,7:3,8:4,9:4,10:5,11:5,12:5,13:6,14:7,15:7,16:7}
    SYM={}
    for fn in glob.glob('symfam_*.json'):
        for k,v in json.load(open(fn)).items(): SYM.setdefault(k,[]).extend(v)
    EXT=json.load(open('ext4.json')); EXT.update(json.load(open('ext5.json')))
    R12=json.load(open('r12sum.json')) if os.path.exists('r12sum.json') else {}
    EXACT=json.load(open('exact_cells.json'))   # (10,5), (10,6), (11,6), (12,6): complete lists from the pattern enumeration + real root isolation
    EXACT_INFO={
      '10,5':('‡','Exactly 19 ten-point five-distance sets. Every ten-point set with at most six distances has one of the 228 candidate patterns of Szöllősi–Östergård (their count, re-derived here by an independent implementation of the orderly generation); the 19 patterns with five colours are all rigid and their real solutions (msolve) give 19 configurations up to similarity: fifteen pieces of the triangular lattice, the regular nonagon with its centre, the regular decagon, the regular 11-gon minus a vertex, and the regular pentagon with the inner pentagon of its pentagram. No published classification is known to us; the search that populated this cell earlier had found 8.',('SO20','pat')),
      '10,6':('‡','Exactly 238 ten-point six-distance sets, from the 209 six-colour patterns among the 228 candidates (all rigid; 29 patterns have more than one real realization): 184 pieces of the triangular lattice, 37 subsets of regular polygons (R₁₃ minus three vertices, R₁₂ or R₁₂ with centre minus two or three, R₁₁ with centre minus two), three twelvefold configurations found only through the candidate list (six similarity classes, all subsets of R₁₂ + R₁₂), and a pentagonal one. No published classification is known to us.',('SO20','pat')),
      '11,6':('†','Exactly 64 eleven-point six-distance sets: Amit–Astra (2026, Lemma 6) classify the eleven-point sets with at most six distances from the 62 candidate patterns (52 lattice and 10 circular configurations, counting the subsets of one regular polygon together); counted up to similarity, the 58 six-colour patterns give 49 pieces of the triangular lattice and 15 polygon subsets (R₁₃ minus two vertices in six ways, R₁₂ with centre minus two in six ways, R₁₂ minus a vertex, R₁₁ with centre minus a vertex, R₁₀ with centre). Re-verified here: the 62 patterns regenerated, certified rigid, and solved.',('AA26','SO20','pat')),
      '12,6':('‡','Exactly 14 twelve-point six-distance sets, from the 14 six-colour patterns among Szöllősi–Östergård’s 15 candidates at twelve points (all rigid, one real solution each): ten pieces of the triangular lattice (the twelve-point cores of Amit–Astra), R₁₃ minus a vertex, R₁₂, R₁₂ with centre minus a vertex, and R₁₁ with its centre. No published classification is known to us.',('SO20','AA26','pat')),
      '16,7':('‡','g(7) = 16 (Amit–Astra 2026), and exactly two sixteen-point seven-distance sets exist: Erdős–Fishburn’s two lattice pieces. Finiteness and the count are our computation (https://github.com/alonamit1729/affinemess-applets/blob/main/tools/fixed-distances/g7count/REPRODUCING.md): every sixteen-point seven-distance set contains a ten-point subset with at most six distances (a diameter-graph argument) or is R₁₄, R₁₅ or R₁₅ minus a vertex plus one or two points; the 257 ten-point sets were each extended by six points with at most two new distance values (circle intersections and circumcentres, then a clique search), the polygon cases checked directly, and the two resulting sixteen-point patterns certified rigid with a unique realization each.',('EF','AA26','SO20','pat')),
    }
    R12_RAD={0:'2',1:'2cos 15°',2:'√3',3:'√2',4:'1',5:'2sin 15°'}
    def r12_label(shells):
        parts=[]
        for k,c in shells:
            if k==6: parts.append('the centre')
            else: parts.append('%d at radius %s'%(c,R12_RAD[k]))
        return 'piece of the twelvefold set R₁₂ + R₁₂: '+', '.join(parts)
    finite_cells={}
    for n in range(6,17):
        for d in range(1,8):
            key=f'{n},{d}'
            if key in cells or d<fmin[n] or n>gmax[d]: continue
            pool=Pool()
            if key in EXACT:   # classified cells: the exact list first (polygon / lattice labels), generic-label entries after the other sources
                for e in EXACT[key]:
                    if not e['label'].startswith('distances'): pool.add(e['pts'],e['label'],'pat')
            if (n,d)==(16,7):
                for P in uniq(json.load(open('lattice_16_7.json'))):
                    rows=sorted(Counter(round(y*2/S3) for x,y in P).values(),reverse=True)
                    pool.add(P,'piece of the triangular lattice with distances 1, √3, 2, √7, 3, √12, √13 (Erdős–Fishburn, Fig. 4): rows of '+', '.join(map(str,rows))+' points','ef')
            # extension-built rigid sets
            for P in EXT.get(f'{n}{d}',[]): pool.add(P,'built by one-point extension of a smaller set','ext')
            for lab,P,N,S in EX.polygon_examples(n,d,limit=(40 if key in EXACT else 4)): pool.add(P,lab,'poly')
            for lab,P,S in EX.lattice_examples(n,d): pool.add(P,lab,'lat')
            if (n,d)==(16,7): P,_=EX.ef16(); pool.add(P,'hexagon of side 2 in the triangular lattice minus one end of each long diameter (Erdős–Fishburn)','ef')
            for e in SYM.get(key,[]): pool.add(e['pts'],struct_label(e['struct']),'sym')
            if (n,d) in ((9,5),(8,5)):
                # the one flexible 9-point <=5-distance pattern (g7fin/F56.md): T + h·R90(T), an equilateral triangle plus a perpendicular copy scaled by h
                def X3(h):
                    T=[(math.cos(2*math.pi*k/3),math.sin(2*math.pi*k/3)) for k in range(3)]; Rr=[(-y*h,x*h) for x,y in T]
                    return [(t[0]+r[0],t[1]+r[1]) for t in T for r in Rr]
                X9=X3(1.7)
                lab9='three concentric equilateral triangles with parallel sides, T + h·R₉₀(T) (Minkowski sum of an equilateral triangle with a perpendicular copy scaled by h; squared sides 1 + h² ± √3h and 1 + h², cross distances 1 and h): the only flexible nine-point pattern with at most five distances; h = 1 is Erdős–Fishburn’s nine-point four-distance set, h = √3 the lattice piece with rows 1, 3, 3, 2'
                if n==9: pool.add(X9,lab9,'tri3')
                else:
                    pool.add(X9[:4]+X9[5:],'T + h·R₉₀(T) (the flexible nine-point five-distance family) minus a vertex of the middle triangle','tri3')
                    pool.add(X9[:2]+X9[3:],'T + h·R₉₀(T) (the flexible nine-point five-distance family) minus a vertex of the outer (equivalently inner) triangle','tri3')
            n_r12=0
            for e in R12.get(key,[]):
                before=len(pool.rigid)+len(pool.fam)
                pool.add(e['pts'],r12_label(e['shells']),'r12',exact=e['exact']); n_r12+=len(pool.rigid)+len(pool.fam)-before
            if key in EXACT:
                for e in EXACT[key]:
                    if e['label'].startswith('distances'): pool.add(e['pts'],e['label'],'pat')
                assert not pool.fam and len(pool.rigid)==len(EXACT[key]),(key,len(pool.rigid),len(EXACT[key]),len(pool.fam))
                for P,lab,src in pool.rigid: assert any(B.same(P,e['pts']) for e in EXACT[key]),(key,lab)
            if key=='16,7':
                assert not pool.fam and len(pool.rigid)==2
            if key in EXACT_INFO:
                mark,note,srcs=EXACT_INFO[key]
                items=[item(P,lab,'rigid',0) for P,lab,src in pool.rigid]
                cell(n,d,{0:E(len(items))},'exact',marks=mark,note=note,finite=True)
                cat(n,d,items,True,f'The {len(items)} {n}-point {d}-distance sets' if len(items)>2 else 'The two sixteen-point seven-distance sets','All configurations shown are rigid; the list is complete.',sources=srcs)
                continue
            c=pool.counts()
            fin=(n,d) in finite_cells
            dims={k:X(v) for k,v in c.items()}
            if fin:
                dims={0:X(c.get(0,0))}
            e=2*n-4-math.comb(n,2)+d
            note,srcs=finite_cells.get((n,d),('',()))
            if (n,d)==(16,7): srcs=('EF','AA26','SO20')
            if not fin:
                note='Exact counts unknown; the bounds come from the configurations shown (each verified numerically to have the stated local dimension). '
                if e>=1: note+=f'Every component has dimension ≥ {e}, so there are no rigid configurations. '
                else: note+='Finiteness is not known to us, so higher-dimensional families cannot be excluded. '
                if key in R12: note+=f'All {len(R12[key])} subsets with these parameters of the twelvefold set R₁₂ + R₁₂ — the 73 sums of two twelfth roots of unity: the centre and points at angles that are multiples of 15° on circles of radii 2, 2cos 15°, √3, √2, 1 and 2sin 15° (circumradius 1 for the 12-gon) — are included, from an exhaustive sweep up to similarity; {n_r12} of them were not previously in the list. '
                if (n,d)==(16,7): note+='Sixteen points is the maximum for seven distances: Amit and Astra (2026) prove g(7) = 16 by a computer-assisted diameter-graph argument whose one computational input taken from the literature is the Szöllősi–Östergård exhaustive candidate list of eleven-point patterns with at most six distances. The two Erdős–Fishburn lattice sets are the only known sixteen-point seven-distance sets; the extremal sets are not classified, so the count is a lower bound. '
            items=[item(P,lab,'rigid',0) for P,lab,src in pool.rigid]+[item(G,lab,'family',dim) for G,lab,src,h,dim in pool.fam]
            for idx,ex in pool.exact.items():
                it=items[idx]
                if it['exact'] is None or any(x is None for x in it['exact']):
                    P=pool.rigid[idx][0]
                    it['exact']=ex; it['frame']='twelvefold frame: the sum set R₁₂ + R₁₂ of the regular 12-gon of circumradius 1, centre at the origin, a vertex of the 12-gon at angle 0'
                    it['coords']=[[('%.10g'%x) if abs(x)>5e-11 else '0',('%.10g'%y) if abs(y)>5e-11 else '0'] for x,y in P]
            cell(n,d,dims,'bounds',note=note,finite=True if fin else None)
            cat(n,d,items,False,f'{n} points, {d} distances: known configurations',listing_note(len(pool.rigid),len(pool.fam),False),sources=srcs)
    for key,c in cats.items():
        dd=cells[key]['d']
        for it in c['items']: assert it['d']==dd,(key,it['label'],it['d'])
        c['items'].sort(key=lambda it: it['dim'])   # stable: rigid first, then 1d, 2d, ... ; original order within a dimension
    meta={'rows':[str(k) for k in range(3,17)],'cols':list(range(1,8))}
    add_stratification(cells,cats)
    json.dump({'cells':cells,'cats':cats,'meta':meta},open('applet_data.json','w'))
    print('cells',len(cells),'items',sum(len(c['items']) for c in cats.values()))
    for k,c in cells.items(): print(k,c['status'],c['dims'],c['marks'])
if __name__=='__main__': main()
