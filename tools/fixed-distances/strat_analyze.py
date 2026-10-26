"""Shape-space incidence graph from strat_{d}.json (geometric arcs per type).  Usage: python3 strat_analyze.py d
Vertices: rigid (5,<d)-sets (identified in the applet's certified catalogs); edges: geometric arcs of the one-parameter
families; half-open edges end at degenerate configurations (two points coincide) or at infinity."""
import json, sys, math, itertools, collections
import numpy as np, networkx as nx
d=int(sys.argv[1]); n=int(sys.argv[2]) if len(sys.argv)>2 else 5
S=json.load(open(f'strat_{n}_{d}.json' if n!=5 else f'strat_{d}.json'))
def dmat(P): return np.array([[math.dist(a,b) for b in P] for a in P])
PERMS=list(itertools.permutations(range(n)))
def canon(P):
    D=dmat(P); m=D[D>1e-9].min(); D=D/m; best=None
    for p in PERMS:
        t=tuple(round(D[p[i],p[j]],6) for i in range(n) for j in range(i+1,n))
        if best is None or t<best: best=t
    return best
def close(c1,c2,tol=2e-5): return max(abs(a-b) for a,b in zip(c1,c2))<tol
A=json.load(open('applet_data.json'))
ref=[]
for dd in range(2,d):
    key=f'{n},{dd}'
    if key in A['cats']:
        for k,it in enumerate(A['cats'][key]['items']):
            if it.get('kind','rigid')=='rigid': ref.append((canon([[float(x) for x in c] for c in it['coords']]),(dd,k,it.get('label',''))))
def lookup(c):
    for rc,lab in ref:
        if close(rc,c): return lab
    return None
V=[]   # list of (canon, label)
def vid(c,tol=2e-5):
    for i,(k,_) in enumerate(V):
        if close(k,c,tol): return i
    V.append((c,lookup(c))); return len(V)-1
Edges=[]; stats=collections.Counter()
for rec in S:
    for g in rec['arcs']:
        ends=[]
        for e in g['ends']:
            if e.get('end')=='fold': ends.append(('v',vid(tuple(e['canon']),5e-3)))
            elif 'end' in e: ends.append(('end',e['end']))
            else: ends.append(('v',vid(tuple(e['canon']))))
        Edges.append({'type':rec['type'],'ends':ends,'mult':g['mult'],'inv':(g['ratio_max'],g['ratio_min'])})
        kinds=tuple(sorted(e[0] if e[0]=='v' else e[1] for e in ends)); stats[kinds]+=1
print('types:',len(S),' arcs by end kinds:',dict(stats))
unid=[i for i,(c,l) in enumerate(V) if l is None]
print('vertices:',len(V),' unidentified:',len(unid))
# classify unidentified vertices (points on lower-dimensional families, not rigid sets) by colouring type
if unid:
    E=list(itertools.combinations(range(n),2))
    perms=[[E.index(tuple(sorted((p[a],p[b])))) for a,b in E] for p in itertools.permutations(range(n))]
    def ccol(rgs):
        best=None
        for pe in perms:
            c=[rgs[pe[i]] for i in range(len(E))]; mp={}; out=[]
            for x in c:
                if x not in mp: mp[x]=len(mp)
                out.append(mp[x])
            t=tuple(out)
            if best is None or t<best: best=t
        return best
    def coltype(P,tol=1e-6):
        ds=[math.dist(P[i],P[j]) for i,j in E]; mx=max(ds); vals=[]; col=[]
        for x in ds:
            for k,v in enumerate(vals):
                if abs(x-v)<tol*mx: col.append(k); break
            else: vals.append(x); col.append(len(vals)-1)
        return ccol(col),len(vals)
    famtypes={}
    for dd in range(2,d):
        try: R=json.load(open(f'cert{n}{dd}.json'))
        except FileNotFoundError: continue
        for r in R:
            if any((q.get('real') or q.get('real_slices')) and q['dim']>0 for q in r['primes']) or r.get('unresolved'): famtypes[(dd,ccol(tuple(r['col'])))]=r['idx']
    rootpts={}; folds=set()
    for rec in S:
        for q in rec['roots']:
            if not q['coincident']: rootpts.setdefault(tuple(q['canon']),q['v'])
        for g in rec['arcs']:
            for e in g['ends']:
                if e.get('end')=='fold': rootpts.setdefault(tuple(e['canon']),e['v']); folds.add(tuple(e['canon']))
    kinds=collections.Counter()
    for i in unid:
        c=V[i][0]; v=None
        for k,vv in rootpts.items():
            if close(k,c): v=vv; break
        if v is None: kinds['no point']+=1; continue
        P=[(0,0),(v[0],0)]+[(v[2*k-3],v[2*k-2]) for k in range(2,n)]
        ct,nv=coltype(P); t=famtypes.get((nv,ct))
        if any(close(k,c) for k in folds): V[i]=(c,(nv,None,f'fold point: a mirror-symmetric {nv}-distance configuration where the family meets its own reflection')); kinds['fold']+=1; continue
        if t is not None: V[i]=(c,(nv,None,f'interior point of a {nv}-distance family (type {t})')); kinds[f'on (5,{nv}) family']+=1
        else: kinds[f'unknown d={nv}']+=1
    print('unidentified vertices classified:',dict(kinds))
print('vertices by cell:',collections.Counter(l[0] for c,l in V if l))
val=collections.Counter()
for e in Edges:
    for end in e['ends']:
        if end[0]=='v': val[end[1]]+=1
print('valence histogram:',sorted(collections.Counter(val.values()).items()))
reached={l[:2] for c,l in V if l}
missing=[lab for c,lab in ref if lab[0]==d-1 and lab[:2] not in reached]
print(f'({n},{d-1}) rigid sets not reached by any (5,{d}) family: {len(missing)}',[m[1] for m in missing][:50])
G=nx.MultiGraph(); leaf=[0]
for i in range(len(V)): G.add_node(('v',i))
for e in Edges:
    a,b=e['ends']
    if a[0]=='end': leaf[0]+=1; a=('leaf',leaf[0])
    if b[0]=='end': leaf[0]+=1; b=('leaf',leaf[0])
    G.add_edge(a,b)
nV=G.number_of_nodes(); nE=G.number_of_edges(); c=nx.number_connected_components(G)
print('graph with ends as leaves: nodes',nV,'edges',nE,'components',c,'b1 =',nE-nV+c)
H=G.copy(); H.remove_nodes_from([v for v in list(H.nodes) if v[0]=='leaf'])
cH=nx.number_connected_components(H)
print('vertex subgraph (rigid sets only): nodes',H.number_of_nodes(),'edges',H.number_of_edges(),'components',cH,'b1 =',H.number_of_edges()-H.number_of_nodes()+cH)
print('mult histogram of arcs:',sorted(collections.Counter(e['mult'] for e in Edges).items()))
json.dump({'vertices':[{'id':i,'canon':list(c),'cell':l[0] if l else None,'item':l[1] if l else None,'label':l[2] if l else None} for i,(c,l) in enumerate(V)],
           'edges':[{'type':e['type'],'ends':[list(x) for x in e['ends']],'mult':e['mult'],'inv':list(e['inv'])} for e in Edges],
           'types':{str(rec['type']):rec['col'] for rec in S}},open(f'strat_graph_{n}_{d}.json' if n!=5 else f'strat_graph_{d}.json','w'))
