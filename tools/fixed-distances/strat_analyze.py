"""Shape-space incidence graph from strat_{d}.json (geometric arcs per type).  Usage: python3 strat_analyze.py d
Vertices: rigid (5,<d)-sets (identified in the applet's certified catalogs); edges: geometric arcs of the one-parameter
families; half-open edges end at degenerate configurations (two points coincide) or at infinity."""
import json, sys, math, itertools, collections
import numpy as np, networkx as nx
d=int(sys.argv[1]); n=5
S=json.load(open(f'strat_{d}.json'))
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
    key=f'5,{dd}'
    if key in A['cats']:
        for k,it in enumerate(A['cats'][key]['items']):
            if it.get('kind','rigid')=='rigid': ref.append((canon([[float(x) for x in c] for c in it['coords']]),(dd,k,it.get('label',''))))
def lookup(c):
    for rc,lab in ref:
        if close(rc,c): return lab
    return None
V=[]   # list of (canon, label)
def vid(c):
    for i,(k,_) in enumerate(V):
        if close(k,c): return i
    V.append((c,lookup(c))); return len(V)-1
Edges=[]; stats=collections.Counter()
for rec in S:
    for g in rec['arcs']:
        ends=[]
        for e in g['ends']:
            if 'end' in e: ends.append(('end',e['end']))
            else: ends.append(('v',vid(tuple(e['canon']))))
        Edges.append({'type':rec['type'],'ends':ends,'mult':g['mult'],'inv':(g['ratio_max'],g['ratio_min'])})
        kinds=tuple(sorted(e[0] if e[0]=='v' else e[1] for e in ends)); stats[kinds]+=1
print('types:',len(S),' arcs by end kinds:',dict(stats))
unid=[i for i,(c,l) in enumerate(V) if l is None]
print('vertices:',len(V),' unidentified:',len(unid))
print('vertices by cell:',collections.Counter(l[0] for c,l in V if l))
val=collections.Counter()
for e in Edges:
    for end in e['ends']:
        if end[0]=='v': val[end[1]]+=1
print('valence histogram:',sorted(collections.Counter(val.values()).items()))
reached={l[:2] for c,l in V if l}
missing=[lab for c,lab in ref if lab[0]==d-1 and lab[:2] not in reached]
print(f'(5,{d-1}) rigid sets not reached by any (5,{d}) family: {len(missing)}',[m[1] for m in missing][:50])
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
           'types':{str(rec['type']):rec['col'] for rec in S}},open(f'strat_graph_{d}.json','w'))
