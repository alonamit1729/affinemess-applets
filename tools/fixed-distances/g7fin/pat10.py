"""Patterns (edge colourings of K_n up to vertex and colour permutations) of all known n-point sets with at most 6 distances,
for n = 10 (and 9, 11 for comparison), from: the 13-point sets R13, R12+O, the 11 twelve-point lattice cores, R11+O, R10+O,
and the applet's catalogs.  Counts distinct pattern classes via a Weisfeiler-Lehman style hash, with VF2 confirmation of merges."""
import json, math, itertools, collections, hashlib, sys
import numpy as np, networkx as nx
from networkx.algorithms import isomorphism
n=int(sys.argv[1]) if len(sys.argv)>1 else 10
S3=math.sqrt(3)
def poly(k,ids): return [(math.cos(2*math.pi*i/k),math.sin(2*math.pi*i/k)) for i in ids]
sets=[]
sets.append(('R13',poly(13,range(13))))
sets.append(('R12+O',[(0,0)]+poly(12,range(12))))
sets.append(('R11+O',[(0,0)]+poly(11,range(11))))
sets.append(('R10+O',[(0,0)]+poly(10,range(10))))
cores=json.load(open('g7/cores.json'))
for i,L in enumerate(cores): sets.append((f'L{i}',[(a+b/2,b*S3/2) for a,b in L]))
import os
if os.path.exists(f'g7/latt_tri_4_25_{n}.json'):
    for i,P in enumerate(json.load(open(f'g7/latt_tri_4_25_{n}.json'))[str(n)]):
        sets.append((f'lat{n}#{i}',[(a+b/2,b*S3/2) for a,b in P]))
if n==10 and os.path.exists('g7fin/pent_ends.json'):
    for i,o in enumerate(json.load(open('g7fin/pent_ends.json'))): sets.append((f'pent{o["family"]}#{i}',[tuple(p) for p in o['pts']]))
A=json.load(open('applet_data.json'))
for key in (f'{n},5',f'{n},6',f'{n},4'):
    if key in A['cats']:
        for k,it in enumerate(A['cats'][key]['items']):
            if it['kind']=='rigid': sets.append((f'cat{key}#{k+1}',[(float(x),float(y)) for x,y in it['coords']]))
def colouring(P,tol=1e-7):
    m=len(P); ds={}
    vals=[]; col={}
    for i,j in itertools.combinations(range(m),2):
        x=math.dist(P[i],P[j])
        for c,v in enumerate(vals):
            if abs(x-v)<tol*max(1,x): col[(i,j)]=c; break
        else: vals.append(x); col[(i,j)]=len(vals)-1
    return col,len(vals)
def wl_hash(col,m,rounds=4):
    # colour-class sizes are a colour invariant; use them to canonicalise colours partially
    csize=collections.Counter(col.values())
    lab={i:'v' for i in range(m)}
    for _ in range(rounds):
        new={}
        for i in range(m):
            prof=sorted((csize[col[tuple(sorted((i,j)))]],lab[j]) for j in range(m) if j!=i)
            new[i]=hashlib.md5((lab[i]+str(prof)).encode()).hexdigest()[:12]
        lab=new
    # also include the edge-colour multiset of pair-labels per colour class size
    edges=sorted((csize[c],tuple(sorted((lab[i],lab[j])))) for (i,j),c in col.items())
    return hashlib.md5(str((sorted(lab.values()),edges)).encode()).hexdigest()
def iso(colA,colB,m):
    GA=nx.complete_graph(m); GB=nx.complete_graph(m)
    for (i,j),c in colA.items(): GA[i][j]['c']=c
    for (i,j),c in colB.items(): GB[i][j]['c']=c
    kA=max(colA.values())+1; kB=max(colB.values())+1
    if kA!=kB: return False
    for perm in itertools.permutations(range(kA)):
        GM=isomorphism.GraphMatcher(GA,GB,edge_match=lambda e1,e2: perm[e1['c']]==e2['c'])
        if GM.is_isomorphic(): return True
    return False
classes={}   # hash -> list of (name, col)
sources=collections.defaultdict(set)
for name,P in sets:
    if len(P)<n: continue
    for sub in itertools.combinations(range(len(P)),n):
        Q=[P[i] for i in sub]; col,k=colouring(Q)
        if k>6: continue
        h=wl_hash(col,n)
        reps=classes.setdefault(h,[])
        for r in reps:
            if iso(r[1],col,n): r[3].add(name.split('#')[0]); break
        else: reps.append((name,col,k,{name.split('#')[0]}))
allc=[r for reps in classes.values() for r in reps]
print(f'n={n}: {len(allc)} pattern classes ({len(classes)} hashes) from {len(sets)} source sets')
byk=collections.Counter(v[2] for v in allc); print('by number of distances:',dict(sorted(byk.items())))
bysrc=collections.Counter()
for r in allc:
    s=r[3]
    if all(x.startswith('cat') for x in s): bysrc['catalog only']+=1
    elif all(x.startswith('lat') for x in s): bysrc['lattice enumeration only']+=1
    elif all(x.startswith('pent') for x in s): bysrc['pentagon endpoints only']+=1
    elif all(x.startswith('L') for x in s): bysrc['lattice cores only']+=1
    elif all(x.startswith('R') for x in s): bysrc['circular only']+=1
    else: bysrc['mixed']+=1
print(dict(bysrc))
json.dump([{'name':v[0],'k':v[2],'col':{f'{i},{j}':c for (i,j),c in v[1].items()},'sources':sorted(v[3])} for v in allc],open(f'g7fin/patterns_{n}.json','w'))
