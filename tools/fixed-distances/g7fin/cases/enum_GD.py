#!/usr/bin/env python3
"""Enumerate combinatorial candidates for the diameter graph G_D on m diameter
endpoints, using:
 (W)  Wei: no even cycle, at most one cycle (so each component is a tree or unicyclic
      with an odd cycle); no isolated vertex.
 (L)  Leaf lemma (proved in REPORT.md): if deg v >= 3, every neighbour of v other than
      the two angularly extreme ones is a leaf.  Hence deleting all leaves leaves a
      graph of maximum degree <= 2: every component is a caterpillar, or an odd cycle
      with pendant leaves ('odd sun'), at most one cycle in all.
Filter by independence number alpha.  Output counts and degree data.
"""
import networkx as nx, itertools, sys
from networkx.algorithms import isomorphism

def is_caterpillar(T):
    leaves=[v for v in T if T.degree(v)==1]
    H=T.copy(); H.remove_nodes_from(leaves)
    return H.number_of_nodes()==0 or (nx.is_connected(H) and max(dict(H.degree()).values(),default=0)<=2 and nx.is_tree(H))

def caterpillars(k):
    if k==1: return []          # no isolated vertices
    if k==2: return [nx.path_graph(2)]
    return [T for T in nx.nonisomorphic_trees(k) if is_caterpillar(T)]

def suns(k):
    """odd cycle C_c plus k-c pendant leaves distributed on cycle vertices, up to iso"""
    out=[]
    for c in range(3,k+1,2):
        extra=k-c
        seen=[]
        for comp in itertools.product(range(extra+1),repeat=c):
            if sum(comp)!=extra: continue
            G=nx.cycle_graph(c); n=c
            for i,l in enumerate(comp):
                for _ in range(l): G.add_edge(i,n); n+=1
            if not any(nx.is_isomorphic(G,H) for H in seen): seen.append(G)
        out+=seen
    return out

def partitions(m,maxpart=None):
    if maxpart is None: maxpart=m
    if m==0: yield []; return
    for p in range(min(m,maxpart),1,-1):   # parts >=2
        for rest in partitions(m-p,p): yield [p]+rest

def alpha(G):
    return len(nx.max_weight_clique(nx.complement(G),weight=None)[0])

def dedup(graphs):
    buckets={}
    out=[]
    for G in graphs:
        h=nx.weisfeiler_lehman_graph_hash(G)
        if not any(nx.is_isomorphic(G,H) for H in buckets.get(h,[])):
            buckets.setdefault(h,[]).append(G); out.append(G)
    return out

def candidates(m):
    cat={k:caterpillars(k) for k in range(2,m+1)}
    sun={k:suns(k) for k in range(3,m+1)}
    res=[]
    for part in partitions(m):
        # choose which part (if any) is the sun
        options=[None]+[i for i in range(len(part)) if part[i]>=3]
        for si in options:
            lists=[]
            for i,p in enumerate(part):
                lists.append(sun[p] if i==si else cat[p])
            for combo in itertools.product(*lists):
                res.append(nx.disjoint_union_all(combo))
    return dedup(res)

if __name__=='__main__':
    m=int(sys.argv[1]); target=int(sys.argv[2])
    C=candidates(m)
    good=[G for G in C if alpha(G)==target]
    print(f"m={m}: {len(C)} graphs satisfy (W)+(L); {len(good)} have alpha={target}")
    forests=[G for G in good if nx.is_forest(G)]
    print(f"  forests: {len(forests)}, unicyclic: {len(good)-len(forests)}")
    from collections import Counter
    print("  by (cycle length, #components, max degree):")
    cnt=Counter()
    for G in good:
        cyc=[len(c) for c in nx.cycle_basis(G)]
        cnt[(cyc[0] if cyc else 0, nx.number_connected_components(G), max(dict(G.degree()).values()))]+=1
    for k,v in sorted(cnt.items()): print("   ",k,v)
    print("  number of maximum independent sets per graph (min/max):",
          min(sum(1 for S in itertools.combinations(G.nodes,target) if not any(G.has_edge(a,b) for a,b in itertools.combinations(S,2))) for G in good),
          max(sum(1 for S in itertools.combinations(G.nodes,target) if not any(G.has_edge(a,b) for a,b in itertools.combinations(S,2))) for G in good))
    import json
    json.dump([sorted(map(list,G.edges())) for G in good],open(f'GD_m{m}_alpha{target}.json','w'))
