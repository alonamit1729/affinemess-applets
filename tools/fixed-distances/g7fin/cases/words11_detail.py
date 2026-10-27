#!/usr/bin/env python3
"""List every class assignment for the surviving words of words11.py."""
import json, itertools
exec(open('words11.py').read().split("surv=[]")[0].replace('print("3-letter','#print("3-letter'))
def all_assignments(w):
    ch=chords(w); types=sorted(set(ch.values())); idx={t:i for i,t in enumerate(types)}
    dom=[list(WIN[t[0]]) for t in types]; lt=set()
    for t1 in types:
        for t2 in types:
            if t1[0]==t2[0] and less(t1[0],t1[1],t2[1]): lt.add((idx[t1],idx[t2]))
    for j in range(N):
        for sgn in (1,-1):
            seq=[idx[ch[frozenset((j,(j+sgn*s)%N))]] for s in (2,4,6,8)]
            for a,b in zip(seq,seq[1:]): lt.add((a,b))
    sols=[]
    for assign in itertools.product(*dom):
        if set(assign)!={1,2,3,4,5,6}: continue
        if all(assign[a]<assign[b] for a,b in lt): sols.append(assign)
    return types,sols
L='abc'
for w in json.load(open('words11_survivors.json')):
    print("word",''.join(L[x] for x in w))
    types,sols=all_assignments(w)
    print("  chord types present:",[(r,''.join(L[x] for x in t)) for r,t in types])
    print("  number of class assignments:",len(sols))
    for s in sols:
        byclass={}
        for t,c in zip(types,s): byclass.setdefault(c,[]).append((t[0],''.join(L[x] for x in t[1])))
        print("   ", {c:byclass[c] for c in sorted(byclass)})
