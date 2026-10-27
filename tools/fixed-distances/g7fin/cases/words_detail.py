#!/usr/bin/env python3
import sys, json, itertools
N=int(sys.argv[1]); sys.argv=[sys.argv[0],str(N)]
src=open('words_odd.py').read().split("words=sorted(")[0]
exec(src)
L='abcde'
def all_assignments(w):
    ch=chords(w); types=sorted(set(ch.values())); idx={t:i for i,t in enumerate(types)}
    dom=[list(WIN[t[0]]) for t in types]; lt=set()
    for t1 in types:
        for t2 in types:
            if t1[0]==t2[0] and less(t1[0],t1[1],t2[1]): lt.add((idx[t1],idx[t2]))
    for j in range(N):
        for sgn in (1,-1):
            seq=[idx[ch[frozenset((j,(j+sgn*2*s)%N))]] for s in range(1,q+1)]
            for a,b in zip(seq,seq[1:]): lt.add((a,b))
    sols=[a for a in itertools.product(*dom) if all(a[x]<a[y] for x,y in lt)]
    return types,sols
for w in json.load(open(f'words{N}_survivors.json')):
    if len(set(w))<3: continue
    types,sols=all_assignments(w)
    # coincidences common to all assignments: pairs of types always in the same class
    common=[(types[i],types[j]) for i in range(len(types)) for j in range(i+1,len(types)) if all(s[i]==s[j] for s in sols)]
    print(''.join(L[x] for x in w), "assignments:",len(sols),"forced coincidences:",len(common),
          [(f"{t1[0]}:{''.join(L[x] for x in t1[1])}",f"{t2[0]}:{''.join(L[x] for x in t2[1])}") for t1,t2 in common])
