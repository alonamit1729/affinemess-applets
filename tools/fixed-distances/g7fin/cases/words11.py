#!/usr/bin/env python3
"""Case (a): X_D is a convex 11-gon p_0..p_10 with diameter graph the step-5 star, and
exactly 7 distance classes e_1<...<e_6<e_7=D.  The star s_j=p_{5j} is an equilateral
closed 11-gon whose angle theta_j is determined by the convex side it subtends, so the
side-length word (over alpha<beta<gamma, i.e. 0<1<2) is a cyclic word of length 11.
Chord dictionary (star r-chord <-> convex step): r=2 side (window e1..e3),
r=3 step 4 (e4..e6), r=4 step 2 (e2..e4), r=5 step 3 (e3..e5).
Convexity chains: for each j, c(j,j+2)<c(j,j+4)<c(j,j+6)<c(j,j+8), both directions.
Monotonicity (proved for all angles < pi/3):
  r=2: increasing in the angle; r=3: f strictly decreasing in both angles;
  r=4: g strictly increasing in the two outer angles, strictly decreasing in the middle one;
  r=5: no inequality used (only 'same type => same class').
Count the cyclic 3-letter words (all three letters present, up to dihedral symmetry)
for which a class assignment exists."""
import itertools, sys
N=11
def canon(w):
    return min(tuple(w[i:]+w[:i]) for w in (w,w[::-1]) for i in range(N))
words=sorted({canon(list(w)) for w in itertools.product(range(3),repeat=N) if len(set(w))==3})
print("3-letter cyclic words up to dihedral symmetry:",len(words))

def chords(w):
    """return dict chord(frozenset{j,k}) -> (r, type) with r in 2..5"""
    ch={}
    for j in range(N):
        for r in (2,3,4,5):
            k=(j+r)%N
            ang=tuple(w[(j+t)%N] for t in range(1,r))
            typ=min(ang,ang[::-1])
            ch[frozenset((j,k))]=(r,typ)
    assert len(ch)==44
    return ch
WIN={2:(1,2,3),3:(4,5,6),4:(2,3,4),5:(3,4,5)}
def less(r,t1,t2):
    """True if chord type t1 is provably strictly shorter than t2 (same r)"""
    if t1==t2: return False
    if r==2: return t1[0]<t2[0]
    if r==3:   # decreasing in each: t1 shorter iff t1 >= t2 componentwise (sorted)
        a,b=sorted(t1); c,d=sorted(t2); return a>=c and b>=d
    if r==4:   # increasing outer, decreasing middle, up to reversal
        for u in (t1,t1[::-1]):
            if u[0]<=t2[0] and u[2]<=t2[2] and u[1]>=t2[1]: return True
        return False
    return False
def feasible(w):
    ch=chords(w)
    types=sorted(set(ch.values()))
    idx={t:i for i,t in enumerate(types)}
    dom=[list(WIN[t[0]]) for t in types]
    # constraints: strict less between type variables
    lt=set()
    for t1 in types:
        for t2 in types:
            if t1[0]==t2[0] and less(t1[0],t1[1],t2[1]): lt.add((idx[t1],idx[t2]))
    for j in range(N):
        for sgn in (1,-1):
            seq=[idx[ch[frozenset((j,(j+sgn*s)%N))]] for s in (2,4,6,8)]
            for a,b in zip(seq,seq[1:]): lt.add((a,b))
    n=len(types); sols=[]
    order=sorted(range(n),key=lambda i:len(dom[i]))
    assign=[None]*n
    def bt(pos):
        if pos==n:
            if set(assign)=={1,2,3,4,5,6}: sols.append(tuple(assign)); return True
            return False
        i=order[pos]
        for v in dom[i]:
            ok=True
            for (a,b) in lt:
                if a==i and assign[b] is not None and not v<assign[b]: ok=False;break
                if b==i and assign[a] is not None and not assign[a]<v: ok=False;break
            if ok:
                assign[i]=v
                if bt(pos+1): return True
                assign[i]=None
        return False
    return bt(0), types
surv=[]
for w in words:
    ok,types=feasible(w)
    if ok: surv.append(w)
print("words admitting a 7-class assignment:",len(surv))
from collections import Counter
print("by letter counts (n_alpha,n_beta,n_gamma):",sorted(Counter(tuple(w.count(l) for l in range(3)) for w in surv).items()))
print("palindromic (mirror-symmetric) survivors:",sum(1 for w in surv if any((w[i:]+w[:i])==(w[i:]+w[:i])[::-1] for i in range(N))))
import json; json.dump([list(w) for w in surv],open('words11_survivors.json','w'))
print("examples:",surv[:8])
