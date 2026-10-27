#!/usr/bin/env python3
"""Generalisation of words11.py to an odd diameter cycle C_N (step k=(N-1)/2 star) whose
vertex set Y has at most 7 distance classes e_1<..<e_6<D.  Convex chain of q=(N-3)/2
non-diameter steps: step t in {e_t..e_{t+6-q}}; sides use letters 0..(6-q).
Star r-chord <-> convex step k*r mod N (folded).  Monotonicity used for r=2,3,4 only.
Option --all7: require all 6 non-diameter classes present (only valid when |Y|>=11).
Reports, by number of letters, the words admitting a class assignment."""
import itertools, sys, json
from collections import Counter
N=int(sys.argv[1]); ALL7='--all7' in sys.argv
k=(N-1)//2; q=(N-3)//2; LETTERS=7-q
kinv=pow(k,-1,N)
def canon(w): return min(tuple(w[i:]+w[:i]) for w in (w,w[::-1]) for i in range(N))
def fold(t): t%=N; return min(t,N-t)
RS=[r for r in range(2,N) if fold(k*r)>=1 and r<=N//2+1]
# each convex step t in 1..q corresponds to exactly one r in 2..N-2 up to the complementary chain; use r<=N-r
RS=sorted({min(r,N-r) for r in range(2,N-1)})
STEP={r:fold(k*r) for r in RS}
WIN={r:tuple(range(STEP[r],STEP[r]+7-q)) for r in RS}
def chords(w):
    ch={}
    for j in range(N):
        for r in RS:
            kk=(j+r)%N; ang=tuple(w[(j+t)%N] for t in range(1,r))
            ch.setdefault(frozenset((j,kk)),(r,min(ang,ang[::-1])))
    assert len(ch)==N*(N-1)//2-N, (len(ch))
    return ch
def less(r,t1,t2):
    if t1==t2: return False
    if r==2: return t1[0]<t2[0]
    if r==3:
        a,b=sorted(t1); c,d=sorted(t2); return a>=c and b>=d
    if r==4:
        for u in (t1,t1[::-1]):
            if u[0]<=t2[0] and u[2]<=t2[2] and u[1]>=t2[1]: return True
    return False
def feasible(w):
    ch=chords(w); types=sorted(set(ch.values())); idx={t:i for i,t in enumerate(types)}
    dom=[list(WIN[t[0]]) for t in types]; lt=set()
    for t1 in types:
        for t2 in types:
            if t1[0]==t2[0] and less(t1[0],t1[1],t2[1]): lt.add((idx[t1],idx[t2]))
    for j in range(N):
        for sgn in (1,-1):
            seq=[idx[ch[frozenset((j,(j+sgn*2*s)%N))]] for s in range(1,q+1)]
            for a,b in zip(seq,seq[1:]): lt.add((a,b))
    n=len(types); order=sorted(range(n),key=lambda i:len(dom[i])); assign=[None]*n
    def bt(pos):
        if pos==n:
            return (not ALL7) or set(assign)=={1,2,3,4,5,6}
        i=order[pos]
        for v in dom[i]:
            if all((not(a==i and assign[b] is not None and not v<assign[b])) and (not(b==i and assign[a] is not None and not assign[a]<v)) for a,b in lt):
                assign[i]=v
                if bt(pos+1): return True
                assign[i]=None
        return False
    return bt(0)
words=sorted({canon(list(w)) for w in itertools.product(range(LETTERS),repeat=N)})
# identify words up to relabelling of letters that preserves order: a word uses letters 0..L-1 consecutively
words=[w for w in words if sorted(set(w))==list(range(len(set(w))))]
print(f"N={N}: star step {k}, chord dictionary r->convex step {STEP}, windows {WIN}, letters<= {LETTERS}")
cnt=Counter(); surv=[]
for w in words:
    if feasible(w): cnt[len(set(w))]+=1; surv.append(w)
tot=Counter(len(set(w)) for w in words)
for L in sorted(tot): print(f"  {L} letters: {cnt[L]} of {tot[L]} cyclic words admit an assignment")
json.dump([list(w) for w in surv],open(f'words{N}_survivors{"_all7" if ALL7 else ""}.json','w'))
