#!/usr/bin/env python3
"""Rigorous 2-D interval certificate (mpmath.iv, outward-rounded interval arithmetic)
that the chord coincidences forced on each surviving word of words11.py have NO common
solution on the closed admissible domain
      0 <= a <= b <= c <= pi/3,   u a + v b + s c = pi   (u,v,s = letter counts).
Parametrise by (a,b), c = (pi - u a - v b)/s.  A box is discarded if it provably misses
the domain, accepted if some residual interval excludes 0, else bisected.
Squared chord of r unit edges with interior angles th_1..th_{r-1}:
      |sum_k e_k|^2 = sum_{j,k} cos(phi_j - phi_k),  phi_k = k pi - sum_{i<k} th_i.
"""
from mpmath import iv, mp
import sys, json
iv.prec=60
PI=iv.pi
def chord(*th):
    r=len(th)+1
    phi=[k*PI-sum(th[:k],iv.mpf(0)) for k in range(r)]
    return sum((iv.cos(phi[j]-phi[k]) for j in range(r) for k in range(r)),iv.mpf(0))
SYS={
 'aaaaaaaabcb':((8,2,1),lambda a,b,c:[chord(b)-chord(a,a,a), chord(c)-chord(a,a,b), chord(a,a,b)-chord(a,a,b,c),
   chord(a,a,b,c)-chord(a,b,c,b), chord(b,c)-chord(a,b,c), chord(a,b,c)-chord(a,a,a,b), chord(a,b)-chord(a,a,a,a)]),
 'ababababccb':((4,5,2),lambda a,b,c:[chord(b)-chord(a,b,a), chord(c)-chord(b,a,b), chord(b,a,b)-chord(a,b,c,c),
   chord(c,c)-chord(a,b,c), chord(a,b,c)-chord(b,c,c), chord(b,c,c)-chord(b,a,b,c), chord(b,a,b,c)-chord(b,c,c,b),
   chord(b,c)-chord(a,b,a,b)]),
}
def excludes0(x): return x.a>0 or x.b<0
summary={}
for w,((u,v,s),res) in SYS.items():
    T=PI/3
    stack=[(iv.mpf([0,T.b]),iv.mpf([0,T.b]),0)]
    accepted=0; discarded=0; maxdepth=0; used=[0]*len(res(iv.mpf(0),iv.mpf(0),iv.mpf(0)))
    while stack:
        A,B,d=stack.pop(); maxdepth=max(maxdepth,d)
        C=(PI-u*A-v*B)/s
        # closed-domain test: a<=b, b<=c, c<=pi/3 must be satisfiable in the box
        if B.b<A.a or C.b<B.a or C.a>T.b: discarded+=1; continue
        R=res(A,B,C)
        hit=[i for i,x in enumerate(R) if excludes0(x)]
        if hit: accepted+=1; used[hit[0]]+=1; continue
        if d>=40: raise RuntimeError(('inconclusive',w,A,B,d))
        if (A.b-A.a)>=(B.b-B.a):
            m=(A.a+A.b)/2; stack+= [(iv.mpf([A.a,m]),B,d+1),(iv.mpf([m,A.b]),B,d+1)]
        else:
            m=(B.a+B.b)/2; stack+= [(A,iv.mpf([B.a,m]),d+1),(A,iv.mpf([m,B.b]),d+1)]
    summary[w]=dict(accepted=accepted,discarded=discarded,max_depth=maxdepth,equation_used=used)
    print("PROVED word",w,": no common zero on the closed domain;",accepted,"boxes accepted,",discarded,"outside domain, max depth",maxdepth,"; which equation excluded:",used,flush=True)
json.dump(summary,open('words11_certificate.json','w'),indent=1)
