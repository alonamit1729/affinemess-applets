#!/usr/bin/env python3
"""For odd N and every nonconstant binary word b, check exactly that the closure
polynomial P_b(z)=sum_k (-1)^k omega^{B_k} z^{nu k - N B_k} (coefficients in
Z[omega]=Z[T]/Phi_{2nu}) is not identically zero.  Then its zeros on the arc are
isolated, so a two-side-length star C_N polygon is rigid (finitely many shapes per word)."""
import sys
N=int(sys.argv[1])
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def idiv(a,b):
    a=a[:];q=[0]*max(1,len(a)-len(b)+1)
    for k in range(len(a)-len(b),-1,-1):
        q[k]=a[k+len(b)-1]
        for j,x in enumerate(b):a[k+j]-=q[k]*x
    assert all(x==0 for x in a); return trim(q)
cyclo={}
for n in range(1,2*N+2):
    p=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:p=idiv(p,cyclo[d])
    cyclo[n]=p
def reduce(a,mod):
    d=len(mod)-1; a=list(a)+[0]*max(0,d-len(a))
    for k in range(len(a)-1,d-1,-1):
        c=a[k]
        for j in range(d+1):a[k-d+j]-=c*mod[j]
    return tuple(a[:d])
bad=[]
for x in range(1,2**N-1):
    w=format(x,'0%db'%N); v=w.count('1'); mod=cyclo[2*v]; d=len(mod)-1
    B=[0]
    for b in w[:-1]:B.append(B[-1]+int(b))
    P={}
    for k in range(N):
        e=v*k-N*B[k]; c=[0]*(2*v); c[B[k]%(2*v)]=(-1)**k
        P[e]=tuple(p+q for p,q in zip(P.get(e,(0,)*d),reduce(c,mod)))
    if all(all(c==0 for c in coef) for coef in P.values()): bad.append(w)
print(f"N={N}: {2**N-2} nonconstant binary words; closure polynomial identically zero for {len(bad)} of them", bad[:10])
