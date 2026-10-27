#!/usr/bin/env python3
"""Exact integer-arithmetic certificate, generalised from Amit-Astra certify_cycle13.py
(Appendix A) to an odd diameter cycle C_N (N odd), i.e. the step-(N-1)/2 star on N
convexly ordered diameter endpoints.

Claim certified: if the boundary sides of the convex N-gon take exactly TWO lengths
a<b (angles alpha<beta at the star vertices, u*alpha+nu*beta=pi, 1<=nu<=N-1), the
closure polynomial of the equilateral star has no zero on the open arc
0<theta<pi/(N nu), theta=alpha/nu.  Hence the star cannot close: the two-length
case is impossible, and with one length the star is regular.

Derivation (identical to Section 6 of the paper with 13 -> N): edge k has direction
k*pi - sum_{j<k} theta_j = k*pi - B_k*pi/nu - theta*(nu*k - N*B_k), so
P_b(z) = sum_k (-1)^k omega^{B_k} z^{nu k - N B_k}, omega=e^{-i pi/nu}, z=e^{-i theta}.
F_nu(z) = (z^N - omega)/(z - eta), eta = omega^r, r = N^{-1} mod 2nu, divides P~_b;
its only zero on the closed arc is the regular endpoint theta = pi/(N nu).
Usage: python3 certify_cycle_odd.py N [--all-words]
"""
from fractions import Fraction
from math import factorial
from pathlib import Path
from functools import lru_cache
import json,time,sys
if not __debug__:
 raise RuntimeError("Run without -O: proof checks require assertions.")
N=int(sys.argv[1]); assert N%2==1 and N>=5
S=10**50
assert 4**121*S < factorial(121)

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def idiv(a,b):
 a=a[:];q=[0]*max(1,len(a)-len(b)+1)
 for k in range(len(a)-len(b),-1,-1):
  assert b[-1]==1
  q[k]=a[k+len(b)-1]
  for j,x in enumerate(b):a[k+j]-=q[k]*x
 assert all(x==0 for x in a)
 return trim(q)
cyclo={}
for n in range(1,2*N+2):
 p=[-1]+[0]*(n-1)+[1]
 for d in range(1,n):
  if n%d==0:p=idiv(p,cyclo[d])
 cyclo[n]=p

def atan_bounds(inv, terms=100):
 p=Fraction(0)
 for k in range(terms):p+=Fraction((-1)**k,(2*k+1)*inv**(2*k+1))
 e=Fraction(1,(2*terms+1)*inv**(2*terms+1))
 return (p,p+e) if terms%2==0 else (p-e,p)
a,b=atan_bounds(5);c,d=atan_bounds(239)
pi_lo,pi_hi=16*a-4*d,16*b-4*c
PI=(pi_lo.numerator*S//pi_lo.denominator,-((-pi_hi.numerator*S)//pi_hi.denominator))
assert 3*S<PI[0]<=PI[1]<4*S

def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 p=[x*y for x in a for y in b]
 return(min(p)//S,-((-max(p))//S))
def divi(a,n):return(a[0]//n,-((-a[1])//n))
def scale(a,n):return(a[0]*n,a[1]*n) if n>=0 else(a[1]*n,a[0]*n)

@lru_cache(None)
def trig_pi(num,den):
 x=divi(scale(PI,num),den);assert 0<=x[0]<=x[1]<=4*S
 xx=mul(x,x);st=x;ct=(S,S);sn=x;cs=ct
 for k in range(1,61):
  st=divi(mul(st,xx),(2*k)*(2*k+1))
  ct=divi(mul(ct,xx),(2*k-1)*(2*k))
  sn=add(sn,st if k%2==0 else neg(st))
  cs=add(cs,ct if k%2==0 else neg(ct))
 return((sn[0]-1,sn[1]+1),(cs[0]-1,cs[1]+1))

def cmul(a,b):
 return(sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0])))
def cadd(a,b):return(add(a[0],b[0]),add(a[1],b[1]))
def zero_excluded(a):return a[0][0]>0 or a[0][1]<0 or a[1][0]>0 or a[1][1]<0

class Ring:
 def __init__(self,v):
  self.v=v;self.mod=cyclo[2*v];self.d=len(self.mod)-1;self.zero=(0,)*self.d;self.one=(1,)+(0,)*(self.d-1)
  self.omega=self.reduce([0,1]);self.powers=[self.one]
  for k in range(1,2*v):self.powers.append(self.mul(self.powers[-1],self.omega))
 def reduce(self,a):
  a=list(a)+[0]*max(0,self.d-len(a))
  for k in range(len(a)-1,self.d-1,-1):
   c=a[k]
   for j in range(self.d+1):a[k-self.d+j]-=c*self.mod[j]
  return tuple(a[:self.d])
 def add(self,a,b):return tuple(x+y for x,y in zip(a,b))
 def neg(self,a):return tuple(-x for x in a)
 def mul(self,a,b):
  out=[0]*(2*self.d-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return self.reduce(out)
 def interval(self,a):
  re=(0,0);im=(0,0)
  for j,x in enumerate(a):
   if x:
    sn,cs=trig_pi(j,self.v)
    re=add(re,scale(cs,x));im=sub(im,scale(sn,x))
  return(re,im)

def canon(s):return min(s[i:]+s[:i] for s in(s,s[::-1]) for i in range(N))
FULL = '--all-words' in sys.argv
fmt='0%db'%N
words=([format(x,fmt) for x in range(1,2**N-1)] if FULL else
       sorted({canon(format(x,fmt)) for x in range(1,2**N-1)}))
allcert=[];start=time.time()
for iw,w in enumerate(words):
 v=w.count('1');ring=Ring(v);B=[0]
 for b in w[:-1]:B.append(B[-1]+int(b))
 exps=[v*k-N*B[k] for k in range(N)];lo=min(exps);hi=max(exps)
 P=[ring.zero]*(hi-lo+1)
 for k,e in enumerate(exps):
  c=ring.powers[B[k]%(2*v)]
  if k%2:c=ring.neg(c)
  P[e-lo]=ring.add(P[e-lo],c)
 r=pow(N,-1,2*v)
 F=[ring.powers[(r*(N-1-j))%(2*v)] for j in range(N)]   # F_nu = sum_j eta^{N-1-j} z^j, eta=omega^r
 rem=P[:];Q=[ring.zero]*(len(P)-(N-1))
 for k in range(len(P)-N,-1,-1):
  c=rem[k+N-1];Q[k]=c
  for j,f in enumerate(F):rem[k+j]=ring.add(rem[k+j],ring.neg(ring.mul(c,f)))
 assert all(a==ring.zero for a in rem),(w,'regular factor failure')
 endpoint_order=0
 while True:
  val=ring.zero
  for a in Q:val=ring.add(val,a)
  if val!=ring.zero:break
  assert len(Q)>1
  q=[ring.zero]*(len(Q)-1);q[-1]=Q[-1]
  for k in range(len(q)-2,-1,-1):q[k]=ring.add(Q[k+1],q[k+1])
  assert ring.add(Q[0],q[0])==ring.zero
  Q=q;endpoint_order+=1
 coeff=[ring.interval(a) for a in Q]
 stack=[(0,1,0)];accepted=[];maxdepth=0
 while stack:
  left,right,depth=stack.pop();den=2**depth
  sl,cl=trig_pi(left,N*v*den);sr,cr=trig_pi(right,N*v*den)
  z=((cr[0],cl[1]),(-sr[1],-sl[0]))
  val=coeff[-1]
  for c in reversed(coeff[:-1]):val=cadd(cmul(val,z),c)
  if zero_excluded(val):
   accepted.append((left,right,depth));maxdepth=max(maxdepth,depth)
  else:
   if depth>=30:raise RuntimeError((w,'inconclusive',left,right,depth,val))
   stack.append((left*2+1,right*2,depth+1));stack.append((left*2,left*2+1,depth+1))
 allcert.append({'word':w,'v':v,'degree_P':len(P)-1,'degree_Q':len(Q)-1,'removed_z1':endpoint_order,'intervals':len(accepted),'max_depth':maxdepth})
out=Path(__file__).with_name('cycle%d_%s_certificate.json'%(N,'allwords' if FULL else 'dihedral'))
out.write_text(json.dumps({'N':N,'scale':str(S),'words':len(words),'total_intervals':sum(c['intervals'] for c in allcert),'max_depth':max(c['max_depth'] for c in allcert),'certificates':allcert},indent=1))
print('PROVED N=%d:'%N,len(words),'cases;',sum(c['intervals'] for c in allcert),'intervals; max depth',max(c['max_depth'] for c in allcert),'elapsed',round(time.time()-start,1),flush=True)
