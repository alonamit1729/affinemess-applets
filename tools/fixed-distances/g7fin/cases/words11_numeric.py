#!/usr/bin/env python3
"""Floating-point reconnaissance for the two surviving words: minimise the sum of squared
residuals of the forced coincidences over the admissible angle domain (heuristic)."""
import numpy as np
from scipy.optimize import minimize
from math import pi
def chord(*th):
    r=len(th)+1; phi=[k*pi-sum(th[:k]) for k in range(r)]
    return sum(np.cos(phi[j]-phi[k]) for j in range(r) for k in range(r))
def side(x): return chord(x)
SYS={
 'aaaaaaaabcb':(lambda a,b,c:[side(b)-chord(a,a,a), side(c)-chord(a,a,b), chord(a,a,b)-chord(a,a,b,c),
   chord(a,a,b,c)-chord(a,b,c,b), chord(b,c)-chord(a,b,c), chord(a,b,c)-chord(a,a,a,b), chord(a,b)-chord(a,a,a,a)],
   (8,2,1)),
 'ababababccb':(lambda a,b,c:[side(b)-chord(a,b,a), side(c)-chord(b,a,b), chord(b,a,b)-chord(a,b,c,c),
   chord(c,c)-chord(a,b,c), chord(a,b,c)-chord(b,c,c), chord(b,c,c)-chord(b,a,b,c), chord(b,a,b,c)-chord(b,c,c,b),
   chord(b,c)-chord(a,b,a,b)],(4,5,2)),
}
rng=np.random.default_rng(1)
for w,(res,(u,v,s)) in SYS.items():
    best=(1e9,None)
    for _ in range(3000):
        a,b=rng.uniform(0,pi/3,2)
        if a>b: a,b=b,a
        f=lambda x:(lambda a,b,c:sum(r*r for r in res(a,b,c)))(x[0],x[1],(pi-u*x[0]-v*x[1])/s)
        cons=[{'type':'ineq','fun':lambda x:x[1]-x[0]},{'type':'ineq','fun':lambda x:(pi-u*x[0]-v*x[1])/s-x[1]},
              {'type':'ineq','fun':lambda x:pi/3-(pi-u*x[0]-v*x[1])/s},{'type':'ineq','fun':lambda x:x[0]}]
        r=minimize(f,[a,b],constraints=cons,method='SLSQP',options={'ftol':1e-14,'maxiter':500})
        if r.fun<best[0]: best=(r.fun,r.x)
    a,b=best[1]; c=(pi-u*a-v*b)/s
    print(w,"min sum of squared residuals over admissible (a<b<c<pi/3):",best[0],"at a,b,c=",a,b,c)
    # also the single pair of 'cheapest' equations
