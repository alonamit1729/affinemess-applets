#!/usr/bin/env python3
"""Star chords of an equilateral polygon.  Edge k has direction k*pi - sum_{j<k} theta_j.
chord_r(theta_1..theta_{r-1}) = |e_0+...+e_{r-1}|^2 (r unit edges, r-1 interior angles).
For the step-5 star on 11 points: r=2 -> convex side, r=3 -> convex step 4,
r=4 -> convex step 2, r=5 -> convex step 3.  Check signs of partials on (0,pi/3)^(r-1)."""
import sympy as sp, itertools, numpy as np
def chord_sq(th):
    r=len(th)+1
    phi=[k*sp.pi - sum(th[:k]) for k in range(r)]
    return sp.simplify(sum(sp.cos(phi[j]-phi[k]) for j in range(r) for k in range(r)))
for r in (3,4,5):
    th=sp.symbols('t1:%d'%r)
    F=chord_sq(th)
    print("r=%d:"%r, sp.simplify(F))
    grads=[sp.lambdify(th,sp.diff(F,t),'numpy') for t in th]
    pts=np.linspace(0.001,np.pi/3-0.001,25)
    for i,gfun in enumerate(grads):
        vals=np.array([gfun(*p) for p in itertools.product(pts,repeat=r-1)])
        print("   d/dtheta_%d: min %.4f max %.4f -> sign %s"%(i+1,vals.min(),vals.max(),'+' if vals.min()>0 else ('-' if vals.max()<0 else 'mixed')))
