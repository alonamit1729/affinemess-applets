#!/usr/bin/env python3
"""Heuristic test: for 12 points in convex position whose diameter graph is a perfect
matching of 6 pairwise-crossing unit segments p_i p_{i+6}, is |p_0p_1|<|p_0p_2|<...<|p_0p_5|
forced?  Random perturbations of the regular 12-gon's six diameters (direction, lateral
offset, slide along the line), keeping only samples whose diameter graph is the matching."""
import numpy as np
rng=np.random.default_rng(0)
def convex_cyclic(P):
    n=len(P)
    for i in range(n):
        a,b,d=P[i],P[(i+1)%n],P[(i+2)%n]
        if (b[0]-a[0])*(d[1]-a[1])-(b[1]-a[1])*(d[0]-a[0])<=0: return False
    return True
viol=0; tot=0; worst=None; fails={'dist':0,'convex':0}
for scale in (0.02,0.05,0.1,0.2):
  for trial in range(100000):
    th=np.array([i*np.pi/6 for i in range(6)])+rng.uniform(-scale,scale,6)
    off=rng.uniform(-scale,scale,6); sh=rng.uniform(-scale,scale,6)
    P=np.zeros((12,2))
    for i in range(6):
        d=np.array([np.cos(th[i]),np.sin(th[i])]); nrm=np.array([-d[1],d[0]])
        P[i+6]=off[i]*nrm+(0.5+sh[i])*d; P[i]=off[i]*nrm-(0.5-sh[i])*d
    D=np.linalg.norm(P[:,None]-P[None],axis=2)
    mask=np.ones((12,12),bool); np.fill_diagonal(mask,False)
    for i in range(6): mask[i,i+6]=mask[i+6,i]=False
    if D[mask].max()>=1: fails['dist']+=1; continue
    if not convex_cyclic(P): fails['convex']+=1; continue
    tot+=1
    for i in range(12):
        seq=[D[i,(i+t)%12] for t in range(1,6)]; seq2=[D[i,(i-t)%12] for t in range(1,6)]
        if any(seq[j]>=seq[j+1] for j in range(4)) or any(seq2[j]>=seq2[j+1] for j in range(4)):
            viol+=1; worst=(scale,i,np.round(seq,3),np.round(seq2,3)); break
print("valid samples:",tot,"rejected:",fails,"samples with a chain violation:",viol,"example:",worst)
