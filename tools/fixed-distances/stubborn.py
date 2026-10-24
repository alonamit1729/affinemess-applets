"""Fallback for types whose primary decomposition did not finish: dimension of the saturated ideal (Rabinowitsch variable t
in the ring, so the inequations are enforced), certified real points by random slices, and the Jacobian rank at the
real points (exact rational arithmetic at the isolation-box midpoint is not on the variety, so rank is reported numerically).
Usage: python3 stubborn.py d idx [timeout]"""
import sys, json, subprocess, time, itertools, math
import numpy as np
d=int(sys.argv[1]); idx=int(sys.argv[2]); tmo=int(sys.argv[3]) if len(sys.argv)>3 else 3600
src=open('certify5d.py').read().split('import os')[0]
sys.argv=[sys.argv[0],str(d)]; exec(src)
col=classes[idx]; eqs=system(col)+[expand(ineq)]; VT=V+['t']
ring='ring r=0,('+','.join(VT)+'),dp;'
script=ring+'\nideal i='+',\n'.join(eqs)+';\nideal G=std(i);\nprint("DIM"); print(dim(G));\nquit;\n'
open(f'sb_{d}_{idx}.sing','w').write(script)
t0=time.time()
try:
    r=subprocess.run([SG,'-q',f'sb_{d}_{idx}.sing'],capture_output=True,text=True,timeout=tmo)
    out=[l.strip() for l in r.stdout.splitlines() if l.strip()]
    dimv=int(out[out.index('DIM')+1])
except subprocess.TimeoutExpired:
    dimv=None
print('type',idx,'col',col,'dim of saturated ideal (t included):',dimv,'%.0fs'%(time.time()-t0),flush=True)
e=2*5-4-10+d
# slices: e random hyperplanes in the point variables, solve with the inequation included
random.seed(idx); hits=0; pts=[]
for k in range(6):
    lins=['+'.join(f'({random.randint(-5,5)})*{v}' for v in XV)+f'-({random.randint(-3,3)})' for _ in range(e)]
    st,md=parse_real(run_msolve(eqs+[expand(l) for l in lins],VT))
    if st=='zerodim' and md: hits+=1; pts+=[m for m in md]
    print('  slice',k,st,len(md),flush=True)
# Jacobian rank at real points (numerical, full-precision midpoints)
syms=[SYM[v] for v in V]
F=[sp.sympify(q.replace('^','**'),locals=SYM) for q in system(col)]
J=sp.Matrix(F).jacobian(syms)
ranks=[]
for m in pts[:8]:
    sub={SYM[v]:sp.Rational(m[i]) for i,v in enumerate(V)}
    Jn=np.array(J.subs(sub).evalf(30).tolist(),dtype=float); sv=np.linalg.svd(Jn,compute_uv=False)
    ranks.append((int((sv>1e-9*sv[0]).sum()), float(sv[-1]/sv[0]) if sv[-1]>0 else 0.0))
print('  real slice hits',hits,'of 6; Jacobian ranks at real points (rank, min/max singular value):',ranks)
print('  expected dim',e,'; codim for a smooth point on an e-dim component =',len(V)-e)
json.dump({'idx':idx,'col':list(col),'dim_saturated':dimv,'slice_hits':hits,'ranks':ranks,'points':[points_from(m) for m in pts[:4]]},open(f'stubborn_{d}_{idx}.json','w'))
