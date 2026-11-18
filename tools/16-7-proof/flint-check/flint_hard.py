"""Hard patterns: try partial saturation polynomials (any partial product of the 'non-zero / distinct' factors is enough
for a completeness proof: if V(planarity) and V(B) agree where q_part != 0, they agree where the full q != 0)."""
import sys, time, faulthandler
src=open('flint_all_origins.py').read().split("def one(i):")[0]
exec(src)
i=int(sys.argv[2]); mode=sys.argv[3]; faulthandler.dump_traceback_later(int(sys.argv[4]),exit=True)
xs,w,gens,_=setup(P[i]); C=[sp.Integer(1),*xs]
B=sp.sympify(saved[i],locals={str(x):x for x in xs})
t=time.time()
if mode=='none':
    G=gb(gens,xs); G2=gb(B,xs); syms=xs
else:
    if mode=='zeros': fac=list(xs)
    elif mode=='zeros1': fac=[*xs,*[x-1 for x in xs]]
    elif mode.startswith('x'): fac=[sp.Symbol(mode)]
    qp=sp.Integer(1)
    for f in fac: qp*=f
    sat=sp.expand(w*qp-1); syms=(w,*xs)
    G=gb([*gens,sat],syms); G2=gb([*B,sat],syms)
print('pattern',i,'mode',mode,'equality',sorted(map(str,G))==sorted(map(str,G2)),'solutions',nsol(G,len(syms)),'vs',nsol(G2,len(syms)),'%.1fs'%(time.time()-t),flush=True)
