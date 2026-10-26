"""Trace the one-dimensional components of the (5,5) types that also have two-dimensional components (226, 250, 270):
decompose the type's equations (Singular minAssGTZ), run the tracer on each non-degenerate one-dimensional prime with x1 = 1."""
import sys, json, os
sys.argv=['strat.py','5']
c5=open('certify5d.py').read().split('import os')[0]
exec(c5)   # classes, system, decompose (over V with s1..s4), ...
_classes,_system,_decompose=classes,system,decompose
exec(open('strat.py').read().split("if os.path.exists(f'sat5{d}.json'):")[0])
classes,system,decompose=_classes,_system,_decompose
for idx in (226,250,270):
    if any(o['type']==idx for o in out): continue
    col=classes[idx]; primes=decompose(system(col)+['x1-1'])
    print(idx,'primes',[(p['dim'],p['degenerate']) for p in primes],flush=True)
    for pi,p in enumerate(primes):
        if p['dim']!=1 or p['degenerate']: continue
        process(idx,list(col),pi,p['gens'])
print('done')
