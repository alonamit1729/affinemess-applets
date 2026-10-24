"""Real-isolated points on the 1-dimensional components of (5,4) types.
A real point of a complex curve that is isolated in the real locus must be a singular point of the curve, so:
singular locus (Singular, sing.lib) -> real singular points (msolve, with the inequations) -> sphere test:
if the curve has no real point on a small sphere around the point, the point is real-isolated (rigorous for that radius)."""
import json, subprocess, re, itertools, math
from fractions import Fraction
exec(open('certify54.py').read().split('results=[]')[0])   # reuse helpers: run_msolve, parse_real, expand, ineq, VARS...
C=json.load(open('cert54.json'))
V=[v for v in VARS if v!='t']
def singular_locus(gens):
    ring='ring r=0,('+','.join(V)+'),dp;'
    script=ring+'\nLIB "sing.lib";\nideal P='+',\n'.join(gens)+';\nideal S=std(slocus(P));\nprint(dim(S));\nprint(S);\nquit;\n'
    open('sl_in.sing','w').write(script)
    r=subprocess.run([SG,'-q','sl_in.sing'],capture_output=True,text=True,timeout=1800)
    lines=[l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]
    dim=int(lines[0]); gens=[l.rstrip(',') for l in lines[1:]]
    return dim,gens
found=[]
for r in C:
    for pi,p in enumerate(r['primes']):
        if p['dim']!=1 or p['degenerate']: continue
        sdim,sg=singular_locus(p['gens'])
        if sdim<0: print('type',r['idx'],'prime',pi,'smooth'); continue
        if sdim>0: print('type',r['idx'],'prime',pi,'singular locus of dim',sdim,'!?'); continue
        o=run_msolve(sg+[expand(ineq)])
        st,md=parse_real(o)
        print('type',r['idx'],'prime',pi,'singular points: status',st,'real',len(md))
        for mid in md:
            # sphere test at two radii around the (rational, 200-bit) midpoint
            isolated=True
            for eps in (Fraction(1,100),Fraction(1,1000)):
                sph='+'.join(f'({v}-({mid[k]}))^2' for k,v in enumerate(V))+f'-({eps*eps})'
                o2=run_msolve(p['gens']+[expand(ineq),expand(sph)])
                st2,md2=parse_real(o2)
                if st2!='zerodim': print('   sphere test status',st2); isolated=False; break
                if md2: isolated=False; break
            P=points_from(mid)
            print('   real singular point',[(round(x,4),round(y,4)) for x,y in P],'isolated:',isolated)
            if isolated: found.append({'idx':r['idx'],'prime':pi,'pts':P})
json.dump(found,open('isolated54.json','w'))
print('real-isolated points found:',len(found))
