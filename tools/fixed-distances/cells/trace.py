import sys,os,json,random
n=int(sys.argv[1]); d=int(sys.argv[2]); idx=sys.argv[3]
sys.argv=['cellcert.py',str(n),str(d)]
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,'cellcert.py')).read(); exec(src[:src.index('results = json.load')])
R=json.load(open(os.path.join(HERE,f'cert_{n}_{d}.json'))); col=tuple(R[idx]['col'])
V,polys=system_gram(n,col); SV=V[:-1]; rab=polys[-1]
facs=[f's{a}-s{b}' for a,b in itertools.combinations(range(1,d),2)]+[f's{c}-1' for c in range(1,d)]+[f's{c}' for c in range(1,d)]
primes=decompose(SV,polys[:-1],facs)
p=[p for p in primes if not p['degenerate'] and p['dim']>0][0]
print(p['gens']); print('slocus',p['slocus'],p['sdim'])
random.seed(3)
lins=[linform(SV)]
print(lins)
out=msolve_real_raw(V,p['gens']+lins+[rab],precision=200)
print(out[:500])
s2=parse_real(out)
import sympy as sp
syms={v:sp.Symbol(v) for v in V}
for sol in (s2 if isinstance(s2,list) else []):
    vals=sol_vals(V,sol); print({k:float(v) for k,v in vals.items()})
    for g in p['gens']+lins: print('  ',float(sp.sympify(g.replace('^','**'),locals=syms).subs({syms[k]:v for k,v in vals.items()})))
