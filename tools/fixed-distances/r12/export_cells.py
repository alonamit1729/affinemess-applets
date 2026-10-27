"""Export selected cells of the R12+R12 sweep for the applet builder: r12sum.json = {"n,d": [{pts, exact, shells}]}.
Also checks the sweep against the exactly classified cells (every subset must be one of the catalog's configurations)."""
import json, math, cmath, itertools, sys, collections
R=json.load(open('r12_subsets_7.json'))
SEL={(8,5),(9,5),(10,5),(11,5),(12,5),(9,6),(10,6),(11,6),(12,6),(13,6),(13,7),(14,7),(15,7),(16,7)}
# exact coordinates: every point is z^a+z^b, z=e^{i pi/6}; coordinates are (p+q*sqrt3)/2 with small integers p,q
def exact_coord(x):
    def term(c,rad):   # c*rad with c integer
        if c==0: return ''
        sgn='−' if c<0 else '+'; c=abs(c)
        return sgn+(('' if c==1 and rad else str(c))+rad)
    def expr(p,q):
        s=(term(p,'')+term(q,'√3')).lstrip('+')
        return s if s else '0'
    for q in range(-4,5):
        p=2*x-q*math.sqrt(3)
        if abs(p-round(p))<1e-9:
            p=int(round(p))
            if p%2==0 and q%2==0: return expr(p//2,q//2)
            e=expr(p,q)
            return f'({e})/2' if ('+' in e[1:] or '−' in e[1:]) else e+'/2'
    raise ValueError(x)
out={}
for key,L in R.items():
    n,d=map(int,key.split(','))
    if (n,d) not in SEL: continue
    recs=[]
    for rec in L:
        P=[tuple(p) for p in rec['pts']]
        ex=[[exact_coord(x),exact_coord(y)] for x,y in P]
        shells=collections.Counter()
        for x,y in P:
            r=math.hypot(x,y)
            k=min(range(7),key=lambda k: abs(2*math.cos(math.pi*k/12)-r)); shells[k]+=1
        recs.append({'pts':P,'exact':ex,'shells':sorted(shells.items())})
    out[key]=recs
json.dump(out,open('r12sum.json','w'))
print({k:len(v) for k,v in out.items()})
