"""10-point sets with <=6 distances from the two-concentric-pentagon families: scan the parameter for coincidences of
distance values (endpoints of the (10,7) families), refine by bisection, and output the sets."""
import math, itertools, json, numpy as np
def pent(r,phi): return [(r*math.cos(phi+2*math.pi*j/5), r*math.sin(phi+2*math.pi*j/5)) for j in range(5)]
def vals(P,tol=1e-9):
    ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)); v=[ds[0]]
    for x in ds[1:]:
        if x-v[-1]>tol*ds[-1]: v.append(x)
    return v
def mingap(P):
    v=vals(P,1e-12); return min(b-a for a,b in zip(v,v[1:]))/v[-1]
fams={'aligned':lambda t: pent(1,0)+pent(t,0),'offset':lambda t: pent(1,0)+pent(t,math.pi/5),'samer':lambda t: pent(1,t)+pent(1,-t)}
out=[]
for name,f in fams.items():
    ts=np.linspace(0.02,6,60000) if name!='samer' else np.linspace(0.001,math.pi/5-0.001,60000)
    g=[mingap(f(t)) for t in ts]
    for i in range(1,len(ts)-1):
        if g[i]<g[i-1] and g[i]<g[i+1] and g[i]<1e-3:
            # refine: golden-section on mingap
            a,b=ts[i-1],ts[i+1]
            for _ in range(80):
                m1=a+(b-a)/3; m2=b-(b-a)/3
                if mingap(f(m1))<mingap(f(m2)): b=m2
                else: a=m1
            t=(a+b)/2; P=f(t); k=len(vals(P,1e-7))
            if k<=6 and min(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(10),2))>1e-6:
                out.append({'family':name,'t':t,'k':k,'pts':P})
print(len(out),'sets;',sorted((o['family'],round(o['t'],5),o['k']) for o in out))
json.dump(out,open('g7fin/pent_ends.json','w'))
