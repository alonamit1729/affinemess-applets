import math, itertools, json
import numpy as np
from scipy.optimize import least_squares, brentq, minimize_scalar
S3=math.sqrt(3)
# ---- 4-point families (from fam.py) ----
def A(t):  return [(0,0),(1,0),(.5,S3/2),(math.cos(t),math.sin(t))]
def B(al): return [(0,0),(1,0),(math.cos(al),math.sin(al)),(1+math.cos(al),math.sin(al))]
def C(ph): return [(0,0),(1,0),(math.cos(ph),math.sin(ph)),(math.cos(ph),-math.sin(ph))]
def D(t):  return [(0,0),(1,0),(.5,S3/2),(.5,t)]
def Et(th):return [(0,0),(1,0),(math.cos(th),math.sin(th)),(1-math.cos(th),math.sin(th))]
def Ep(b):
    x2=(2-b*b)/2; y2=math.sqrt(max(1-x2*x2,0)); return [(0,0),(1,0),(x2,y2),(1-x2,-y2)]
def F(h,s=1):
    P1=(0,0);P2=(1,0);P0=(.5,h); b=math.hypot(.5,h); mid=(.25,h/2); d=(-h,.5); n=math.hypot(*d)
    t=math.sqrt(max(1-b*b/4,0)); P3=(mid[0]+s*t*d[0]/n, mid[1]+s*t*d[1]/n); return [P0,P1,P2,P3]
def G(s):  return [(0,s),(-.5,0),(.5,0),(0,s+1)]
def H(ph,s=1):
    P2=(0,0);P1=(1,0);P3=(math.cos(2*ph),math.sin(2*ph)); b=2*math.sin(ph); dx=math.cos(2*ph)-.5; r=b*b-dx*dx
    if r<0: return None
    return [(.5,math.sin(2*ph)+s*math.sqrt(r)),P1,P2,P3]
def I(r):  return [(0,0),(1,0),(0,r),(1,r)]
phmin=0.2087
families=[
 ('1',A,1e-3,2*math.pi-1e-3),('2',B,1e-3,math.pi-1e-3),('3',C,1e-3,math.pi-1e-3),('4',D,-8,8),
 ('5i',Et,1e-3,math.pi-1e-3),('5ii',Ep,1e-3,2-1e-3),('6+',lambda h:F(h,1),0,math.sqrt(15)/2),('6-',lambda h:F(h,-1),0,math.sqrt(15)/2),
 ('7',G,-0.499,8),('8+',lambda p:H(p,1),phmin+1e-4,math.pi/2),('8-',lambda p:H(p,-1),phmin+1e-4,math.pi/2),('9',I,1e-3,1-1e-3)]
# ---- 2-distance sets ----
k1=math.sqrt(2-S3);
two=[
 ('square',[(0,0),(1,0),(1,1),(0,1)]),
 ('rhombus60',[(0,0),(1,0),(.5,S3/2),(1.5,S3/2)]),
 ('tri+centre',[(0,0),(1,0),(.5,S3/2),(.5,S3/6)]),
 ('pent-trap',[(math.cos(2*math.pi*k/5),math.sin(2*math.pi*k/5)) for k in range(4)]),
 ('kite',[(0,0),(1,0),(.5,S3/2),(S3/2,.5)]),
 ('dart',[(0,0),(1,0),(.5,S3/2),(-S3/2,-.5)]),
]
def pdist(P): return [math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)]
def distinct(vals,tol):
    v=sorted(vals); out=[v[0]]
    for x in v[1:]:
        if x-out[-1]>tol: out.append(x)
    return out
def circ2(c1,r1,c2,r2):
    d=math.dist(c1,c2)
    if d<1e-12: return []
    a=(r1*r1-r2*r2+d*d)/(2*d); h2=r1*r1-a*a
    if h2<-1e-9: return []
    h=math.sqrt(max(h2,0)); ex=((c2[0]-c1[0])/d,(c2[1]-c1[1])/d); m=(c1[0]+a*ex[0],c1[1]+a*ex[1])
    return [(m[0]-h*ex[1],m[1]+h*ex[0]),(m[0]+h*ex[1],m[1]-h*ex[0])]
def circumcenter(p,q,r):
    ax,ay=p;bx,by=q;cx,cy=r
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-10: return None
    ux=((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d
    uy=((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d
    return (ux,uy)
# ---- polish a 5-point candidate on its colouring, return polished points or None ----
def polish(P):
    ds=pdist(P); vals=distinct(ds,0.02)
    if len(vals)!=3: return None
    E=list(itertools.combinations(range(5),2))
    cls=[min(range(3),key=lambda k:abs(vals[k]-d)) for d in ds]
    if len(set(cls))!=3: return None
    def res(z):
        Q=[(0,0),(z[0],0)]+[(z[1+2*i],z[2+2*i]) for i in range(3)]
        v=z[7:10]
        return [math.dist(Q[i],Q[j])-v[cls[k]] for k,(i,j) in enumerate(E)]
    # initial: align P so that P0 at origin, P1 on x-axis
    p0=P[0]; p1=P[1]; ang=math.atan2(p1[1]-p0[1],p1[0]-p0[0]); c,s=math.cos(-ang),math.sin(-ang)
    Q=[((x-p0[0])*c-(y-p0[1])*s,(x-p0[0])*s+(y-p0[1])*c) for x,y in P]
    z0=[Q[1][0]]+[w for q in Q[2:] for w in q]+vals
    r=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15)
    if max(abs(x) for x in r.fun)>1e-9: return None
    z=r.x; Q=[(0,0),(z[0],0)]+[(z[1+2*i],z[2+2*i]) for i in range(3)]
    ds=pdist(Q)
    if min(ds)<1e-6: return None
    if len(distinct(ds,1e-6))!=3: return None
    return Q
def dm(P): return [[math.dist(P[i],P[j]) for j in range(len(P))] for i in range(len(P))]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(5)):
        s=D2[p[0]][p[1]]/D1[0][1]
        if all(abs(D1[i][j]*s-D2[p[i]][p[j]])<1e-6 for i in range(5) for j in range(i+1,5)): return True
    return False
found=[]
def add(P,src):
    Q=polish(P)
    if Q is None: return
    for q,_ in found:
        if congruent(q,Q): return
    found.append((Q,src)); print('NEW',len(found),src,[round(v,5) for v in distinct(pdist(Q),1e-6)])
# ---- Case A: extend 2-distance sets ----
for name,Qs in two:
    vals=distinct(pdist(Qs),1e-9); a,b=vals
    for pat in itertools.product([0,1,2],repeat=4):
        nc=pat.count(2)
        if nc==0: continue
        known=[i for i in range(4) if pat[i]!=2]
        cands=[]
        if nc<=2:
            i,j=known[0],known[1]
            cands=circ2(Qs[i],vals[pat[i]],Qs[j],vals[pat[j]])
        elif nc==3:
            cp=[i for i in range(4) if pat[i]==2]
            cc=circumcenter(Qs[cp[0]],Qs[cp[1]],Qs[cp[2]])
            if cc: cands=[cc]
        else:
            cc=circumcenter(Qs[0],Qs[1],Qs[2])
            if cc: cands=[cc]
        for P in cands:
            ok=True; cval=None
            for i in range(4):
                d=math.dist(P,Qs[i])
                if pat[i]==2:
                    if cval is None: cval=d
                    elif abs(d-cval)>1e-7: ok=False
                else:
                    if abs(d-vals[pat[i]])>1e-7: ok=False
            if ok and cval is not None and min(abs(cval-a),abs(cval-b))>1e-6 and min(math.dist(P,q) for q in Qs)>1e-6:
                add(Qs+[P],name)
# ---- Case B: extend 3-distance families ----
N=3000
for fname,f,lo,hi in families:
    ts=[lo+(hi-lo)*i/N for i in range(N+1)]
    for pat in itertools.product([0,1,2],repeat=4):
        for br in (0,1):
            def g(t):
                Q=f(t)
                if Q is None: return float('nan')
                vals=distinct(pdist(Q),1e-9)
                if len(vals)!=3: return float('nan')
                c=circ2(Q[0],vals[pat[0]],Q[1],vals[pat[1]])
                if len(c)<2: return float('nan')
                P=c[br]
                return math.dist(P,Q[2])-vals[pat[2]]
            gs=[g(t) for t in ts]
            roots=[]
            for i in range(1,N+1):
                a_,b_=gs[i-1],gs[i]
                if math.isnan(a_) or math.isnan(b_): continue
                if a_*b_<0:
                    try: roots.append(brentq(g,ts[i-1],ts[i]))
                    except Exception: pass
                elif i<N and not math.isnan(gs[i+1]) and abs(b_)<abs(a_) and abs(b_)<=abs(gs[i+1]) and abs(b_)<2e-3:
                    r=minimize_scalar(lambda t:abs(g(t)) if not math.isnan(g(t)) else 1.0,bounds=(ts[i-1],ts[i+1]),method='bounded',options={'xatol':1e-12})
                    if r.fun<1e-7: roots.append(r.x)
            for t in roots:
                Q=f(t)
                if Q is None: continue
                vals=distinct(pdist(Q),1e-9)
                if len(vals)!=3: continue
                c=circ2(Q[0],vals[pat[0]],Q[1],vals[pat[1]])
                if len(c)<2: continue
                P=c[br]
                if abs(math.dist(P,Q[3])-vals[pat[3]])<1e-5 and min(math.dist(P,q) for q in Q)>1e-5:
                    add(Q+[P],fname+str(pat))
print('TOTAL',len(found))
json.dump([{'pts':Q,'src':s} for Q,s in found],open('five.json','w'))
