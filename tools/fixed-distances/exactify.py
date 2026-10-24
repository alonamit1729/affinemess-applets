"""Put a configuration into a natural exact frame (regular polygon / triangular lattice / square lattice /
collinear) and produce exact coordinate strings where the structure allows; decimals otherwise."""
import math, itertools
from fractions import Fraction
S3=math.sqrt(3)
def circumcenter(p,q,r):
    ax,ay=p;bx,by=q;cx,cy=r
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-10: return None
    ux=((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d
    uy=((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d
    return (ux,uy)
def transform(P,c,scale,ang):
    co,si=math.cos(-ang),math.sin(-ang)
    out=[]
    for x,y in P:
        x-=c[0]; y-=c[1]; x/=scale; y/=scale
        out.append((x*co-y*si,x*si+y*co))
    return out
def fmt(x):
    if abs(x)<5e-11: return '0'
    s=('%.10g'%x)
    return s
def trig_exact(k,N):
    """exact strings for cos(2πk/N), sin(2πk/N) when N | 12 or N=8 ... else symbolic"""
    k%=N
    table={}
    # angles in degrees with exact values
    ex={0:('1','0'),30:('√3/2','1/2'),45:('√2/2','√2/2'),60:('1/2','√3/2'),90:('0','1'),120:('−1/2','√3/2'),135:('−√2/2','√2/2'),150:('−√3/2','1/2'),180:('−1','0'),
        210:('−√3/2','−1/2'),225:('−√2/2','−√2/2'),240:('−1/2','−√3/2'),270:('0','−1'),300:('1/2','−√3/2'),315:('√2/2','−√2/2'),330:('√3/2','−1/2')}
    deg=360*k/N
    if abs(deg-round(deg))<1e-9 and int(round(deg)) in ex: return ex[int(round(deg))]
    g=math.gcd(k,N); kk,NN=k//g,N//g
    if kk==0: return ('1','0')
    a='2π·%d/%d'%(kk,NN) if kk!=1 else '2π/%d'%NN
    return ('cos(%s)'%a,'sin(%s)'%a)
def polygon_frame(P,Nhint=None):
    """find circle containing most points at angles that are multiples of 2π/N (N≤36)"""
    n=len(P); best=None
    m=min(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(n),2))
    Ns=[Nhint] if Nhint else list(range(3,37))
    for T in itertools.combinations(range(n),3):
        c=circumcenter(P[T[0]],P[T[1]],P[T[2]])
        if c is None: continue
        R=math.dist(c,P[T[0]])
        on=[i for i in range(n) if abs(math.dist(c,P[i])-R)<1e-7*m]
        if len(on)<3: continue
        angs=[math.atan2(P[i][1]-c[1],P[i][0]-c[0]) for i in on]
        for N in Ns:
            base=angs[0]; ok=True; ks=[]
            for a in angs:
                t=(a-base)/(2*math.pi/N); kr=round(t)
                if abs(t-kr)>1e-6: ok=False;break
                ks.append(kr%N)
            if ok:
                cand=(len(on),-N,c,R,base,on,ks)
                if best is None or cand[:2]>best[:2]: best=cand
    if best is None: return None
    cnt,negN,c,R,base,on,ks=best
    N=-negN
    # rotate so that a vertex is at angle 0 and, if possible, the set has a vertical mirror axis... keep simple: vertex 0 at angle 0
    Q=transform(P,c,R,base)
    exact=[None]*n
    for i,k in zip(on,ks): exact[i]=trig_exact(k,N)
    for i in range(n):
        if exact[i] is None and abs(Q[i][0])<1e-9 and abs(Q[i][1])<1e-9: exact[i]=('0','0')
    return {'frame':'regular %d-gon frame: circumradius 1, centre at the origin, a vertex at angle 0'%N,'N':N,'pts':Q,'exact':exact,'on':len(on)}
def lattice_frame(P):
    n=len(P); m=min(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(n),2))
    # choose a shortest edge as unit vector
    for i,j in itertools.combinations(range(n),2):
        if abs(math.dist(P[i],P[j])-m)<1e-9:
            base=(i,j); break
    i,j=base
    ang=math.atan2(P[j][1]-P[i][1],P[j][0]-P[i][0])
    for extra in (0.0,):
        Q=transform(P,P[i],m,ang+extra)
        # triangular: x = p/2, y = q*sqrt3/2, p ≡ q mod 2
        ok=True; ex=[]
        for x,y in Q:
            p=2*x; q=2*y/S3
            if abs(p-round(p))>1e-7 or abs(q-round(q))>1e-7 or (round(p)-round(q))%2: ok=False;break
            p=int(round(p)); q=int(round(q))
            sx=str(p//2) if p%2==0 else ('%d/2'%p)
            if q==0: sy='0'
            elif q%2==0: sy=('√3' if q//2==1 else ('−√3' if q//2==-1 else '%d√3'%(q//2)))
            else: sy=('√3/2' if q==1 else ('−√3/2' if q==-1 else '%d√3/2'%q))
            ex.append((sx,sy))
        if ok: return {'frame':'triangular lattice, unit spacing','pts':Q,'exact':ex}
        ok=True; ex=[]
        for x,y in Q:
            if abs(x-round(x))>1e-7 or abs(y-round(y))>1e-7: ok=False;break
            ex.append((str(int(round(x))),str(int(round(y)))))
        if ok: return {'frame':'square lattice, unit spacing','pts':Q,'exact':ex}
    return None
def collinear_frame(P):
    n=len(P)
    xs=[p[0] for p in P]; ys=[p[1] for p in P]
    # check collinearity
    (x0,y0),(x1,y1)=P[0],P[1]
    for x,y in P[2:]:
        if abs((x1-x0)*(y-y0)-(y1-y0)*(x-x0))>1e-9: return None
    m=min(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(n),2))
    ang=math.atan2(y1-y0,x1-x0)
    Q=transform(P,P[0],m,ang)
    xs=sorted(q[0] for q in Q); shift=xs[0]
    Q=[(x-shift,0.0) for x,_ in Q]
    ex=[]
    for x,_ in Q:
        f=Fraction(x).limit_denominator(1000)
        if abs(float(f)-x)<1e-9: ex.append((str(f),'0'))
        else: ex.append(None)
    return {'frame':'on a line, shortest distance 1','pts':Q,'exact':ex if all(ex) else None}
def exactify(P,Nhint=None):
    P=[tuple(map(float,p)) for p in P]
    col=collinear_frame(P)
    if col: return col
    lat=lattice_frame(P)
    if lat: return lat
    pf=polygon_frame(P,Nhint)
    if pf and pf['on']>=3 and (pf['N']<=36): return pf
    # generic frame: P0 origin, P1 on x-axis, shortest distance 1
    n=len(P); m=min(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(n),2))
    ang=math.atan2(P[1][1]-P[0][1],P[1][0]-P[0][0])
    Q=transform(P,P[0],m,ang)
    return {'frame':'shortest distance 1, first point at the origin, second on the positive x-axis','pts':Q,'exact':None}
if __name__=='__main__':
    import json
    hexa=[(math.cos(math.pi*k/3),math.sin(math.pi*k/3)) for k in range(6)]+[(0,0)]
    r=exactify(hexa); print(r['frame'],r['exact'])
    pent=[(math.cos(2*math.pi*k/5+0.3)*2.2+1,math.sin(2*math.pi*k/5+0.3)*2.2) for k in range(4)]+[(1,0)]
    r=exactify(pent); print(r['frame'],r['exact'])
    r=exactify([(0,0),(1,0),(2,0),(4,0)]); print(r['frame'],r['exact'])
