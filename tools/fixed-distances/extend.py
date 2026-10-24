"""Complete one-point extensions: given an n-point d-distance set X, find every point P such that X+P has
exactly d or d+1 distances (P's distances lie in the old values plus at most one new value)."""
import math, itertools
def dvals(P,tol=1e-7):
    ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)); v=[ds[0]]
    for x in ds[1:]:
        if x-v[-1]>tol*ds[-1]: v.append(x)
    return v
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
def extensions(X,allow_new=True):
    vals=dvals(X); D=max(vals); n=len(X)
    cands=[]
    for i,j in itertools.combinations(range(n),2):
        for ri in vals:
            for rj in vals: cands+=circ2(X[i],ri,X[j],rj)
    if allow_new:
        for T in itertools.combinations(range(n),3):
            c=circumcenter(X[T[0]],X[T[1]],X[T[2]])
            if c: cands.append(c)
    out=[]
    for P in cands:
        if min(math.dist(P,x) for x in X)<1e-6*D: continue
        new=[]
        ok=True
        for x in X:
            dd=math.dist(P,x)
            if min(abs(dd-v) for v in vals)<1e-6*D: continue
            new.append(dd)
        if new:
            if not allow_new or max(new)-min(new)>1e-6*D: continue
        Y=X+[P]
        if any(same(Y,Z) for Z in out): continue
        out.append(Y)
    return out
def simkey(P):
    ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)); m=ds[0]
    return tuple(round(x/m,6) for x in ds)

def same(P,Q,tol=1e-6):
    if len(P)!=len(Q): return False
    a=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2))
    b=sorted(math.dist(Q[i],Q[j]) for i,j in itertools.combinations(range(len(Q)),2))
    s=b[0]/a[0]
    return all(abs(x*s-y)<tol*b[-1] for x,y in zip(a,b))
