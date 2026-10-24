"""Known example configurations for cells whose full list is not known: subsets of regular polygons,
triangular-lattice pieces (hexagonal), and the Erdős–Fishburn 16-point 7-distance set."""
import math, itertools, json, random
S3=math.sqrt(3)
def ndist(P):
    ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2))
    vals=[ds[0]]
    for x in ds[1:]:
        if x-vals[-1]>1e-7*ds[-1]: vals.append(x)
    return vals
def poly(N,S,R=1.0):
    return [(R*math.cos(2*math.pi*k/N),R*math.sin(2*math.pi*k/N)) for k in S]
def chord_classes(N,S):
    return set(min((a-b)%N,(b-a)%N) for a,b in itertools.combinations(S,2))
def polygon_examples(n,d,limit=3):
    out=[]
    for N in (2*d,2*d+1):
        if n>N: continue
        if n==N:
            out.append(('regular %d-gon'%N,poly(N,range(N)),N,list(range(N)))); continue
        # search subsets of size n with exactly d chord classes: prefer "consecutive-ish" nice ones
        found=None
        cands=[]
        # consecutive block
        S=list(range(n))
        if len(chord_classes(N,S))==d: cands.append(S)
        # polygon minus a few vertices (remove N-n vertices spread out)
        rem=N-n
        for R_ in itertools.combinations(range(N),rem):
            S=[k for k in range(N) if k not in R_]
            if len(chord_classes(N,S))==d:
                cands.append(S)
                if len(cands)>=40: break
        # pick canonical: minimise lexicographic "gap pattern" up to rotation/reflection
        def key(S):
            best=None
            for r in range(N):
                for sgn in (1,-1):
                    T=sorted(((sgn*k+r)%N) for k in S)
                    if best is None or T<best: best=T
            return best
        seen=set(); picks=[]
        for S in cands:
            k=tuple(key(S))
            if k in seen: continue
            seen.add(k); picks.append(list(k))
        picks.sort(key=lambda S:(sum(S),S))
        for S in picks[:2]:
            rem_=[k for k in range(N) if k not in S]
            lab='regular %d-gon minus vertex%s %s'%(N,'' if len(rem_)==1 else 'es',', '.join(map(str,rem_))) if len(rem_)<=4 else 'vertices %s of a regular %d-gon'%(', '.join(map(str,S)),N)
            out.append((lab,poly(N,S),N,S))
    return out[:limit]
def tri(a,b): return (a+b/2,b*S3/2)
def hexagon_side(r):
    # lattice points with hex-norm <= r
    pts=[]
    for a in range(-r,r+1):
        for b in range(-r,r+1):
            c=-a-b
            if max(abs(a),abs(b),abs(c))<=r: pts.append((a,b))
    return pts
def lattice_examples(n,d):
    """a few lattice pieces with exactly d distances, found by greedy removal from hexagons"""
    out=[]
    random.seed(n*100+d)
    base=hexagon_side(2) if n<=19 else hexagon_side(3)
    tries=0; seen=set()
    while tries<4000 and len(out)<2:
        tries+=1
        S=random.sample(base,n)
        P=[tri(*p) for p in S]
        if len(ndist(P))==d:
            key=tuple(sorted(round(x,6) for x in [math.dist(P[i],P[j]) for i,j in itertools.combinations(range(n),2)]))
            if key in seen: continue
            seen.add(key); out.append(('piece of the triangular lattice',P,S))
    return out
def ef16():
    # hexagon of side 2 (19 points) minus one endpoint of each of the three long diameters -> 16 points, 7 distances
    pts=hexagon_side(2)
    outer=[(a,b) for a,b in pts if max(abs(a),abs(b),abs(a+b))==2]
    # remove (2,0),(0,2),(-2,2): one from each antipodal pair among the 6 corners
    rem={(2,0),(0,2),(-2,2)}
    S=[p for p in pts if p not in rem]
    P=[tri(*p) for p in S]
    assert len(P)==16 and len(ndist(P))==7, (len(P),len(ndist(P)))
    return P,S
if __name__=='__main__':
    P,S=ef16(); print('EF16 ok',[round(v,4) for v in ndist(P)])
    for n,d in [(8,4),(7,4),(6,4),(10,5),(11,5),(12,6),(14,7),(16,7),(6,7),(7,7)]:
        ex=polygon_examples(n,d)
        print(n,d,[(e[0],len(ndist(e[1]))) for e in ex])
        le=lattice_examples(n,d); print('   lattice',len(le))
