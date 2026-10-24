import json, math, itertools
E5=list(itertools.combinations(range(5),2))
def pdist(P): return [math.dist(P[i],P[j]) for i,j in E5]
def dm(P): return [[math.dist(P[i],P[j]) for j in range(5)] for i in range(5)]
def congruent(P,Q,tol=1e-6):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(5)):
        s=D2[p[0]][p[1]]/D1[0][1]
        if all(abs(D1[i][j]*s-D2[p[i]][p[j]])<tol for i in range(5) for j in range(i+1,5)): return True
    return False
A=[x['pts'] for x in json.load(open('five.json'))]
B=[x['pts'] for x in json.load(open('five_brute.json')) if min(pdist(x['pts']))>1e-3]
print(len(A),len(B))
allc=list(A)
for Q in B:
    if not any(congruent(P,Q) for P in allc): allc.append(Q); print('brute-only',sorted(pdist(Q)))
for P in A:
    if not any(congruent(P,Q) for Q in B): print('ext-only',sorted(pdist(P)))
print('merged',len(allc))
def symorder(P):
    D=dm(P); n=0
    for p in itertools.permutations(range(5)):
        if all(abs(D[i][j]-D[p[i]][p[j]])<1e-6 for i in range(5) for j in range(i+1,5)): n+=1
    return n
def circumcenter(p,q,r):
    ax,ay=p;bx,by=q;cx,cy=r
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-9: return None
    ux=((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d
    uy=((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d
    return (ux,uy)
def info(P):
    ds=pdist(P); m=min(ds); vals=sorted(set(round(d/m,5) for d in ds))
    cnt=[sum(1 for d in ds if abs(d/m-v)<1e-4) for v in vals]
    # collinear triples
    col=[T for T in itertools.combinations(range(5),3) if abs((P[T[1]][0]-P[T[0]][0])*(P[T[2]][1]-P[T[0]][1])-(P[T[1]][1]-P[T[0]][1])*(P[T[2]][0]-P[T[0]][0]))<1e-7*m*m]
    # largest concyclic subset
    best=None
    for S in (list(itertools.combinations(range(5),5))+list(itertools.combinations(range(5),4))):
        c=circumcenter(P[S[0]],P[S[1]],P[S[2]])
        if c and all(abs(math.dist(c,P[i])-math.dist(c,P[S[0]]))<1e-7*m for i in S):
            best=(S,c); break
    ang=None
    if best:
        S,c=best; angs=sorted((math.degrees(math.atan2(P[i][1]-c[1],P[i][0]-c[0]))%360) for i in S)
        gaps=[round(((angs[(i+1)%len(angs)]-angs[i])%360),3) for i in range(len(angs))]
        R=math.dist(c,P[S[0]])/m
        others=[i for i in range(5) if i not in S]
        od=[round(math.dist(c,P[i])/m,4) for i in others]
        ang=(len(S),gaps,round(R,4),od)
    return vals,cnt,symorder(P),col,ang
out=[]
for k,P in enumerate(allc):
    vals,cnt,so,col,ang=info(P)
    print(k+1,'ratios',vals,'mult',cnt,'sym',so,'collinear',col,'concyclic',ang)
    out.append({'pts':P,'ratios':vals,'mult':cnt,'sym':so,'collinear':col,'concyclic':ang})
json.dump(out,open('five_all.json','w'))
