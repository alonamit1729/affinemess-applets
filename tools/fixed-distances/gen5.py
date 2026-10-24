import json, math, itertools
import numpy as np
five=[d['pts'] for d in json.load(open('five_all.json'))]
six=json.load(open('six.json'))
S3=math.sqrt(3)
def cc(p,q,r):
    ax,ay=p;bx,by=q;cx,cy=r
    dd=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(dd)<1e-9: return None
    return (((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/dd,((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/dd)
def pd(P): return [math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)]
def find_circle(P,n):
    # largest concyclic subset whose angles are multiples of 360/n; return centre, R, base angle (deg)
    m=min(pd(P)); best=None
    for k in range(len(P),2,-1):
        for S in itertools.combinations(range(len(P)),k):
            o=cc(P[S[0]],P[S[1]],P[S[2]])
            if not o: continue
            R=math.dist(o,P[S[0]])
            if not all(abs(math.dist(o,P[i])-R)<1e-6*m for i in S): continue
            if k==3 and min(math.dist(o,p) for p in P)>1e-6*m: continue
            angs=[math.degrees(math.atan2(P[i][1]-o[1],P[i][0]-o[0])) for i in S]
            base=angs[0]
            if all(abs(((a-base)/(360/n))-round((a-base)/(360/n)))<1e-4 for a in angs):
                if best is None or k>best[0]: best=(k,o,R,base)
        if best: break
    return best
def symmetry_axis(P):
    n=len(P); D=[[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
    X=np.array(P); c=X.mean(0); X0=X-c
    for p in itertools.permutations(range(n)):
        if list(p)==list(range(n)): continue
        if not all(abs(D[i][j]-D[p[i]][p[j]])<1e-6 for i in range(n) for j in range(i+1,n)): continue
        Y=X0[list(p)]
        H=X0.T@Y; U,s,Vt=np.linalg.svd(H); Rm=(U@Vt).T
        if np.linalg.det(Rm)<0:
            w,v=np.linalg.eigh(Rm); ax=v[:,np.argmax(w)]
            return c,ax
    return c,None
def prep(P,n,ghost_extra=None,rot_hint=None):
    P=[tuple(p) for p in P]
    c,ax=symmetry_axis(P)
    if ax is not None:
        ang=math.atan2(ax[1],ax[0])-math.pi/2
    else:
        ang=0.0
    circ=find_circle(P,n) if n else None
    if ax is None and circ: ang=math.radians(circ[3])-math.pi/2  # put a polygon vertex at top
    def rot(p):
        x,y=p[0]-c[0],p[1]-c[1]; co,si=math.cos(-ang),math.sin(-ang); return (x*co-y*si, x*si+y*co)
    Q=[rot(p) for p in P]
    ghost=[]; circle=None
    if circ:
        k,o,R,base=circ; oc=rot(o); circle=(oc,R)
        for j in range(n):
            a=math.radians(base)+2*math.pi*j/n-ang
            g=(oc[0]+R*math.cos(a),oc[1]+R*math.sin(a))
            if min(math.dist(g,q) for q in Q)>1e-6: ghost.append(g)
    if ghost_extra:
        ghost+= [rot(g) for g in ghost_extra]
    # flip so that more mass is at the bottom? keep as is. normalise
    xs=[q[0] for q in Q]+[g[0] for g in ghost]; ys=[q[1] for q in Q]+[g[1] for g in ghost]
    if circle:
        xs+=[circle[0][0]-circle[1],circle[0][0]+circle[1]]; ys+=[circle[0][1]-circle[1],circle[0][1]+circle[1]]
    cx=(max(xs)+min(xs))/2; cy=(max(ys)+min(ys))/2; sc=max(max(xs)-min(xs),max(ys)-min(ys))
    N=lambda p:(round((p[0]-cx)/sc,5),round((p[1]-cy)/sc,5))
    ds=pd(P); vals=sorted(set(round(d,7) for d in ds)); vals=[v for i,v in enumerate(vals) if i==0 or v-vals[i-1]>1e-6]
    cls=[min(range(3),key=lambda k:abs(vals[k]-d)) for d in ds]
    return {'pts':[N(q) for q in Q],'ghost':[N(g) for g in ghost],'circle':[N(circle[0]),round(circle[1]/sc,5)] if circle else None,'cls':cls,'ratio':[round(v/vals[0],3) for v in vals]}
def tri_ghost(P):
    # ghost = missing midpoint of the side-2 triangle
    m=min(pd(P)); T=[S for S in itertools.combinations(range(5),3) if all(abs(math.dist(P[a],P[b])-2*m)<1e-6 for a,b in itertools.combinations(S,2))][0]
    mids=[((P[a][0]+P[b][0])/2,(P[a][1]+P[b][1])/2) for a,b in itertools.combinations(T,2)]
    return [g for g in mids if min(math.dist(g,p) for p in P)>1e-6]
F=lambda k:five[k-1]
groups=[
 ('Square','Four vertices of a square plus a point on its axis of symmetry. Side 1.',4,[
  (5,'the centre','1/√2 : 1 : √2'),(1,'inner apex of the equilateral triangle on a side','√(2−√3) : 1 : √2'),(2,'outer apex of that triangle','1 : √2 : √(2+√3)'),
  (3,'axis point at distance √2 from the far vertices','√(3−√7) : 1 : √2'),(4,'axis point at distance √2 from the near vertices','1 : √2 : √(3+√7)')]),
 ('Regular pentagon','Four vertices of a regular pentagon plus one point. The missing vertex is “the gap”, the side opposite it “the far side”. Side 1, diagonal φ, circumradius R ≈ 0.851.',5,[
  (14,'the centre','1 : φ : R'),(8,'crossing point of the two diagonals','1/φ : 1 : φ'),(13,'pentagram tip over the far side','1 : φ : φ²'),
  (11,'inner equilateral apex on the far side','0.813 : 1 : φ'),(12,'outer equilateral apex on the far side','1 : φ : 1.989'),
  (9,'outer equilateral apex on the open diagonal','1 : φ : 2.405'),(10,'axis point beyond the far side, at diagonal-distance from the gap vertices','0.673 : 1 : φ')]),
 ('Triangular lattice','Five points of the triangular lattice; every one has distances 1 : √3 : 2.',6,[
  (30,'regular hexagon minus a vertex','1 : √3 : 2'),(6,'hexagon + centre minus two adjacent vertices','1 : √3 : 2'),(7,'hexagon + centre minus two vertices √3 apart','1 : √3 : 2'),
  (34,'hexagon + centre minus two opposite vertices','1 : √3 : 2'),(26,'side‑2 lattice triangle minus a side midpoint','1 : √3 : 2')]),
 ('Regular heptagon','The regular heptagon minus two vertices, in the three inequivalent ways. Side 1, diagonals 1.802 and 2.247.',7,[
  (28,'adjacent vertices removed','1 : 1.802 : 2.247'),(27,'vertices two apart removed','1 : 1.802 : 2.247'),(25,'vertices three apart removed','1 : 1.802 : 2.247')]),
 ('Regular 12‑gon','Points of a regular dodecagon with centre O and consecutive vertices V₀, V₁, …; side 1, circumradius R = (√6+√2)/2 ≈ 1.932.',12,[
  (17,'O + V₀V₁V₂V₃','1 : R : 1+√3'),(20,'V₀V₁V₂V₃ + outer equilateral apex on V₁V₂','1 : R : 1+√3'),(29,'V₀V₁V₂V₃ + reflection of O in V₁V₂','1 : R : 1+√3'),
  (18,'O + V₀V₁V₂ + outer equilateral apex on V₀V₁','1 : R : 1+√3'),(16,'O + V₀V₁V₂ + inner equilateral apex on V₀V₁, which lands at distance 1 from O','1 : √2 : R'),
  (19,'O + V₀V₂V₅V₇','R : 1+√3 : 2+√3'),(15,'V₀V₂V₇V₉ + fourth vertex of the square O V₀ · V₉','R : 1+√3 : 2+√3'),(24,'V₀V₃V₅V₁₀ + reflection of O in the chord V₅V₁₀','R : 1+√3 : 2+√3')]),
 ('Regular 18‑gon','Four vertices of a regular 18‑gon plus its centre O; circumradius 1, chords cₖ = 2 sin 10k°, so c₃ = 1.',18,[
  (23,'O + four consecutive vertices','c₁ : c₂ : 1'),(21,'O + V₀V₄V₇V₁₁','1 : c₄ : c₇'),(22,'O + V₀V₅V₈V₁₃','1 : c₅ : c₈')]),
 ('Lattice pieces','The points 0, u, v, u+v, u−v of a lattice; the parallelogram 0, u, v, u+v is one of the 2 : 1 parallelograms, completed by the reflection of v through u/2.',None,[
  (31,'|u| = 2, |v| = 1, |u+v| = √6','1 : 2 : √6'),(32,'|u| = √2, |v| = 1, |u+v| = 2','1 : √2 : 2')]),
 ('One loner','A rectangle with sides 1 and w = √((1+√13)/2) ≈ 1.517 and the apex of the isosceles triangle with legs 1 erected on a long side; w is exactly what makes the apex sit at diagonal‑distance from the far vertices.',None,[
  (33,'rectangle + isosceles apex','1 : w : √((3+√13)/2)')]),
]
out=[]; num=0
for title,desc,n,items in groups:
    g={'title':title,'desc':desc,'items':[]}
    for k,label,dist in items:
        num+=1
        extra=tri_ghost(F(k)) if k==26 else None
        it=prep(F(k),None if k==26 else n,extra); it.update({'num':num,'label':label,'dist':dist}); g['items'].append(it)
    out.append(g)
# six- and seven-point sets
sixdefs=[(5,'regular pentagon + centre',5,'1 : φ : R'),(8,'regular heptagon minus a vertex',7,'1 : 1.802 : 2.247'),(9,'regular hexagon',6,'1 : √3 : 2'),(4,'hexagon + centre minus a vertex',6,'1 : √3 : 2'),(3,'the side‑2 lattice triangle',None,'1 : √3 : 2'),
 (1,'12‑gon: O + V₀V₁V₂V₃ + outer apex on V₁V₂',12,'1 : R : 1+√3'),(7,'12‑gon: V₀V₁V₂V₃ + outer apex on V₁V₂ + reflection of O in V₁V₂',12,'1 : R : 1+√3'),
 (2,'12‑gon: O + V₀V₂V₇V₉ + fourth vertex of the square O V₀ · V₉',12,'R : 1+√3 : 2+√3'),(6,'12‑gon: V₀V₃V₅V₁₀ + reflection of O in V₅V₁₀ + fourth vertex of the square O V₀ · V₃',12,'R : 1+√3 : 2+√3')]
g={'title':'Six points','desc':'The nine six‑point three‑distance sets. Each contains several of the 34 above; the 12‑gon ones are unions of neighbouring entries.','items':[]}
for k,label,n,dist in sixdefs:
    it=prep(six[k-1],n); it.update({'num':'','label':label,'dist':dist}); g['items'].append(it)
out.append(g)
hept=[(math.cos(2*math.pi*j/7),math.sin(2*math.pi*j/7)) for j in range(7)]
hexc=[(0,0)]+[(math.cos(2*math.pi*j/6),math.sin(2*math.pi*j/6)) for j in range(6)]
g={'title':'Seven points','desc':'The maximum. Erdős–Fishburn showed no three‑distance set in the plane has more than seven points; these are the only two with seven.','items':[]}
for P,label,n,dist in [(hept,'regular heptagon',7,'1 : 1.802 : 2.247'),(hexc,'regular hexagon + centre',6,'1 : √3 : 2')]:
    it=prep(P,n); it.update({'num':'','label':label,'dist':dist}); g['items'].append(it)
out.append(g)
json.dump(out,open('five_page.json','w'))
print(num, sum(len(g['items']) for g in out))
