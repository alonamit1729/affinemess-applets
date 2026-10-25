"""Reconstruct the 42 seven-point four-distance sets from Lan–Wei (2013), Fig. 4 and the proof of Theorem 8:
 701-703, 708: R9 minus two vertices (four dissimilar versions); 707: R8 minus a vertex; 704: R7+ minus a vertex;
 705, 709-716: nine 7-subsets of the 9-point set of three concentric triangles (Erdős–Fishburn's third 9-set);
 706, 721-739: twenty pieces of the triangular lattice; 717-718: 7-subsets of the square with four apexes;
 719-720: 7-subsets of the 10-point 5-distance set (R5 doubled about its centre);
 740: pattern 1:√2:2:√5; 741: 1:√(2+√3):√(4+2√3):√(5+2√3); 742: 1:√2:√(2+√3):√(4+√3)."""
import math, itertools, json, numpy as np, networkx as nx
s3=math.sqrt(3); phi=(1+math.sqrt(5))/2
def poly(n,center=False):
    P=[(math.cos(2*math.pi*k/n),math.sin(2*math.pi*k/n)) for k in range(n)]
    return P+([(0.0,0.0)] if center else [])
def dvals(P):
    D=sorted(math.dist(a,b) for a,b in itertools.combinations(P,2)); D=[x/D[0] for x in D]; cl=[D[0]]
    for x in D[1:]:
        if x-cl[-1]>1e-7: cl.append(x)
    return cl
def key(P):
    D=np.array([[math.dist(a,b) for b in P] for a in P]); D/=D[D>0].min()
    return tuple(sorted(tuple(np.round(sorted(r),6)) for r in D))
def congruent(P,Q):
    if key(P)!=key(Q): return False
    DP=np.array([[math.dist(a,b) for b in P] for a in P]); DQ=np.array([[math.dist(a,b) for b in Q] for a in Q])
    DP/=DP[DP>0].min(); DQ/=DQ[DQ>0].min(); n=len(P)
    for p in itertools.permutations(range(n)):
        if all(abs(DP[i][j]-DQ[p[i]][p[j]])<1e-6 for i in range(n) for j in range(i+1,n)): return True
    return False
found=[]
def add(label,P,src):
    if len(P)!=7 or len(dvals(P))!=4: return False
    for f in found:
        if congruent(P,f['pts']): return False
    found.append({'label':label,'pts':[(float(x),float(y)) for x,y in P],'src':src}); return True
def subsets(P,label,src):
    n=0
    for drop in itertools.combinations(range(len(P)),len(P)-7):
        Q=[p for i,p in enumerate(P) if i not in drop]
        if add(f'{label} minus {len(drop)} point(s)',Q,src): n+=1
    return n
# polygons
print('R9-2:',subsets(poly(9),'regular nonagon','poly'))
print('R8-1:',subsets(poly(8),'regular octagon','poly'))
print('R7+-1:',subsets(poly(7,True),'regular heptagon with centre','poly'))
# three concentric triangles (Erdős–Fishburn)
FOUR=[(1,0),(0,1),(-1,0),(0,-1),(0,-s3),((s3-1)/2,-(s3-1)/2),(-(s3-1)/2,-(s3-1)/2),((s3+1)/2,-(s3+1)/2),(-(s3+1)/2,-(s3+1)/2)]
print('three triangles:',subsets(FOUR,'three concentric equilateral triangles (9-point set)','ef9'))
# square with four apexes (Shinohara fig 1(e))
h=.5+s3/2; FIGE=[(.5,.5),(-.5,.5),(-.5,-.5),(.5,-.5),(0,h),(0,-h),(h,0),(-h,0)]
print('square+apexes:',subsets(FIGE,'square with four outward equilateral apexes (8-point set)','sh8'))
# double pentagon: R5 and the inner pentagon of its pentagram (concentric, rotated by pi/5, ratio phi^2): 10 points, 5 distances
P10=poly(5)+[(math.cos(2*math.pi*k/5+math.pi/5)/phi**2,math.sin(2*math.pi*k/5+math.pi/5)/phi**2) for k in range(5)]
assert len(dvals(P10))==5
print('double pentagon:',subsets(P10,'the 10-point five-distance set (a regular pentagon and the inner pentagon of its pentagram)','dp10'))
# lattice pieces: triangular lattice, all 4-subsets of squared distances
def lattice_search(pts,d2,sqs,label,src):
    n=0
    for D4 in itertools.combinations(sqs,4):
        G=nx.Graph(); G.add_nodes_from(range(len(pts)))
        for i,j in itertools.combinations(range(len(pts)),2):
            if d2(pts[i],pts[j]) in D4: G.add_edge(i,j)
        for c in nx.find_cliques(G):
            if len(c)<7: continue
            for sub in itertools.combinations(c,7):
                P=[pts[i] for i in sub]
                if len({d2(a,b) for a,b in itertools.combinations(P,2)})==4 and add(label,P,src): n+=1
    return n
tri=[(a+b/2,b*s3/2) for a in range(-6,7) for b in range(-6,7)]; tri=[p for p in tri if math.hypot(*p)<=4.6]
print('triangular lattice:',lattice_search(tri,lambda p,q: round((p[0]-q[0])**2+(p[1]-q[1])**2),[1,3,4,7,9,12,13,16,19,21],'piece of the triangular lattice','lat'))
sq=[(a,b) for a in range(-5,6) for b in range(-5,6)]
print('square lattice:',lattice_search(sq,lambda p,q: (p[0]-q[0])**2+(p[1]-q[1])**2,[1,2,4,5,8,9,10,13],'piece of the square lattice','sqlat'))
# prescribed distance sets (741, 742): backtracking construction
def build_with(vals,label,src,maxpts=7):
    vals=sorted(vals); n=0
    def circ(c1,r1,c2,r2):
        dx,dy=c2[0]-c1[0],c2[1]-c1[1]; dd=math.hypot(dx,dy)
        if dd<1e-12 or dd>r1+r2+1e-12 or dd<abs(r1-r2)-1e-12: return []
        a=(r1*r1-r2*r2+dd*dd)/(2*dd); hh=math.sqrt(max(r1*r1-a*a,0)); mx,my=c1[0]+a*dx/dd,c1[1]+a*dy/dd
        return [(mx+hh*dy/dd,my-hh*dx/dd),(mx-hh*dy/dd,my+hh*dx/dd)]
    def ok(P,q):
        return all(any(abs(math.dist(p,q)-v)<1e-9 for v in vals) for p in P) and all(math.dist(p,q)>1e-9 for p in P)
    seen=set()
    def rec(P):
        nonlocal n
        if len(P)==maxpts:
            if len({round(math.dist(a,b),6) for a,b in itertools.combinations(P,2)})==4 and add(label,P,src): n+=1
            return
        cand=[]
        for i,j in itertools.combinations(range(len(P)),2):
            for r1 in vals:
                for r2 in vals:
                    for q in circ(P[i],r1,P[j],r2):
                        if ok(P,q): cand.append(q)
        uniq=[]
        for q in cand:
            if not any(math.dist(q,u)<1e-7 for u in uniq): uniq.append(q)
        for q in uniq:
            Q=P+[q]; k=key(Q)
            if k in seen: continue
            seen.add(k); rec(Q)
    for d1 in vals: rec([(0.0,0.0),(d1,0.0)])
    return n
print('741:',build_with([1,math.sqrt(2+s3),math.sqrt(4+2*s3),math.sqrt(5+2*s3)],'distances 1 : √(2+√3) : 1+√3 : √(5+2√3)','lw741'))
print('742:',build_with([1,math.sqrt(2),math.sqrt(2+s3),math.sqrt(4+s3)],'distances 1 : √2 : √(2+√3) : √(4+√3)','lw742'))
print('TOTAL',len(found))
json.dump(found,open('seven4.json','w'))
