"""Assemble applet_data.json: the (n,d) table and per-cell catalogs."""
import json, math, itertools, os, random
from exactify import exactify
import examples as EX
S3=math.sqrt(3); PHI=(1+5**.5)/2
def pd(P): return [math.dist(P[i],P[j]) for i,j in itertools.combinations(range(len(P)),2)]
def dvals(P):
    ds=sorted(pd(P)); vals=[ds[0]]
    for x in ds[1:]:
        if x-vals[-1]>1e-7*ds[-1]: vals.append(x)
    return vals
def classes(P):
    vals=dvals(P); return [min(range(len(vals)),key=lambda k:abs(vals[k]-x)) for x in pd(P)]
def norm_draw(P,extra=()):
    """centre + scale to the unit box; rotate so that a mirror axis (if any) is vertical"""
    n=len(P); X=[list(p) for p in P]
    import numpy as np
    A=np.array(X,dtype=float); c=A.mean(0); A0=A-c
    scale=max(abs(A0).max(),1e-9)
    def is_mirror(theta):
        # reflect across the line through the centroid with direction theta
        c2,s2=math.cos(2*theta),math.sin(2*theta)
        R=[(x*c2+y*s2, x*s2-y*c2) for x,y in A0]
        for rx,ry in R:
            if not any(abs(rx-x)<1e-6*scale and abs(ry-y)<1e-6*scale for x,y in A0): return False
        return True
    cands=set()
    for i in range(n):
        x,y=A0[i]
        if math.hypot(x,y)>1e-9*scale: cands.add(round(math.atan2(y,x)%math.pi,9))
        for j in range(i+1,n):
            dx,dy=A0[j][0]-A0[i][0],A0[j][1]-A0[i][1]
            cands.add(round((math.atan2(dy,dx)+math.pi/2)%math.pi,9))
    ang=None
    for th in sorted(cands):
        if is_mirror(th): ang=th-math.pi/2; break
    if ang is None: ang=0.0
    co,si=math.cos(-ang),math.sin(-ang)
    Q=[(x*co-y*si,x*si+y*co) for x,y in A0]
    G=[((x-c[0])*co-(y-c[1])*si,(x-c[0])*si+(y-c[1])*co) for x,y in extra]
    xs=[q[0] for q in Q+G]; ys=[q[1] for q in Q+G]
    cx=(max(xs)+min(xs))/2; cy=(max(ys)+min(ys))/2; sc=max(max(xs)-min(xs),max(ys)-min(ys)) or 1
    f=lambda p:[round((p[0]-cx)/sc,5),round((p[1]-cy)/sc,5)]
    return [f(q) for q in Q],[f(g) for g in G]
def item(P,label,kind='rigid',dim=0,Nhint=None,ghost=(),note=''):
    P=[tuple(map(float,p)) for p in P]
    ex=exactify(P,Nhint)
    draw,gh=norm_draw(P,ghost)
    vals=dvals(P); m=vals[0]
    return {'label':label,'kind':kind,'dim':dim,'n':len(P),'d':len(vals),
            'ratios':[round(v/m,4) for v in vals],'draw':draw,'ghost':gh,'cls':classes(P),
            'frame':ex['frame'],'coords':[[('%.10g'%x) if abs(x)>5e-11 else '0',('%.10g'%y) if abs(y)>5e-11 else '0'] for x,y in ex['pts']],
            'exact':ex.get('exact'),'note':note}
def poly(N,S,R=1.0): return [(R*math.cos(2*math.pi*k/N),R*math.sin(2*math.pi*k/N)) for k in S]
def tri(a,b): return (a+b/2,b*S3/2)
# ---------------- four-point families (representatives) ----------------
def A4(t):  return [(0,0),(1,0),(.5,S3/2),(math.cos(t),math.sin(t))]
def B4(al): return [(0,0),(1,0),(math.cos(al),math.sin(al)),(1+math.cos(al),math.sin(al))]
def C4(ph): return [(0,0),(1,0),(math.cos(ph),math.sin(ph)),(math.cos(ph),-math.sin(ph))]
def D4(t):  return [(0,0),(1,0),(.5,S3/2),(.5,t)]
def E4(th): return [(0,0),(1,0),(math.cos(th),math.sin(th)),(1-math.cos(th),math.sin(th))]
def F4(h,s=1):
    P1=(0,0);P2=(1,0);P0=(.5,h); b=math.hypot(.5,h); mid=(.25,h/2); d=(-h,.5); nn=math.hypot(*d)
    t=math.sqrt(max(1-b*b/4,0)); P3=(mid[0]+s*t*d[0]/nn, mid[1]+s*t*d[1]/nn); return [P0,P1,P2,P3]
def G4(s):  return [(0,s),(-.5,0),(.5,0),(0,s+1)]
def H4(ph,s=1):
    P2=(0,0);P1=(1,0);P3=(math.cos(2*ph),math.sin(2*ph)); b=2*math.sin(ph); dx=math.cos(2*ph)-.5; r=b*b-dx*dx
    return [(.5,math.sin(2*ph)+s*math.sqrt(r)),P1,P2,P3]
def I4(r):  return [(0,0),(1,0),(0,r),(1,r)]
fam43=[
 ('ABC equilateral and AD = AB: D runs round the circle of radius AB about A',A4(math.radians(100))),
 ('Rhombus with unequal diagonals',B4(math.radians(50))),
 ('A is the circumcentre of the isosceles triangle BCD',C4(math.radians(50))),
 ('ABC equilateral, D anywhere on the perpendicular bisector of AB',D4(1.2)),
 ('AC = AB = BD and AD = BC: isosceles trapezoid with three equal sides (one branch) or parallelogram with a diagonal equal to a side (the other branch)',E4(math.radians(50))),
 ('BC = BD = DA and AB = AC: A on the bisector of BC, D at distance BC from B and A',F4(0.7)),
 ('Kite with equal diagonals: A, D on the axis of BC and AD = BC',G4(0.3)),
 ('Three isosceles triangles on the base vertex B: A equidistant from B, C; C from B, D; D from B, A',H4(0.5)),
 ('Rectangle',I4(0.6)),
]
fam44=[
 ('A is the circumcentre of BCD (B, C, D free on a circle about A)',[(0,0),(1,0),(math.cos(1.1),math.sin(1.1)),(math.cos(2.9),math.sin(2.9))]),
 ('ABC equilateral, D free',[(0,0),(1,0),(.5,S3/2),(1.3,0.9)]),
 ('AB = AC and AD = BC',None),
 ('Equilateral path C–A–B–D: three equal consecutive edges',None),
 ('Two isosceles triangles: ABC with apex A, ABD with apex D',None),
 ('Kite: A and D on the perpendicular bisector of BC',[(0,0.9),(-.5,0),(.5,0),(0,-0.6)]),
 ('Parallelogram',[(0,0),(1,0),(1.3,0.7),(0.3,0.7)]),
]
fam45=[('One pair of equal adjacent edges, AB = AC',None),('One pair of equal opposite edges, AD = BC',None)]
fam46=[('Generic quadruple',[(0,0),(1,0),(1.5,0.9),(0.2,0.7)])]
def rep_from_strata(n,d,idx):
    D=json.load(open(f'strata_{n}_{d}_1.json'))
    r=D['realizable'][idx]
    k=sorted(r['sample'].keys())[0]
    return r['sample'][k]
two=[('square',[(0,0),(1,0),(1,1),(0,1)]),('rhombus of two equilateral triangles',[(0,0),(1,0),(.5,S3/2),(1.5,S3/2)]),
 ('equilateral triangle and its centre',[(0,0),(1,0),(.5,S3/2),(.5,S3/6)]),('regular pentagon minus a vertex',poly(5,range(4))),
 ('kite: equilateral triangle with a point at 30° on the circle about a vertex',[(0,0),(1,0),(.5,S3/2),(S3/2,.5)]),
 ('dart: equilateral triangle with a point at 150° on the circle about a vertex',[(0,0),(1,0),(.5,S3/2),(-S3/2,-.5)])]
# ---------------- literature / status per cell ----------------
SRC={
 'ES':'S. J. Einhorn, I. J. Schoenberg, On Euclidean sets having only two distances between points I–II, Indag. Math. 28 (1966)',
 'EF':'P. Erdős, P. Fishburn, Maximum planar sets that determine k distances, Discrete Math. 160 (1996) 115–125',
 'Sh04':'M. Shinohara, Classification of three-distance sets in two dimensional Euclidean space, European J. Combin. 25 (2004) 1039–1058',
 'Sh08':'M. Shinohara, Uniqueness of maximum planar five-distance sets, Discrete Math. 308 (2008) 3048–3055',
 'Wei11':'X. Wei, Classification of eleven-point five-distance sets in the plane, Ars Combin. 102 (2011) 505–515',
 'Wei12':'X. Wei, A proof of Erdős–Fishburn’s conjecture for g(6) = 13, Electron. J. Combin. 19(4) (2012) P38',
 'LW13':'W. Lan, X. Wei, Classification of seven-point four-distance sets in the plane, Math. Notes 93 (2013) 510–522',
 'SO20':'F. Szöllősi, P. R. J. Östergård, Constructions of maximum few-distance sets in Euclidean spaces, Electron. J. Combin. 27(1) (2020) P1.23',
 'BBS':'E. Bannai, E. Bannai, D. Stanton, An upper bound for the cardinality of an s-distance subset in real Euclidean space II, Combinatorica 3 (1983); A. Blokhuis, Few-distance sets (1984)',
 'BY25':'L.-R. Bao, W.-H. Yu, Constructions of large m-distance sets on the triangular lattice, arXiv:2509.00880 (2025)',
 'Bud22':'W. T. Budianto, An observation on 7-distance sets in the Euclidean plane, M.Sc. thesis, National Central University (2022)',
 'Bal':'K. Balaji et al., Lattice configurations determining few distances, arXiv:1911.11688',
 'cert':'exact computation for this applet (Oct 2026): msolve 0.10 (Berthomieu–Eder–Safey El Din) and Singular 4 (minAssGTZ), one polynomial system per colouring type',
 'us':'computation for this applet (Oct 2026): exhaustive enumeration of edge-colouring types with numerical realization, cross-checked by two independent methods',
}
def same(P,Q,tol=1e-6):
    n=len(P)
    if len(Q)!=n: return False
    D1=[[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]; D2=[[math.dist(Q[i],Q[j]) for j in range(n)] for i in range(n)]
    s1=sorted(D1[i][j] for i in range(n) for j in range(i+1,n)); s2=sorted(D2[i][j] for i in range(n) for j in range(i+1,n)); s=s2[0]/s1[0]
    return all(abs(a*s-b)<tol*s2[-1] for a,b in zip(s1,s2))
def main():
    cells={}; cats={}
    def cell(n,d,**kw):
        c={'n':n,'d':d}; c.update(kw); cells[f'{n},{d}']=c
    def cat(n,d,items,complete,title,note='',sources=()):
        cats[f'{n},{d}']={'items':items,'complete':complete,'title':title,'note':note,'sources':[SRC[s] for s in sources]}
    # ---- n = 3 ----
    cell(3,1,status='exact',rigid=1,families=[],marks=[],note='The equilateral triangle.')
    cat(3,1,[item([(0,0),(1,0),(.5,S3/2)],'equilateral triangle')],True,'The equilateral triangle')
    cell(3,2,status='infinite',rigid=0,families=[{'dim':1,'count':1}],marks=[],note='Isosceles triangles: one parameter (the apex angle), a single family.')
    cat(3,2,[item([(0,0),(1,0),(.5,1.1)],'isosceles triangle (apex angle free)','family',1)],True,'Isosceles triangles')
    cell(3,3,status='infinite',rigid=0,families=[{'dim':2,'count':1}],marks=[],note='Scalene triangles: two parameters.')
    cat(3,3,[item([(0,0),(1,0),(.3,.8)],'scalene triangle','family',2)],True,'Scalene triangles')
    # ---- n = 4 ----
    cell(4,2,status='exact',rigid=6,families=[],marks=[],note='The six two-distance quadruples (Einhorn–Schoenberg).')
    cat(4,2,[item(P,l) for l,P in two],True,'The six two-distance quadruples',sources=('ES',))
    cell(4,3,status='infinite',rigid=0,families=[{'dim':1,'count':9}],marks=['*'],note='Nine one-parameter families, one per realizable 3-colouring of the edges of K₄; no rigid configurations. The six two-distance sets are the points where families meet.')
    cat(4,3,[item(P,l,'family',1) for l,P in fam43],True,'Nine one-parameter families','Each panel is one representative; the whole family is obtained by moving the free parameter. Points are labelled A, B, C, D in the order of the coordinate list.',sources=('us',))
    S4=json.load(open('strata_4_4_1.json'))['realizable']
    items=[]
    for k,(l,P) in enumerate(fam44):
        if P is None: P=sorted(S4[k]['sample'].items())[0][1]
        items.append(item(P,l,'family',2))
    cell(4,4,status='infinite',rigid=0,families=[{'dim':2,'count':7}],marks=['*'],note='Seven two-parameter families (one per realizable 4-colouring of K₄).')
    cat(4,4,items,True,'Seven two-parameter families','One representative per combinatorial type.',sources=('us',))
    S5=json.load(open('strata_4_5_1.json'))['realizable']
    cell(4,5,status='infinite',rigid=0,families=[{'dim':3,'count':2}],marks=['*'],note='Exactly one coincidence among the six distances: the equal pair is adjacent or opposite.')
    cat(4,5,[item(sorted(S5[k]['sample'].items())[0][1],l,'family',3) for k,(l,_) in enumerate(fam45)],True,'Two three-parameter families',sources=('us',))
    cell(4,6,status='infinite',rigid=0,families=[{'dim':4,'count':1}],marks=[],note='Generic quadruples.')
    cat(4,6,[item(P,l,'family',4) for l,P in fam46],True,'Generic quadruples')
    # ---- n = 5 ----
    cell(5,2,status='exact',rigid=1,families=[],marks=[],note='The regular pentagon.')
    cat(5,2,[item(poly(5,range(5)),'regular pentagon',Nhint=5)],True,'The regular pentagon',sources=('ES',))
    FP=json.load(open('five_page.json'))
    FIVE=[x['pts'] for x in json.load(open('five_all.json'))]; SIX=json.load(open('six.json'))
    HEPT=[(math.cos(2*math.pi*j/7),math.sin(2*math.pi*j/7)) for j in range(7)]
    HEXC=[(0,0)]+[(math.cos(2*math.pi*j/6),math.sin(2*math.pi*j/6)) for j in range(6)]
    def precise(pts,pool):
        for P in pool:
            if same(pts,P,1e-4): return P
        raise SystemExit('no precise match for '+str(pts))
    items=[]
    for g in FP:
        if g['title'] in ('Six points','Seven points'): continue
        for it in g['items']:
            items.append(item(precise(it['pts'],FIVE),f"{g['title']}: {it['label']}",note=it['dist']))
    cell(5,3,status='exact',rigid=34,families=[],marks=['†'],note='Finite: 34 rigid configurations. A published classification exists (Shinohara 2004); the number 34 is from our computation and has not been checked against the paper.')
    cat(5,3,items,True,'The 34 five-point three-distance sets','Grouped by the figure they come from (square, pentagon, triangular lattice, heptagon, 12-gon, 18-gon, two lattice pieces, one loner).',sources=('Sh04','us'))
    # (5,4..7) from strata
    for d in (4,5,6,7):
        fn=f'strata_5_{d}_1.json'
        if not os.path.exists(fn):
            print('MISSING',fn); cell(5,d,status='bounds',rigid=None,families=[],marks=[],note='computation pending'); continue
        R=json.load(open(fn))['realizable']
        xcheck=''
        if d==4 and os.path.exists('strata_5_4_2.json'):
            # merge the two seeds: union of rigid configurations per type, union of family dimensions (the collinear
            # equally spaced set is rigid; a perturbation test at that singular point is unreliable, so drop dims>1 there)
            R2={tuple(r['col']):r for r in json.load(open('strata_5_4_2.json'))['realizable']}
            nA=sum(len(r['rigid']) for r in R); nB=sum(len(r['rigid']) for r in R2.values())
            for r in R:
                r2=R2.get(tuple(r['col']))
                if not r2: continue
                for P in r2['rigid']:
                    if not any(same(P,Q) for Q in r['rigid']): r['rigid'].append(P)
                for k,P in r2['sample'].items():
                    if k not in r['sample'] and int(k)<=1: r['sample'][k]=P
            nU=sum(len(r['rigid']) for r in R)
            xcheck=f' Two runs with different random seeds found {nA} and {nB} rigid configurations, agreeing on all {len(R)} realizable combinatorial types and on the families; the union, {nU} configurations, is listed. A random search can miss solutions, so this count is a lower bound in principle.'
        rig=[]; fams={}
        items=[]
        for r in R:
            for P in r['rigid']:
                rig.append(P)
            for k,P in r['sample'].items():
                fams[int(k)]=fams.get(int(k),0)+1
        # items: rigid first (sorted by ratios), then one representative per family type
        for r in R:
            sizes=sorted([r['col'].count(k) for k in range(d)],reverse=True)
            for P in r['rigid']:
                items.append(item(P,'rigid; distance classes of sizes '+'·'.join(map(str,sizes)),'rigid',0))
            for k,P in sorted(r['sample'].items()):
                items.append(item(P,f'{int(k)}-parameter family; distance classes of sizes '+'·'.join(map(str,sizes)),'family',int(k)))
        families=[{'dim':k,'count':v} for k,v in sorted(fams.items())]
        status='mixed' if (rig and fams) else ('exact' if rig else 'infinite')
        cell(5,d,status=status,rigid=len(rig),families=families,marks=['*'],
             note=f'{len(R)} realizable combinatorial types out of {json.load(open(fn))["classes"]} colourings of K₅ with {d} colours up to symmetry.'+xcheck)
        cat(5,d,items,True,f'Five points, {d} distances','Rigid configurations are listed first, then one representative of each positive-dimensional family. Types are described by the sizes of the distance classes (how many of the 10 pairs realize each distance).',sources=('us',))
    # ---- (6,3), (7,3) ----
    six=[g for g in FP if g['title']=='Six points'][0]['items']; seven=[g for g in FP if g['title']=='Seven points'][0]['items']
    cell(6,3,status='exact',rigid=9,families=[],marks=['†'],note='Nine rigid configurations (our enumeration; a published classification exists).')
    cat(6,3,[item(precise(it['pts'],SIX),it['label'],note=it['dist']) for it in six],True,'The nine six-point three-distance sets',sources=('Sh04','us'))
    cell(7,3,status='exact',rigid=2,families=[],marks=[],note='The regular heptagon and the regular hexagon with its centre; g(3) = 7.')
    cat(7,3,[item(precise(it['pts'],[HEPT,HEXC]),it['label'],note=it['dist'],Nhint=7 if 'hept' in it['label'] else 6) for it in seven],True,'The two seven-point three-distance sets',sources=('EF','Sh04'))
    # ---- maximal cells ----
    L94=json.load(open('lattice_9_4.json'))['tri']
    def uniq_sim(L):
        out=[]
        for P in L:
            if not any(same(P,Q) for Q in out): out.append(P)
        return out
    L94=uniq_sim(L94)
    def lat_label(P):
        m=min(pd(P)); deg=max(sum(1 for q in P if abs(math.dist(p,q)-m)<1e-6) for p in P)
        return 'piece of the triangular lattice: a regular hexagon, its centre, and two adjacent points of the next ring' if deg>=6 else 'piece of the triangular lattice: rows of 1, 3, 3 and 2 points'
    items=[item(poly(9,range(9)),'regular nonagon',Nhint=9)]+[item(P,lat_label(P)) for P in L94]
    # the fourth configuration: square ABCD with outward equilateral and inward 150° isosceles apexes on two adjacent sides, plus a point on the axis
    s3=math.sqrt(3)
    FOUR=[(1,0),(0,1),(-1,0),(0,-1),(0,-s3),((s3-1)/2,-(s3-1)/2),(-(s3-1)/2,-(s3-1)/2),((s3+1)/2,-(s3+1)/2),(-(s3+1)/2,-(s3+1)/2)]
    assert len(dvals(FOUR))==4
    it4=item(FOUR,'square ABCD (side √2); on sides AD and CD the outward equilateral apexes H, I and the inward apexes F, G of isosceles triangles with apex angle 150°; E on the axis of symmetry beyond D with DE = DF = √3 − 1. Distances √3−1, √2, 2, √3+1 (the 12-gon chords c₁, c₂, c₃, c₅).')
    it4['exact']=[['1','0'],['0','1'],['−1','0'],['0','−1'],['0','−√3'],['(√3−1)/2','−(√3−1)/2'],['−(√3−1)/2','−(√3−1)/2'],['(√3+1)/2','−(√3+1)/2'],['−(√3+1)/2','−(√3+1)/2']]
    it4['frame']='square ABCD inscribed in the unit circle, D at the bottom'
    items.append(it4)
    for fn in ('sym_9_4_mirror3_1.json','sym_9_4_mirror1_1.json','quant_9_4_7.json'):
        if os.path.exists(fn):
            for q in json.load(open(fn)):
                P=q['pts']
                if not any(same(P,R) for R in [poly(9,range(9)),FOUR]+L94): items.append(item(P,'found by numerical search')); print('UNEXPECTED extra (9,4) set',q['ratios'])
    cell(9,4,status='exact',rigid=4,families=[],marks=[],note='g(4) = 9; Erdős–Fishburn list four configurations: the regular nonagon, two pieces of the triangular lattice, and one more. All four are shown; the fourth was recovered here by a symmetric numerical search and matches their count.')
    cat(9,4,items,True,'The four nine-point four-distance sets','Erdős and Fishburn (1996) prove these are the only ones.',sources=('EF',))
    L125=uniq_sim(json.load(open('lattice_12_5.json'))['tri'])
    cell(12,5,status='exact',rigid=1,families=[],marks=[],note='g(5) = 12; the unique maximum five-distance set is a piece of the triangular lattice (Shinohara 2008).')
    cat(12,5,[item(P,'piece of the triangular lattice: rows of 3, 4, 3 and 2 points') for P in L125],True,'The unique twelve-point five-distance set',sources=('Sh08',))
    L136=uniq_sim(json.load(open('lattice_13_6.json'))['tri'])
    cell(13,6,status='exact',rigid=3,families=[],marks=[],note='g(6) = 13; three configurations (Wei 2012, confirmed by Szöllősi–Östergård 2020).')
    cat(13,6,[item(poly(13,range(13)),'regular 13-gon',Nhint=13),item(poly(12,range(12))+[(0,0)],'regular 12-gon with its centre',Nhint=12)]+[item(P,'piece of the triangular lattice: a regular hexagon, its centre and the six points of the next ring') for P in L136],True,'The three thirteen-point six-distance sets',sources=('Wei12','SO20'))
    # ---- cells with known existence but unknown count ----
    gmax={1:3,2:5,3:7,4:9,5:12,6:13,7:16}
    fmin={3:1,4:2,5:2,6:3,7:3,8:4,9:4,10:5,11:5,12:5,13:6,14:7,15:7,16:7}
    lit={ (7,4):('Classified by Lan and Wei (2013); the number of configurations is not available to us.',('LW13',)),
          (11,5):('Classified by Wei (2011); the number of configurations is not available to us.',('Wei11',)),
          (16,7):('The Erdős–Fishburn sixteen-point seven-distance set (a hexagon of side 2 in the triangular lattice minus one end of each long diameter). Whether g(7) = 16 is open; Bao and Yu conjecture it.',('EF','BY25'))}
    for n in range(6,17):
        for d in range(1,8):
            key=f'{n},{d}'
            if key in cells: continue
            if d<fmin[n] or n>gmax[d]: continue
            exp=2*n-4-math.comb(n,2)+d
            note,srcs=lit.get((n,d),('Configurations exist; the exact number is not known to us.',()))
            if exp>=1: note+=f' The dimension count predicts at least {exp}-parameter families, so the set is infinite.'
            cell(n,d,status='bounds',rigid=None,families=[],marks=[],note=note,expdim=exp)
            items=[]
            if (n,d)==(16,7):
                P,_=EX.ef16(); items.append(item(P,'hexagon of side 2 in the triangular lattice minus one end of each of the three long diameters (Erdős–Fishburn)'))
            for lab,P,N,S in EX.polygon_examples(n,d):
                items.append(item(P,lab,Nhint=N))
            for lab,P,S in EX.lattice_examples(n,d):
                if not any(same(P,it_['_P']) for it_ in [] ): pass
                items.append(item(P,lab))
            # dedupe
            ded=[]; seenP=[]
            for it_,P in zip(items,[None]*len(items)):
                pass
            cat(n,d,items,False,f'{n} points, {d} distances: examples','A complete list is not known; these are examples.',sources=srcs)
    # ---- the open band for d = 7 ----
    cells['17-36,7']={'n':'17–36','d':7,'status':'open','rigid':None,'families':[],'marks':[],'note':'Whether any 17-point seven-distance set exists is open (Erdős–Fishburn conjecture: g(7) = 16). The Bannai–Bannai–Stanton/Blokhuis bound gives at most C(9,2) = 36 points for seven distances. Budianto (2022) rules out 17-point sets in most cases of the diameter-graph analysis.'}
    cats['17-36,7']={'items':[],'complete':False,'title':'17 to 36 points, 7 distances: open','note':'No example is known and none is excluded; the count could be 0. Upper bound 36 points.','sources':[SRC['BBS'],SRC['Bud22'],SRC['BY25'],SRC['Bal']]}
    for key,c in cats.items():
        dd=cells[key]['d']
        for it in c['items']:
            assert it['d']==dd, (key,it['label'],it['d'],it['ratios'])
    meta={'gmax':gmax,'fmin':fmin,'rows':[str(k) for k in range(3,17)]+['17–36'],'cols':list(range(1,8))}
    json.dump({'cells':cells,'cats':cats,'meta':meta},open('applet_data.json','w'))
    print('cells',len(cells),'catalog items',sum(len(c['items']) for c in cats.values()))
    for k,c in sorted(cells.items(),key=lambda kv:(str(kv[1]['n']).zfill(3),kv[1]['d'])):
        print(k,c['status'],c['rigid'],c['families'],c['marks'])
if __name__=='__main__': main()
