import math, random, itertools
def nd(P,tol=1e-7):
    ds=sorted(math.dist(P[i],P[j]) for i,j in itertools.combinations(range(4),2))
    k=1
    for i in range(1,6):
        if ds[i]-ds[i-1]>tol: k+=1
    return k, ds
S3=math.sqrt(3)
def A(t):  # equilateral 012, 3 on unit circle about 0
    return [(0,0),(1,0),(.5,S3/2),(math.cos(t),math.sin(t))]
def B(al): # rhombus angle al
    return [(0,0),(1,0),(math.cos(al),math.sin(al)),(1+math.cos(al),math.sin(al))]
def C(ph): # 0 center, 1 at 0, 2,3 at +-ph
    return [(0,0),(1,0),(math.cos(ph),math.sin(ph)),(math.cos(ph),-math.sin(ph))]
def D(t): # equilateral 012, 3 on bisector of 01
    return [(0,0),(1,0),(.5,S3/2),(.5,t)]
def E_trap(th):
    return [(0,0),(1,0),(math.cos(th),math.sin(th)),(1-math.cos(th),math.sin(th))]
def E_par(b):
    # parallelogram 0-2-1-3, 02=1,21=b,01=1: triangle 0,2,1 sides 1,b,1 apex 2
    x=b/2; y=math.sqrt(max(1-x*x,0))
    P2=(.5-0.0,0) # placeholder
    # 0=(0,0),1=(1,0); 2 at distance 1 from 0 and b from 1
    x2=(1+1-b*b)/2; y2=math.sqrt(max(1-x2*x2,0))
    return [(0,0),(1,0),(x2,y2),(1-x2,-y2)]
def F(h,s=1):
    # 1=(0,0),2=(1,0),0=(.5,h); 3 at dist 1 from 1 and 0
    P1=(0,0);P2=(1,0);P0=(.5,h)
    b=math.hypot(.5,h); mid=(.25,h/2); d=(-h,.5); n=math.hypot(*d)
    t=math.sqrt(max(1-b*b/4,0))
    P3=(mid[0]+s*t*d[0]/n, mid[1]+s*t*d[1]/n)
    return [P0,P1,P2,P3]
def G(s):
    return [(0,s),(-.5,0),(.5,0),(0,s+1)]
def H(ph,s=1):
    P2=(0,0);P1=(1,0);P3=(math.cos(2*ph),math.sin(2*ph)); b=2*math.sin(ph)
    dx=math.cos(2*ph)-.5
    r=b*b-dx*dx
    if r<0: return None
    y=math.sin(2*ph)+s*math.sqrt(r)
    return [(.5,y),P1,P2,P3]
def I(r):
    return [(0,0),(1,0),(0,r),(1,r)]
random.seed(1)
fams={'A':lambda:A(random.uniform(0,2*math.pi)),'B':lambda:B(random.uniform(0,math.pi/2)),
'C':lambda:C(random.uniform(0,math.pi)),'D':lambda:D(random.uniform(-3,3)),
'Etrap':lambda:E_trap(random.uniform(0,math.pi)),'Epar':lambda:E_par(random.uniform(0,2)),
'F+':lambda:F(random.uniform(0,1.7),1),'F-':lambda:F(random.uniform(0,1.7),-1),
'G':lambda:G(random.uniform(-.5,3)),'H+':lambda:H(random.uniform(0.3,math.pi/2),1),'H-':lambda:H(random.uniform(0.3,math.pi/2),-1),
'I':lambda:I(random.uniform(0,1))}
for k,f in fams.items():
    cnt={}
    for _ in range(2000):
        P=f()
        if P is None: continue
        n,_=nd(P); cnt[n]=cnt.get(n,0)+1
    print(k,cnt)
# H tangency
import scipy.optimize as so
g=lambda ph: 2*math.sin(ph)-abs(math.cos(2*ph)-.5)
phm=so.brentq(g,0.1,0.6); print('H tangency phi deg',math.degrees(phm), nd(H(phm,1)))
print('H phi=90', nd(H(math.pi/2-1e-9,1)))
# H collinear 0,2,3 members (heptagon?)
for s in (1,-1):
    def col(ph):
        P=H(ph,s); 
        if P is None: return 1
        (x0,y0),_,(x2,y2),(x3,y3)=P
        return x0*y3-y0*x3
    xs=[0.3+i*0.001 for i in range(1200)]
    vals=[col(x) for x in xs]
    for i in range(1,len(xs)):
        if vals[i-1]*vals[i]<0:
            r=so.brentq(col,xs[i-1],xs[i]); print('H collinear 0,2,3 at phi',math.degrees(r),'s',s,nd(H(r,s)))
# F endpoints
print('F h=0',nd(F(0,1)),nd(F(0,-1)))
# E junction
print('E collinear',nd(E_trap(math.pi-1e-9)),nd(E_par(2-1e-9)))
