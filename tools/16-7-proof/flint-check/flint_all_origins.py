"""FLINT (independent of SymPy's Groebner engine) check of completeness of the saved length solutions.
Generators: the Cayley-Menger planarity condition for EVERY 4 of the 10 points (3x3 Gram minors over all base points),
plus the Rabinowitsch polynomial w*q - 1 forcing all lengths non-zero and distinct.
Test:  reduced GB of <planarity, w*q-1>  ==  reduced GB of <saved basis B, w*q-1>."""
import sys, time, json, ast, resource, multiprocessing as mp
from itertools import combinations
import sympy as sp, flint
AUD=sys.argv[1]
P=json.load(open(f'{AUD}/L10.json'))['patterns']
pairs=list(combinations(range(10),2))
saved={}
for line in open(f'{AUD}/palettes.log'):
    i,done,k,gs,_=ast.literal_eval(line); saved[i]=gs
saved[12]="[x0-3/13,x1-7/13,x2-1/13,x3-9/13,x4-4/13]"; saved[142]="[x0-7/12,x1-1/12,x2-3/4,x3-1/4,x4-1/3]"
def setup(p):
    m=max(p)+1; xs=sp.symbols('x:%d'%(m-1)); w=sp.Symbol('w'); C=[sp.Integer(1),*xs]
    D=[[0]*10 for _ in range(10)]
    for (a,b),c in zip(pairs,p): D[a][b]=D[b][a]=C[c]
    gens=set()
    for quad in combinations(range(10),4):
        o,t=quad[0],quad[1:]
        M=sp.Matrix(3,3,lambda r,c: D[o][t[r]]+D[o][t[c]]-D[t[r]][t[c]])
        g=sp.expand(M.det())
        if g!=0: gens.add(g)
    q=sp.Integer(1)
    for x in xs: q*=x
    for u,v in combinations(C,2): q*=(u-v)
    return xs,w,sorted(gens,key=str),sp.expand(w*q-1)
def gb(polys,syms):
    ctx=flint.fmpz_mpoly_ctx.get(tuple(str(s) for s in syms),'degrevlex')
    out=[]
    for f in polys:
        pf=sp.Poly(f,*syms); den=sp.ilcm(1,*[c.q for c in pf.coeffs()])
        out.append(ctx.from_dict({tuple(int(e) for e in mon): int(c*den) for mon,c in pf.terms()}))
    G=flint.fmpz_mpoly_vec(out,ctx).buchberger_naive().autoreduction(groebner=True)
    assert G.is_groebner()
    return G
def nsol(G,nv):
    lead=[tuple(int(e) for e in g.monoms()[0]) for g in G]
    std={tuple([0]*nv)}; fr=list(std)
    while fr:
        new=[]
        for mon in fr:
            for k in range(nv):
                mm=list(mon); mm[k]+=1; mm=tuple(mm)
                if mm in std or any(all(a>=b for a,b in zip(mm,l)) for l in lead): continue
                std.add(mm); new.append(mm)
        fr=new; assert len(std)<5000
    return len(std)
def one(i):
    xs,w,gens,sat=setup(P[i]); syms=(w,*xs)
    G=gb([*gens,sat],syms)
    B=sp.sympify(saved[i],locals={str(x):x for x in xs})
    G2=gb([*B,sat],syms)
    return i, sorted(map(str,G))==sorted(map(str,G2)), nsol(G,len(syms)), len(gens)
def run(i,q):
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    try: q.put(one(i))
    except BaseException as e: q.put((i,'ERR',repr(e)[:80],0))
if __name__=='__main__':
    ids=[int(a) for a in sys.argv[2:]] or range(len(P)); res={}; t0=time.time()
    for i in ids:
        q=mp.Queue(); pr=mp.Process(target=run,args=(i,q)); pr.start()
        try: r=q.get(timeout=45)
        except Exception: r=(i,'TIMEOUT',None,0)
        pr.kill(); pr.join(); res[i]=r
        print(r,flush=True)
    ok=[i for i,r in res.items() if r[1] is True]
    from collections import Counter
    print("patterns with FLINT-confirmed equality:",len(ok),"of",len(res))
    print("not confirmed:",[(i,r[1]) for i,r in res.items() if r[1] is not True])
    print("non-degenerate solutions, total:",sum(res[i][2] for i in ok)," per-pattern histogram:",dict(sorted(Counter(res[i][2] for i in ok).items())))
    print("planarity equations per pattern: min %d max %d"%(min(r[3] for r in res.values() if r[3]),max(r[3] for r in res.values())),' time %.0fs'%(time.time()-t0))
