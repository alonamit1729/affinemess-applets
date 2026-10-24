"""Layer 1: exact classification of 5-point 4-distance sets, type by type.
For each of the 371 colouring types of K5 with 4 colours (up to S5):
  variables x1 ; x2,y2 ; x3,y3 ; x4,y4 ; s1,s2,s3 (squared distance values; s0 = 1 fixes the scale) ; t (Rabinowitsch)
  equations: |Pi-Pj|^2 - s_{class} = 0 (10) ; t * prod(s_i - s_j) * s1*s2*s3 - 1 = 0  (distinct, nonzero values)
  msolve decides: no complex solution / zero-dimensional (then isolates the real solutions) / positive-dimensional.
  Positive-dimensional types are decomposed with Singular (minAssGTZ over Q); each minimal prime is handled by dimension.
Output: cert54.json"""
import itertools, subprocess, json, re, sys, os, math, time
from fractions import Fraction
MS='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
SG='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/Singular'
n=5; d=4
E=list(itertools.combinations(range(n),2)); m=len(E)
# ---- types (same enumeration as strata.py) ----
def partitions(items,k):
    def rec(i,cur,mx):
        if i==len(items):
            if mx==k: yield tuple(cur)
            return
        for b in range(min(mx+1,k)):
            cur.append(b); yield from rec(i+1,cur,max(mx,b+1)); cur.pop()
    yield from rec(0,[],0)
perms=[[E.index(tuple(sorted((p[a],p[b])))) for a,b in E] for p in itertools.permutations(range(n))]
def canon(rgs):
    best=None
    for pe in perms:
        c=[rgs[pe[i]] for i in range(m)]; mp={}; out=[]
        for x in c:
            if x not in mp: mp[x]=len(mp)
            out.append(mp[x])
        t=tuple(out)
        if best is None or t<best: best=t
    return best
classes=sorted(set(canon(r) for r in partitions(E,d)))
assert len(classes)==371, len(classes)
VARS=['x1','x2','y2','x3','y3','x4','y4','s1','s2','s3','t']
coords={0:('0','0'),1:('x1','0'),2:('x2','y2'),3:('x3','y3'),4:('x4','y4')}
def sq(i,j):
    (a,b),(c,e)=coords[i],coords[j]
    return f'(({a})-({c}))^2+(({b})-({e}))^2'
S=['1','s1','s2','s3']
ineq='t*'+'*'.join(f'({S[a]}-{S[b]})' for a,b in itertools.combinations(range(4),2))+'*s1*s2*s3-1'
import sympy as sp
SYM={v:sp.Symbol(v) for v in VARS}
def expand(expr):
    return str(sp.expand(sp.sympify(expr,locals=SYM))).replace('**','^').replace(' ','')
def system(col):
    eqs=[f'{sq(i,j)}-({S[col[k]]})' for k,(i,j) in enumerate(E)]
    return [expand(e) for e in eqs+[ineq]]
def run_msolve(eqs,extra_vars=(),timeout=600):
    vars_=VARS+list(extra_vars)
    fn='ms_in.txt'; open(fn,'w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    r=subprocess.run([MS,'-f',fn,'-o','ms_out.txt','-p','200','-t','2'],capture_output=True,text=True,timeout=timeout)
    out=open('ms_out.txt').read()
    return out
def parse_real(out):
    """msolve real-solution output: [0, [1, [ [ [a,b], [c,d], ... ], ... ] ]]  -> list of solutions (midpoints as Fractions)"""
    s=out.replace('\n','').replace(' ','')
    if s.startswith('[-1]'): return 'empty',[]
    if s.startswith('[1,'): return 'posdim',[]
    # zero-dimensional: find the solution list: structure [0,[nbsols? ...]] ; parse all innermost [a,b] pairs grouped by solution
    # msolve format: [0, [1, [ [ [x1lo,x1hi], [x2lo,x2hi], ... ], [ ... ] ] ] ]
    try:
        import ast
        s2=s.rstrip(':')
        # convert rationals a/b to tuples for ast: replace "p/q" by "Fraction('p/q')" is unsafe for ast; parse manually
        tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s2)
        def build(i):
            if tokens[i]=='[':
                lst=[]; i+=1
                while tokens[i]!=']':
                    v,i=build(i)
                    lst.append(v)
                    if tokens[i]==',': i+=1
                return lst,i+1
            tk=tokens[i]
            if '^' in tk:
                num,den=tk.split('/'); base,ex=den.split('^'); return Fraction(int(num),int(base)**int(ex)),i+1
            return Fraction(tk),i+1
        tree,_=build(0)
        sols=tree[1][1] if isinstance(tree[1],list) and len(tree[1])>1 else []
        mids=[]
        for sol in sols:
            mids.append([ (lo+hi)/2 for lo,hi in sol ])
        return 'zerodim',mids
    except Exception as e:
        return 'parse_error:'+str(e)[:80]+' :: '+s[:200],[]
def points_from(mid):
    x1,x2,y2,x3,y3,x4,y4=[float(v) for v in mid[:7]]
    return [(0.0,0.0),(x1,0.0),(x2,y2),(x3,y3),(x4,y4)]
def dm(P): return [[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(n)):
        if all(abs(D1[i][j]-D2[p[i]][p[j]])<1e-9 for i in range(n) for j in range(i+1,n)): return True
    return False
def singular_decompose(eqs):
    """minimal associated primes over Q of the distance equations alone (no Rabinowitsch variable);
    for each prime report its dimension and whether some inequation polynomial vanishes identically on it (degenerate)"""
    V=[v for v in VARS if v!='t']
    ring='ring r=0,('+','.join(V)+'),dp;'
    eqs=[e for e in eqs if not re.search(r'\bt\b',e)]
    facs=[f'({S[a]}-{S[b]})' for a,b in itertools.combinations(range(4),2)]+['s1','s2','s3']
    script=ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(eqs)+';\nlist L=minAssGTZ(i);\nint k; int j; ideal G; int bad;\nfor(k=1;k<=size(L);k++){ G=std(L[k]); print("PRIME"); print(dim(G)); bad=0; '+' '.join(f'if(reduce({f},G)==0){{bad=1;}}' for f in facs)+' print(bad); print(L[k]); }\nquit;\n'
    open('sg_in.sing','w').write(script)
    r=subprocess.run([SG,'-q','sg_in.sing'],capture_output=True,text=True,timeout=1800)
    out=r.stdout
    primes=[]; cur=None
    for line in out.splitlines():
        line=line.strip()
        if line=='PRIME': cur={'dim':None,'degenerate':None,'gens':[]}; primes.append(cur); continue
        if cur is None: continue
        if cur['dim'] is None and re.fullmatch(r'-?\d+',line): cur['dim']=int(line); continue
        if cur['degenerate'] is None and re.fullmatch(r'[01]',line): cur['degenerate']=bool(int(line)); continue
        if line.startswith('//') or not line: continue
        cur['gens'].append(line.rstrip(','))
    return primes,out
results=[]
t0=time.time()
only=[int(a) for a in sys.argv[1:]]
for ci,col in enumerate(classes):
    if only and ci not in only: continue
    eqs=system(col)
    out=run_msolve(eqs)
    status,mids=parse_real(out)
    rec={'idx':ci,'col':list(col),'status':status,'rigid':[],'primes':[]}
    if status=='zerodim':
        confs=[]
        for mid in mids:
            P=points_from(mid)
            if not any(congruent(P,Q) for Q in confs): confs.append(P)
        rec['rigid']=confs; rec['n_real_solutions']=len(mids)
    elif status=='posdim':
        primes,raw=singular_decompose(eqs)
        rec['primes']=[{'dim':p['dim'],'degenerate':p['degenerate'],'gens':p['gens']} for p in primes]
        # zero-dimensional non-degenerate primes: isolate their real points with msolve (inequations re-imposed via t)
        confs=[]
        for p in primes:
            if p['dim']==0 and not p['degenerate']:
                o=run_msolve(p['gens']+[expand(ineq)])
                st,md=parse_real(o)
                for mid in md:
                    P=points_from(mid)
                    if not any(congruent(P,Q) for Q in confs): confs.append(P)
        rec['rigid']=confs
        # one-dimensional non-degenerate primes: probe for real points with random hyperplane slices
        import random
        random.seed(ci)
        for p in rec['primes']:
            if p['dim']==1 and not p['degenerate']:
                hits=0
                for _ in range(4):
                    lin='+'.join(f'({random.randint(-5,5)})*{v}' for v in ['x1','x2','y2','x3','y3','x4','y4'])+f'-({random.randint(-3,3)})'
                    o=run_msolve(p['gens']+[expand(ineq),expand(lin)])
                    st,md=parse_real(o)
                    if st=='zerodim' and md: hits+=1
                p['real_slices']=hits
    results.append(rec)
    print(ci,status,'rigid',len(rec['rigid']),'primes',[p['dim'] for p in rec['primes']],'%.0fs'%(time.time()-t0),flush=True)
    json.dump(results,open('cert54.json' if not only else 'cert54_part.json','w'))
print('done')
