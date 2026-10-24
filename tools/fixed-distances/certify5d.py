"""Exact classification of 5-point d-distance sets for d = 5, 6, 7 (expected dimension >= 1).
Per colouring type: Singular minAssGTZ of the squared-distance equations (scale fixed by s0 = 1);
each minimal prime: dimension, degeneracy test (does an inequation polynomial vanish on it), and for
non-degenerate primes a real-point test by random hyperplane slices (dim hyperplanes -> zero-dim -> msolve
real roots with the inequations). Zero-dimensional primes: real roots directly.  Usage: python3 certify5d.py d"""
import itertools, subprocess, json, re, sys, math, time, random
from fractions import Fraction
import sympy as sp
MS='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
SG='/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/Singular'
n=5; d=int(sys.argv[1])
E=list(itertools.combinations(range(n),2)); m=len(E)
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
XV=['x1','x2','y2','x3','y3','x4','y4']; SV=[f's{k}' for k in range(1,d)]; V=XV+SV
SYM={v:sp.Symbol(v) for v in V+['t']}
coords={0:('0','0'),1:('x1','0'),2:('x2','y2'),3:('x3','y3'),4:('x4','y4')}
S=['1']+SV
def sq(i,j):
    (a,b),(c,e)=coords[i],coords[j]; return f'(({a})-({c}))^2+(({b})-({e}))^2'
facs=[f'({S[a]}-{S[b]})' for a,b in itertools.combinations(range(d),2)]+SV
ineq='t*'+'*'.join(facs)+'-1'
def expand(expr): return str(sp.expand(sp.sympify(expr,locals=SYM))).replace('**','^').replace(' ','')
def system(col): return [expand(f'{sq(i,j)}-({S[col[k]]})') for k,(i,j) in enumerate(E)]
def run_msolve(eqs,vars_):
    open('ms_in5.txt','w').write(','.join(vars_)+'\n0\n'+',\n'.join(eqs)+'\n')
    subprocess.run([MS,'-f','ms_in5.txt','-o','ms_out5.txt','-p','200'],capture_output=True,text=True,timeout=900)
    return open('ms_out5.txt').read()
def parse_real(out):
    s=out.replace('\n','').replace(' ','')
    if s.startswith('[-1]'): return 'empty',[]
    if s.startswith('[1,'): return 'posdim',[]
    tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s.rstrip(':'))
    def build(i):
        if tokens[i]=='[':
            lst=[]; i+=1
            while tokens[i]!=']':
                v,i=build(i); lst.append(v)
                if tokens[i]==',': i+=1
            return lst,i+1
        tk=tokens[i]
        if '^' in tk:
            num,den=tk.split('/'); base,ex=den.split('^'); return Fraction(int(num),int(base)**int(ex)),i+1
        return Fraction(tk),i+1
    tree,_=build(0)
    sols=tree[1][1] if isinstance(tree[1],list) and len(tree[1])>1 else []
    return 'zerodim',[[(lo+hi)/2 for lo,hi in sol] for sol in sols]
def decompose(eqs):
    ring='ring r=0,('+','.join(V)+'),dp;'
    for algo,tmo in (('minAssGTZ',60),('minAssChar',120),('minAssGTZ',900)):
        script=ring+'\nLIB "primdec.lib";\nideal i='+',\n'.join(eqs)+';\nlist L='+algo+'(i);\nint k; ideal G; int bad;\nfor(k=1;k<=size(L);k++){ G=std(L[k]); print("PRIME"); print(dim(G)); bad=0; '+' '.join(f'if(reduce({f},G)==0){{bad=1;}}' for f in facs)+' print(bad); print(L[k]); }\nquit;\n'
        open('sg_in5.sing','w').write(script)
        try:
            r=subprocess.run([SG,'-q','sg_in5.sing'],capture_output=True,text=True,timeout=tmo); break
        except subprocess.TimeoutExpired:
            r=None; print('   decomposition timeout with',algo,tmo,flush=True)
    if r is None: return None
    primes=[]; cur=None
    for line in r.stdout.splitlines():
        line=line.strip()
        if line=='PRIME': cur={'dim':None,'degenerate':None,'gens':[]}; primes.append(cur); continue
        if cur is None or line.startswith('//') or not line: continue
        if cur['dim'] is None and re.fullmatch(r'-?\d+',line): cur['dim']=int(line); continue
        if cur['degenerate'] is None and re.fullmatch(r'[01]',line): cur['degenerate']=bool(int(line)); continue
        cur['gens'].append(line.rstrip(','))
    return primes
def points_from(mid):
    x1,x2,y2,x3,y3,x4,y4=[float(v) for v in mid[:7]]; return [(0.0,0.0),(x1,0.0),(x2,y2),(x3,y3),(x4,y4)]
def dm(P): return [[math.dist(P[i],P[j]) for j in range(n)] for i in range(n)]
def congruent(P,Q):
    D1=dm(P); D2=dm(Q)
    for p in itertools.permutations(range(n)):
        if all(abs(D1[i][j]-D2[p[i]][p[j]])<1e-9 for i in range(n) for j in range(i+1,n)): return True
    return False
import os
results=json.load(open(f'cert5{d}.json')) if os.path.exists(f'cert5{d}.json') else []
done_idx={r['idx'] for r in results}
t0=time.time()
for ci,col in enumerate(classes):
    if ci in done_idx: continue
    eqs=system(col); primes=decompose(eqs)
    rec={'idx':ci,'col':list(col),'primes':[],'rigid':[]}
    if primes is None:
        rec['unresolved']=True
        o=run_msolve(eqs+[expand(ineq)],V+['t']); st,_=parse_real(o)
        rec['msolve_status']=st
        results.append(rec); print(ci,'/',len(classes),'UNRESOLVED',st,'%.0fs'%(time.time()-t0),flush=True)
        json.dump(results,open(f'cert5{d}.json','w')); continue
    random.seed(ci)
    for p in primes:
        q={'dim':p['dim'],'degenerate':p['degenerate'],'real':None}
        if not p['degenerate']:
            if p['dim']==0:
                st,md=parse_real(run_msolve(p['gens']+[expand(ineq)],V+['t']))
                for mid in md:
                    P=points_from(mid)
                    if not any(congruent(P,Q) for Q in rec['rigid']): rec['rigid'].append(P)
                q['real']=len(md)>0
            else:
                hits=0
                for _ in range(4):
                    lins=['+'.join(f'({random.randint(-5,5)})*{v}' for v in XV)+f'-({random.randint(-3,3)})' for _ in range(p['dim'])]
                    st,md=parse_real(run_msolve(p['gens']+[expand(l) for l in lins],V))
                    good=[]
                    for mid in (md if st=='zerodim' else []):
                        vals=[1.0]+[float(v) for v in mid[7:7+d-1]]
                        if min(abs(a-b) for a,b in itertools.combinations(vals,2))>1e-9 and min(vals)>1e-9: good.append(mid)
                    if good:
                        hits+=1
                        if 'point' not in q: q['point']=points_from(good[0])
                q['real']=hits>0; q['slice_hits']=hits
        rec['primes'].append(q)
    results.append(rec)
    dims=sorted(set(q['dim'] for q in rec['primes'] if q['real']))
    print(ci,'/',len(classes),'real dims',dims,'rigid',len(rec['rigid']),'%.0fs'%(time.time()-t0),flush=True)
    json.dump(results,open(f'cert5{d}.json','w'))
print('done')
