import json, math, itertools, re, sys
sys.path.insert(0,'.')
from fractions import Fraction
import realize
from realize import system_coords, msolve_real
pats=[tuple(c) for c in json.load(open('new3/L10.json'))['patterns']]
E=list(itertools.combinations(range(10),2))
out=[]
for idx,col in enumerate(pats):
    V,P=system_coords(10,col)
    txt=msolve_real(V,P,precision=100,timeout=1800)
    s=txt.replace('\n','').replace(' ','').rstrip(':')
    tokens=re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]',s)
    def build(i):
        if tokens[i]=='[':
            lst=[];i+=1
            while tokens[i]!=']':
                v,i=build(i);lst.append(v)
                if tokens[i]==',':i+=1
            return lst,i+1
        tk=tokens[i]
        if '^' in tk:
            num,den=tk.split('/');base,ex=den.split('^');return Fraction(int(num),int(base)**int(ex)),i+1
        return Fraction(tk),i+1
    tree,_=build(0); sols=tree[1][1]
    print('pattern',idx,'real solutions',len(sols))
    recs=[]
    for sol in sols:
        vals=[float((a+b)/2) for a,b in sol]
        vv=dict(zip(V,vals))
        P10=[(0.0,0.0),(vv['x1'],0.0)]+[(vv[f'x{i}'],vv[f'y{i}']) for i in range(2,10)]
        ds=sorted(set(round(math.dist(P10[i],P10[j]),9) for i,j in E))
        m=ds[0]
        recs.append({'pts':P10,'dists':[d/m for d in ds]})
    # distinct shapes
    shapes={}
    for r in recs:
        key=tuple(round(x,6) for x in r['dists'])
        shapes.setdefault(key,r)
    print('  distinct distance spectra:',len(shapes))
    for key,r in shapes.items(): print('   ',[round(x,5) for x in key])
    out.append({'col':list(col),'nsol':len(sols),'shapes':list(shapes.values())})
json.dump(out,open('new3/new3_solutions.json','w'),indent=1)
