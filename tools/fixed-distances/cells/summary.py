"""Summarize cells/cert_{n}_{d}.json: rigid configurations up to similarity (across patterns), families (patterns with a
positive-dimensional non-degenerate real component) by dimension, hidden-rigid candidates, complex-only patterns.
Usage: python3 cells/summary.py n d [--write]  (writes cells/rigid_{n}_{d}.json, cells/fam_{n}_{d}.json)"""
import json, sys, os, collections, math, itertools
HERE=os.path.dirname(os.path.abspath(__file__))
n=int(sys.argv[1]); d=int(sys.argv[2])
sys.argv=[sys.argv[0],str(n),str(d)]
src=open(os.path.join(HERE,'cellcert.py')).read()
exec(src[src.index('def congruent'):src.index('def singular')])
R=json.load(open(os.path.join(HERE,f'cert_{n}_{d}.json')))
rig=[]; fams=collections.Counter(); famlist=[]; hidden=[]; st=collections.Counter(); noreal=[]
buckets={}
def fpkey(P):
    ds=sorted(math.dist(a,b) for a,b in itertools.combinations(P,2)); m=ds[0]
    return tuple(round(x/m,4) for x in ds)
for key,r in sorted(R.items(),key=lambda kv:int(kv[0])):
    st[r['status']]+=1
    for P,rank in r['rigid']:
        b=buckets.setdefault(fpkey(P),[])
        if not any(congruent(P,rig[i][0]) for i in b):
            b.append(len(rig)); rig.append((P,rank,r['idx']))
    dims=set()
    for q in r['primes']:
        if q['degenerate'] or q['dim']==0: continue
        if q.get('real_psd'): dims.add(q['dim']); 
        if not q.get('real_psd'): noreal.append((r['idx'],q['dim']))
        if q.get('singular_real') not in (0,None): hidden.append((r['idx'],q['dim'],q['singular_real']))
    for k in dims:
        fams[k]+=1
        q=[q for q in r['primes'] if not q['degenerate'] and q['dim']==k and q.get('real_psd')][0]
        famlist.append({'idx':r['idx'],'col':r['col'],'dim':k,'point':q['point'],'ranks':q.get('ranks'),'components':sum(1 for q in r['primes'] if not q['degenerate'] and q['dim']==k and q.get('real_psd'))})
print(f'({n},{d}): patterns {len(R)} status {dict(st)}; rigid {len(rig)} (collinear {sum(1 for P,rk,_ in rig if rk==1)}); families by dim {dict(fams)}; posdim primes without real PSD point {len(noreal)} {noreal[:10]}; singular-locus real PSD points {hidden[:10]}')
if '--write' in sys.argv or True:
    json.dump([{'pts':P,'rank':rk,'pat':i} for P,rk,i in rig],open(os.path.join(HERE,f'rigid_{n}_{d}.json'),'w'))
    json.dump(famlist,open(os.path.join(HERE,f'fam_{n}_{d}.json'),'w'))
