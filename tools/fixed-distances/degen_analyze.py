"""Degeneration poset from degen_{d}.json (d = 4..7): for every positive-dimensional type of (5,d), the strata of
(5,d') (d' < d) in its closure, each identified by its canonical colouring and its dimension, with the (5,d') type's
own status (rigid / family) from cert5{d'}.json.  Usage: python3 degen_analyze.py d"""
import json, sys, itertools, collections
d=int(sys.argv[1])
E=list(itertools.combinations(range(5),2)); m=10
perms=[[E.index(tuple(sorted((p[a],p[b])))) for a,b in E] for p in itertools.permutations(range(5))]
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
# status of every type of every smaller d
status={}
for dd in range(2,d):
    try: R=json.load(open(f'cert5{dd}.json'))
    except FileNotFoundError: continue
    for r in R:
        key=canon(tuple(r['col']))
        if r.get('status')=='empty' or (not r['rigid'] and not any((q.get('real') or q.get('real_slices')) and q['dim']>0 for q in r['primes']) and not r.get('unresolved')):
            st='unrealizable'
        elif r['rigid'] and not any((q.get('real') or q.get('real_slices')) and q['dim']>0 for q in r['primes']): st='rigid'
        elif r.get('unresolved'): st='family'
        else: st='family'
        status[(dd,key)]=(st,r['idx'],len(r['rigid']))
D=json.load(open(f'degen_{d}.json'))
SAT=json.load(open(f'sat5{d}.json'))
bytype=collections.defaultdict(list)
for rec in D:
    for comp in rec['components']:
        if not comp['real']: continue
        if comp.get('coincident'): bytype[rec['type']].append(('degenerate',comp['dim'])); continue
        key=canon(tuple(comp['col_type'])); dn=comp['d_new']
        st=status.get((dn,key),('?',None,None))
        bytype[rec['type']].append(((dn,st[1],st[0]),comp['dim']))
summary=collections.Counter(); edges=collections.Counter(); per_type={}
for t,L in bytype.items():
    S=set(L); per_type[t]=sorted(S,key=str)
    dimT=SAT[str(t)]['dim']
    for x,k in S:
        if x=='degenerate': summary[('degenerate',dimT,k)]+=1
        else: summary[('to d=%d %s'%(x[0],x[2]),dimT,k)]+=1; edges[(t,x)]+=1
print(f'(5,{d}): {len(bytype)} positive-dimensional types with boundary data; timeouts {sum(1 for r in D if r.get("timeout"))}')
print('boundary strata by (kind, dim of family, dim of stratum):')
for k,v in sorted(summary.items(),key=str): print('  ',k,v)
# per lower stratum: how many families of (5,d) have it in their closure
incoming=collections.Counter(x for (t,x) in edges)
for dn in range(2,d):
    rig=[x for x in incoming if x[0]==dn and x[2]=='rigid']; fam=[x for x in incoming if x[0]==dn and x[2]=='family']
    tot_r=sum(1 for k,v in status.items() if k[0]==dn and v[0]=='rigid'); tot_f=sum(1 for k,v in status.items() if k[0]==dn and v[0]=='family')
    print(f'  (5,{dn}): rigid types reached {len(rig)} of {tot_r}; family types reached {len(fam)} of {tot_f}')
unreached={dn:[v[1] for k,v in status.items() if k[0]==dn and v[0]!='unrealizable' and not any(x[0]==dn and x[1]==v[1] for x in incoming)] for dn in range(2,d)}
print('  unreached realizable lower types:',{k:v for k,v in unreached.items() if v})
json.dump({'per_type':{str(t):[list(x) if isinstance(x,tuple) else x for x in v] for t,v in per_type.items()}},open(f'degen_poset_{d}.json','w'))
