"""Degeneration poset from exact (degen_{d}.json) and sampled (degen_s_{d}.json) data.  Usage: python3 degen_analyze2.py d
Edges: (type T of (5,d)) -> (type T' of (5,d'), d'<d) when T' occurs on the boundary of T (one coincidence equation);
the stratum T' is labelled rigid/family from cert5{d'}.json.  Exact data lists components with dimensions; sampled data
lists types met by random slices.  Output degen_poset_{d}.json and a summary."""
import json, sys, itertools, collections, os, math
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
status={}
for dd in range(2,d):
    try: R=json.load(open(f'cert5{dd}.json'))
    except FileNotFoundError: continue
    for r in R:
        key=canon(tuple(r['col']))
        fam=any((q.get('real') or q.get('real_slices')) and q['dim']>0 for q in r['primes']) or r.get('unresolved')
        if r.get('status')=='empty' or (not r['rigid'] and not fam): st='unrealizable'
        elif r['rigid'] and not fam: st='rigid'
        else: st='family'
        status[(dd,key)]=(st,r['idx'],len(r['rigid']))
SAT=json.load(open(f'sat5{d}.json'))
exact=collections.defaultdict(set); exact_types=set(); timeouts=set()
if os.path.exists(f'degen_{d}.json'):
    for rec in json.load(open(f'degen_{d}.json')):
        if rec.get('timeout'): timeouts.add((rec['type'],tuple(rec['eq']))); continue
        exact_types.add(rec['type'])
        for comp in rec['components']:
            if not comp['real']: continue
            if comp.get('coincident'): exact[rec['type']].add(('degenerate',comp['dim'])); continue
            key=canon(tuple(comp['col_type'])); st=status.get((comp['d_new'],key),('?',None,None))
            exact[rec['type']].add(((comp['d_new'],st[1],st[0]),comp['dim']))
sampled=collections.defaultdict(set)
if os.path.exists(f'degen_s_{d}.json'):
    for rec in json.load(open(f'degen_s_{d}.json')):
        for s_ in rec['samples']:
            if s_['coincident']: sampled[rec['type']].add('degenerate'); continue
            # re-classify from the coordinates with a tight tolerance (the samples are 200-bit roots)
            P=s_['sample']; ds=[math.dist(P[i],P[j]) for i,j in E]; mx=max(ds); vals=[]; col=[]
            for x in ds:
                for k,v in enumerate(vals):
                    if abs(x-v)<1e-10*mx: col.append(k); break
                else: vals.append(x); col.append(len(vals)-1)
            dn=len(vals)
            if dn>=d: continue
            key=canon(tuple(col)); st=status.get((dn,key),('?',None,None))
            if st[0]=='unrealizable': print('   warning: sample of type',rec['type'],'classified as unrealizable type',st[1]); continue
            sampled[rec['type']].add((dn,st[1],st[0]))
types=[int(k) for k,v in SAT.items() if not v.get('timeout')]
# compare where both exist (types fully exact: no timeouts)
full_exact=[t for t in types if t in exact_types and not any(tt==t for tt,_ in timeouts)]
miss=0; tot=0; extra=0
for t in full_exact:
    ex={x for x,k in exact[t] if x!='degenerate'}; sa={x for x in sampled[t] if x!='degenerate'}
    tot+=len(ex); miss+=len(ex-sa); extra+=len(sa-ex)
print(f'(5,{d}): {len(types)} positive-dimensional types; exact data complete for {len(full_exact)}; sampled data for {len(sampled)}')
if full_exact: print(f'  sampling check on the exact types: {tot} exact incidences, {miss} missed by sampling, {extra} found only by sampling')
# merged poset: exact where complete, else sampled
poset={}
for t in types:
    if t in full_exact: poset[t]={'src':'exact','to':sorted({x for x,k in exact[t]},key=str),'dims':sorted({(str(x),k) for x,k in exact[t]})}
    else: poset[t]={'src':'sampled','to':sorted(sampled[t],key=str)}
summary=collections.Counter()
for t,v in poset.items():
    dimT=SAT[str(t)]['dim']
    for x in v['to']:
        if x=='degenerate': summary[('degenerate',dimT)]+=1
        else: summary[(f'to (5,{x[0]}) {x[2]}',dimT)]+=1
print('  incidences by (kind, dim of family):')
for k,v in sorted(summary.items(),key=str): print('    ',k,v)
incoming=collections.Counter(x for t,v in poset.items() for x in v['to'] if x!='degenerate')
for dn in range(2,d):
    rig=sum(1 for x in incoming if x[0]==dn and x[2]=='rigid'); fam=sum(1 for x in incoming if x[0]==dn and x[2]=='family')
    tot_r=sum(1 for k,v in status.items() if k[0]==dn and v[0]=='rigid'); tot_f=sum(1 for k,v in status.items() if k[0]==dn and v[0]=='family')
    print(f'  (5,{dn}): rigid types reached {rig} of {tot_r}; family types reached {fam} of {tot_f}')
unreached={dn:[v[1] for k,v in status.items() if k[0]==dn and v[0]!='unrealizable' and not any(x[0]==dn and x[1]==v[1] for x in incoming)] for dn in range(2,d)}
print('  unreached realizable lower types:',{k:v for k,v in unreached.items() if v})
json.dump({str(t):{'src':v['src'],'to':[list(x) if isinstance(x,tuple) else x for x in v['to']]} for t,v in poset.items()},open(f'degen_poset_{d}.json','w'))
