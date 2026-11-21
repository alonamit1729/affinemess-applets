import json,collections
rows=json.load(open('final_rows.json')); CH={(r['cell'],r['i']):r for r in json.load(open('chain_out.json'))}
seen={(c,i) for c,i,*_ in rows}
for o in json.load(open('rigid.json')):
  if (o['cell'],o['i']) not in seen: rows.append((o['cell'],o['i'],'1 ruler+compass',[],None,o['label']))
def ch(c,i):
  r=CH[(c,i)]; h=r['how']
  if h is None: return 'Z'
  if h in('lattice-tri','lattice-sq'): return 'A'
  if h=='polygon-subset': return 'B%d'%r['m']
  if h.startswith('lattice') or h=='pair+chain': return 'C'
  return 'D%d'%r['m']
X=collections.defaultdict(collections.Counter)
for c,i,v,ks,_,lab in rows: X[v][ch(c,i)]+=1
for v in sorted(X): print(v, dict(sorted(X[v].items())))
Y=collections.defaultdict(collections.Counter)
for c,i,v,ks,_,lab in rows:
  k=ch(c,i); Y[c][('easy' if k!='Z' else 'HARD')+'/'+v[0]]+=1
key=lambda c:tuple(map(int,c.split(',')))
for c in sorted(Y,key=key): print(c, dict(sorted(Y[c].items())))
json.dump([[c,i,v,ch(c,i),lab] for c,i,v,ks,_,lab in rows],open('cross.json','w'))
