"""Recompute the minimal primes of the types that have a non-degenerate component without real points found by slices,
store the generators, and print the short generators for inspection. Usage: python3 noreal.py d"""
import sys, json
d=int(sys.argv[1])
src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
C=json.load(open(f'cert5{d}.json'))
out=[]
for r in C:
    bad=[pi for pi,p in enumerate(r['primes']) if not p['degenerate'] and not p['real']]
    if not bad: continue
    primes=decompose(system(classes[r['idx']]))
    assert primes is not None and len(primes)==len(r['primes']), (r['idx'],len(primes),len(r['primes']))
    for pi in bad:
        p=primes[pi]; assert p['dim']==r['primes'][pi]['dim'] and not p['degenerate']
        out.append({'d':d,'idx':r['idx'],'prime':pi,'dim':p['dim'],'gens':p['gens']})
        print(f'd={d} type {r["idx"]} prime {pi} dim {p["dim"]} ngens {len(p["gens"])}; short gens:',[g for g in p['gens'] if len(g)<70][:10],flush=True)
json.dump(out,open(f'noreal_{d}.json','w'))
