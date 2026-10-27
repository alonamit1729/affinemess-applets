"""numerical trace along every family with e<=0 and n>=6 (and the (5,4) ones): symmetry types seen at generic
points of the trace vs at the catalog sample."""
import json, os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from famcheck import check_item
rows = json.load(open(os.path.join(HERE, 'step1_families.json')))
out = []
for r in rows:
    if r['e'] > 0: continue
    info = check_item(r['coords'], r['cls'])
    seen = collections.Counter()
    for p in info.get('trace', []):
        for s in p['syms'][1:]:
            if s not in ('boundary', 'lost', 'singular'): seen[s] += 1
    rec = {'cell': r['cell'], 'idx': r['idx'], 'e': r['e'], 'dim_catalog': r['dim'], 'tangent_dim': info['tangent_dim'], 'sample_sym': info['sym0'],
           'trace_syms': dict(seen), 'res': info['res_proj'], 'label': r['label']}
    out.append(rec)
    flag = '' if set(seen) <= {info['sym0']} else '   <-- symmetry changes along the trace'
    print(f"{r['cell']} #{r['idx']} e={r['e']} dim={r['dim']} tdim={info['tangent_dim']} sample={info['sym0']} trace={dict(seen)}{flag}", flush=True)
json.dump(out, open(os.path.join(HERE, 'step3_trace.json'), 'w'), indent=1)
