import json, os, collections, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from symmetry import e_dim
lines = []
summary = {}
for d in (4, 5, 6, 7):
    fn = os.path.join(HERE, f'samples_5{d}.json')
    if not os.path.exists(fn): continue
    S = json.load(open(fn))
    cnt = collections.Counter(); bydim = collections.defaultdict(collections.Counter); mixed = []; nos = []
    for k, r in S.items():
        types = collections.Counter(s['type'] for s in r['samples'])
        if not types: verdict = 'no samples'; nos.append(k)
        elif len(types) == 1: verdict = list(types)[0]
        elif 'trivial' in types: verdict = 'mixed'; mixed.append((k, dict(types)))
        else: verdict = 'mixed-sym'; mixed.append((k, dict(types)))
        r['verdict'] = verdict
        cnt[verdict] += 1; bydim[r['dim']][verdict] += 1
    summary[d] = {'types': len(S), 'verdicts': dict(cnt), 'bydim': {k: dict(v) for k, v in bydim.items()}, 'mixed': mixed, 'nosamples': nos}
    lines.append(f"(5,{d}) e={e_dim(5,d)}: {len(S)} positive-dimensional types; by dim: " + '; '.join(f"dim {k}: " + ', '.join(f'{t}:{c}' for t, c in sorted(v.items(), key=lambda kv: -kv[1])) for k, v in sorted(bydim.items())))
    if mixed: lines.append('   mixed types: ' + str(mixed))
    if nos: lines.append('   no samples: ' + str(nos))
    json.dump(S, open(fn, 'w'))
open(os.path.join(HERE, 'step2_summary.json'), 'w').write(json.dumps(summary, indent=1))
print('\n'.join(lines))
