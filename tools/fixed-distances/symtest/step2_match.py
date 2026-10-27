"""match catalog family items of (5,d) to the sampled types by the canonical form of the edge colouring; compare
the catalog sample's symmetry with the generic symmetry."""
import json, os, sys, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
E = list(itertools.combinations(range(5), 2)); PERMS = list(itertools.permutations(range(5)))
def canon_col(col):
    best = None
    for p in PERMS:
        c = [None] * 10; relabel = {}
        for k, (i, j) in enumerate(E):
            a, b = sorted((p[i], p[j])); c[E.index((a, b))] = col[k]
        # relabel colours by first appearance
        out = []
        for x in c:
            if x not in relabel: relabel[x] = len(relabel)
            out.append(relabel[x])
        t = tuple(out)
        if best is None or t < best: best = t
    return best
rows = json.load(open(os.path.join(HERE, 'step1_families.json')))
for d in (4, 5, 6, 7):
    S = json.load(open(os.path.join(HERE, f'samples_5{d}.json')))
    bycan = collections.defaultdict(list)
    for k, r in S.items(): bycan[canon_col(r['col'])].append(k)
    cmp = collections.Counter(); unmatched = 0; special = []
    for r in rows:
        if r['cell'] != f'5,{d}': continue
        ks = bycan.get(canon_col(r['cls']))
        if not ks: unmatched += 1; continue
        gen = {S[k]['verdict'] for k in ks}
        cmp[(r['type'], tuple(sorted(gen)))] += 1
        if r['type'] not in gen: special.append((r['idx'], r['type'], sorted(gen), ks))
    print(f'(5,{d}): catalog items matched; (catalog sample symmetry, generic symmetry of type) counts:')
    for (a, b), c in sorted(cmp.items(), key=lambda kv: -kv[1]): print(f'    {a:8s} generic={b}: {c}')
    print('   unmatched', unmatched, '; catalog samples more symmetric than generic:', special)
