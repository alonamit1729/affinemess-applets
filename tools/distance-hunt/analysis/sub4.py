import json, itertools, fields, chain
R = json.load(open('rigid.json')); X = json.load(open('cross.json'))
V = {(c, i): v for c, i, v, k, l in X}
imp = [o for o in R if V[(o['cell'], o['i'])].startswith('5')]
print(len(imp), 'impossible')
fields.DPS = 400
def pow2(k): return k & (k - 1) == 0
res = []
for idx, o in enumerate(R):
    if o not in imp: continue
    z = [complex(float(a), float(b)) for a, b in o['coords']]
    for drop in range(5):
        sub = [z[k] for k in range(5) if k != drop]
        how = chain.analyse({'coords': [[str(q.real), str(q.imag)] for q in sub]})
        if how['how'] is not None:
            res.append((idx, o['label'], V[(o['cell'], o['i'])], drop, how))
for r in res: print(r)
