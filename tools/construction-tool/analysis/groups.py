import json, sys
from sympy import Poly, symbols
from sympy.polys.numberfields.galoisgroups import galois_group
x = symbols('x')
res = []
for f in sys.argv[2:]:
    res += json.load(open(f))
polys = {}
for r in res:
    for q in r.get('ratios', []):
        if q.get('poly') and 4 <= q['deg'] <= 6:
            polys[tuple(q['poly'])] = None
print(len(polys), 'polys of degree 4..6')
out = {}
for co in polys:
    P = Poly(list(reversed(co)), x)
    try:
        G, alt = galois_group(P, by_name=True, max_tries=60)
        out[str(list(co))] = [str(G), G.order() if hasattr(G, 'order') else None]
    except Exception as e:
        out[str(list(co))] = ['ERR ' + repr(e)[:80], None]
json.dump(out, open(sys.argv[1], 'w'), indent=0)
import collections
print(collections.Counter(v[0] for v in out.values()))
