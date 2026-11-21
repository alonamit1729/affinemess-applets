"""python3 classify.py <probe dir> <recipes dir> [--write]
Turn galois_probe.py's results into verdicts (the rules of analysis/final.py) and, with --write, fill them into
../verdicts.json where it has null. Per squared-distance ratio:
  every Frobenius cycle order a power of 2            -> straightedge and compass
  abelian (all cycle types of equal length)            -> needs a regular polygon (cyclotomic)
  cubic: three real roots -> trisection; one -> a cube root
  degree 4-6: the Galois group (sympy): S5, S6, A5, A6, PGL2(5), PSL2(5) -> not solvable;
              S4 -> trisection (0 or 4 real roots) or 's4'; S4xC2 -> trisection (6 real roots)
  degree 7 with cycle orders 10 or 12                  -> not solvable (S7/A7)
  anything else with a cycle order that is not a power of 2 -> 'nrc': beyond straightedge and compass (which other
                                                          tool, if any, would do is not decided)
  no minimal polynomial found                          -> undecided ('?')
A set's verdict is the worst over its ratios. A polygon verdict names the m of the polygon its construction uses
when that construction uses one straightedge and compass cannot draw (p<m>), else just 'p'."""
import json, os, sys, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__))
GAUSS = {3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34}
_groups = {}

def group(poly):
    key = tuple(poly)
    if key not in _groups:
        try:
            from sympy import Poly, symbols
            from sympy.polys.numberfields.galoisgroups import galois_group
            x = symbols('x')
            G, _ = galois_group(Poly(list(reversed(poly)), x), by_name=True, max_tries=60)
            _groups[key] = str(G).split('.')[-1]
        except Exception:
            _groups[key] = None
    return _groups[key]

def pow2(k): return k & (k - 1) == 0

def rkind(q, res):
    if res > 1e-200: return 'rc'
    if q.get('deg') is None: return '?'
    if all(pow2(o) for o in q['orders']): return 'rc'
    if q['equal_deg']: return 'p'
    if q['deg'] == 3: return 'tri' if q['real_roots'] == 3 else 'cbrt'
    if 4 <= q['deg'] <= 6 and q.get('poly'):
        name = group(q['poly'])
        if name in ('S5', 'S6', 'A5', 'A6', 'PGL2F5', 'PSL2F5'): return name if name in ('S5', 'S6') else 'S5'
        if name == 'S4': return 'tri' if q['real_roots'] in (0, 4) else 's4'
        if name == 'S4xC2' and q['real_roots'] == 6: return 'tri'
        if name: return 'nrc'
    if q['deg'] == 7 and any(o in q['orders'] for o in (10, 12)): return 'S7'
    return 'nrc'

RANK = {'rc': 0, 'p': 1, 'tri': 2, 's4': 3, 'nrc': 4, '?': 5, 'cbrt': 6, 'S5': 7, 'S6': 7, 'S7': 7}

def verdict(r, poly_m, level=None):
    v = verdict0(r, poly_m)
    # a set built with a regular polygon that straightedge and compass cannot draw is constructible with it; when the
    # analysis rules out straightedge and compass, the polygon is what it needs (square-root extensions of a cyclotomic
    # field are not abelian, so the cycle-type test alone calls them 'nrc')
    if level == 3 and poly_m and poly_m not in GAUSS and v != 'rc': return 'p%d' % poly_m
    return v

def verdict0(r, poly_m):
    if 'ratios' not in r: return '?'
    ks = [rkind(q, r.get('res', 0)) for q in r['ratios']]
    if not ks: return 'rc'
    worst = max(ks, key=lambda k: RANK[k])
    if worst == 'p': return 'p%d' % poly_m if poly_m and poly_m not in GAUSS else 'p'
    return worst

if __name__ == '__main__':
    probe_dir, rec_dir = sys.argv[1:3]
    V = {}
    for f in sorted(glob.glob(os.path.join(probe_dir, '*.json'))):
        cell = os.path.basename(f)[:-5]
        polys, levels = {}, {}
        for g in glob.glob(os.path.join(rec_dir, cell, '[0-9]*.json')):
            for num, it in json.load(open(g))['items'].items():
                m = [a[3] for a in it[6]['n'] if a[0] == 'poly']
                polys[int(num)] = m[0] if m else 0; levels[int(num)] = it[2]
        for num, r in json.load(open(f)).items():
            V[(cell.replace('_', ','), int(num))] = verdict(r, polys.get(int(num), 0), levels.get(int(num)))
    tally = collections.Counter(V.values())
    print(len(V), 'sets:', dict(tally.most_common()))
    if '--write' in sys.argv:
        fn = os.path.join(HERE, '..', 'verdicts.json'); D = json.load(open(fn)); k = 0
        for row in D:
            if row[2] is None and (row[0], row[1]) in V: row[2] = V[(row[0], row[1])]; k += 1
        open(fn, 'w').write(json.dumps(D, separators=(',', ':'))); print('filled', k, 'verdicts into verdicts.json')
