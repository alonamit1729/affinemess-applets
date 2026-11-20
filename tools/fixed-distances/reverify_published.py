"""Re-verify the published cells of the table by exact computation and compare with the catalog.

(6,3), (7,3), (7,4), (8,4), (9,4): every pattern with exactly d colours in the certified candidate list
L_n of the (16, 7) proof (../16-7-proof/certificates/L{n}.json) is solved with fewdist.solve_lengths
(SymPy Groebner basis of the Gram minors with the Rabinowitsch equation, real roots in exact number
fields); a flexible pattern would raise ValueError.  The real Euclidean realisations are counted up to
similarity (fewdist.shape_key).

(11,5), (12,5), (13,6), (11,6): the sets found by msolve real root isolation on the complete candidate
lists L11-L13 (g7count/gram_n.py -> g7count/sets1113_by_d.json) are compared with the catalog.

Each list is matched to the catalog of its cell one to one by an invariant (the sorted list of the
points' sorted squared-distance profiles, scaled to minimum 1, rounded to 4 decimals); the script checks
that the invariant separates the sets within each list.  Writes data/reverify_published.json.

    cd tools/fixed-distances && python3 reverify_published.py        # about 5 seconds
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
PROOF = os.path.join(HERE, '..', '16-7-proof')
sys.path.insert(0, PROOF); sys.path.insert(0, HERE)
import fewdist as fd
import datautil

SOLVE = [(6, 3, 'Shinohara 2004, Theorem 2'), (7, 3, 'Erdős–Fishburn 1996; Shinohara 2004'),
         (7, 4, 'Lan–Wei 2013, Theorem 8'), (8, 4, 'Shinohara 2008, Theorem 1.2'), (9, 4, 'Erdős–Fishburn 1996')]
GRAM = [(11, 5, 'Wei 2011, Theorem 13'), (12, 5, 'Erdős–Fishburn 1996; Shinohara 2008'),
        (13, 6, 'Wei 2012; Szöllősi–Östergård 2020'), (11, 6, 'Amit–Astra 2026, Lemma 6')]


def profile(P, nd=4):
    n = len(P)
    D = [[(P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2 for j in range(n)] for i in range(n)]
    m = min(D[i][j] for i in range(n) for j in range(n) if i != j)
    return tuple(sorted(tuple(sorted(round(D[i][j] / m, nd) for j in range(n) if j != i)) for i in range(n)))


def n_distances(P, nd=6):
    n = len(P)
    D = [(P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2 for i in range(n) for j in range(i + 1, n)]
    m = min(D)
    return len({round(x / m, nd) for x in D})


def compare(computed, catalog):
    a = [profile(P) for P in computed]
    b = [profile(P) for P in catalog]
    return {'computed': len(a), 'catalog': len(b),
            'computed_distinct': len(set(a)) == len(a), 'catalog_distinct': len(set(b)) == len(b),
            'match': sorted(a) == sorted(b)}


def main():
    data = datautil.load(os.path.join(HERE, 'data', 'applet_data.json'))
    catalog = lambda n, d: [[[float(x), float(y)] for x, y in it['coords']] for it in data['cats'][f'{n},{d}']['items']]
    out = {}
    for n, d, src in SOLVE:
        t = time.time()
        L = json.load(open(os.path.join(PROOF, 'certificates', f'L{n}.json')))['patterns']
        pats = [tuple(p) for p in L if max(p) + 1 == d]
        sets, flexible = [], []
        for i, p in enumerate(pats):
            try:
                sols = fd.solve_lengths(p)
            except ValueError:
                flexible.append(i)
                continue
            shapes = {}
            for s in sols:
                for r in s.embeddings:
                    shapes.setdefault(fd.shape_key(s, r), s.realisation(r).xy())
            sets += [[list(map(float, q)) for q in P] for P in shapes.values()]
        assert all(n_distances(P) == d for P in sets)
        rec = {'published': src, 'method': 'certified L%d, %d patterns with %d colours, fewdist.solve_lengths' % (n, len(pats), d),
               'flexible_patterns': flexible, 'seconds': round(time.time() - t, 1)}
        rec.update(compare(sets, catalog(n, d)))
        out[f'{n},{d}'] = rec
        print(n, d, rec, flush=True)
    G = json.load(open(os.path.join(HERE, 'g7count', 'sets1113_by_d.json')))
    for n, d, src in GRAM:
        sets = [x['pts'] for x in G[f'{n},{d}']]
        rec = {'published': src, 'method': 'msolve on the complete candidate list L%d (g7count/gram_n.py)' % n}
        rec.update(compare(sets, catalog(n, d)))
        out[f'{n},{d}'] = rec
        print(n, d, rec, flush=True)
    ok = all(r['match'] and r['computed_distinct'] and r['catalog_distinct'] and not r.get('flexible_patterns')
             for r in out.values())
    json.dump({'all_match': ok, 'cells': out}, open(os.path.join(HERE, 'data', 'reverify_published.json'), 'w'),
              ensure_ascii=False, indent=1)
    print('all match' if ok else 'MISMATCH')


if __name__ == '__main__':
    main()
