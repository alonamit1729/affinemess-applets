"""The 11-point planar sets with at most six distances (Amit-Astra 2026, Lemma 6): one-point deletions of the
eleven 12-point lattice cores (g7/cores.json) with <= 6 distances, and the circular ones (subsets of R13,
R12 ∪ O, R11 ∪ O, R10 ∪ O).  Coordinates to 100 digits.  Also checks their patterns against g7fin/patterns_11.json."""
import json, os, sys, itertools
import mpmath as mp
from extend import mp, d2, dist_values, Classes, pattern, normalise, TOL
SCR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S3 = mp.sqrt(3)


def tri(a, b):
    return (mp.mpf(a) + mp.mpf(b) / 2, mp.mpf(b) * S3 / 2)


def lattice_sets():
    cores = json.load(open(os.path.join(SCR, 'g7', 'cores.json')))
    C = Classes()
    for ci, L in enumerate(cores):
        for i in range(12):
            A = [tri(*p) for j, p in enumerate(L) if j != i]
            if len(dist_values(A)) <= 6:
                C.add(A, ('lattice', ci, i))
    return [(r[0], r[1]) for r in C.reps]


def regular(n, R=1):
    return [(R * mp.cos(2 * mp.pi * k / n), R * mp.sin(2 * mp.pi * k / n)) for k in range(n)]


def circular_sets():
    out = []
    R13 = regular(13)
    for k in range(1, 7):
        out.append(([p for i, p in enumerate(R13) if i not in (0, k)], ('R13-{0,%d}' % k,)))
    R12 = regular(12); O = (mp.mpf(0), mp.mpf(0))
    out.append((R12[1:], ('R12-v',)))
    for k in range(1, 7):
        out.append(([p for i, p in enumerate(R12) if i not in (0, k)] + [O], ('R12+O-{0,%d}' % k,)))
    out.append((regular(10) + [O], ('R10+O',)))
    out.append((regular(11), ('R11',)))
    out.append((regular(11)[1:] + [O], ('R11-v+O',)))
    C = Classes()
    for P, info in out:
        assert len(P) == 11 and len(dist_values(P)) <= 6, info
        C.add(P, info)
    return [(r[0], r[1]) for r in C.reps]


def all_sets():
    return lattice_sets() + circular_sets()


if __name__ == '__main__':
    sys.path.insert(0, os.path.join(SCR, 'g7fin', 'gen'))
    import patlib
    L = lattice_sets(); Cc = circular_sets()
    print('lattice 11-point sets:', len(L), ' circular:', len(Cc))
    pats = json.load(open(os.path.join(SCR, 'g7fin', 'patterns_11.json')))
    canon = {}
    for p in pats:
        col = [p['col'][f'{i},{j}'] for i, j in itertools.combinations(range(11), 2)]
        canon.setdefault(patlib.canon(11, patlib.rgf(col)), []).append(p['name'])
    print('distinct canonical patterns in patterns_11.json:', len(canon))
    seen = set()
    for P, info in L + Cc:
        c = patlib.canon(11, pattern(P))
        assert c in canon, info
        seen.add(c)
    print('patterns realised by these', len(L) + len(Cc), 'sets:', len(seen), 'of', len(canon))
