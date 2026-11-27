"""Planar sets of eleven or more points with at most six distances: all of them, up to similarity.

Every such set contains ten points with at most six distances (any ten of its points), and those
are similar to one of the 282 sets of Proposition 6 of the paper (data/solutions.json). So every
such set is a core K plus some points, with at most 6 - |Delta(K)| new squared distances in all.
fewdist.grow lists every way to add points within that budget, in exact arithmetic (paper,
Section 6; the budget is 0 or 1, so at least eight points of K are at old distances from a new
point unless it is a circumcentre, and no core has more than five collinear points).

Result (paper, Theorem 28): 4, 64, 1, 14 and 3 sets with (points, distances) = (11, 5), (11, 6),
(12, 5), (12, 6), (13, 6), and none with fourteen points; hence g(5) = 12 and g(6) = 13.

    python3 six_distances.py          (under half a minute)
"""
import os
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "fifteen"))
import fewdist as fd                                            # noqa: E402
from fifteen import shape_key                                   # noqa: E402

EXPECTED = {(11, 5): 4, (11, 6): 64, (12, 5): 1, (12, 6): 14, (13, 6): 3}


def main():
    t0 = time.time()
    patterns = fd.load_patterns(os.path.join(ROOT, "data", "L10.json"))
    solutions = fd.load_solutions(os.path.join(ROOT, "data", "solutions.json"), patterns)
    shapes = defaultdict(set)          # (points, distances) -> similarity classes
    largest = Counter()                # largest size reached -> number of cores
    for sols in solutions:
        for sol in sols:
            K = fd.place(sol.pattern, sol.lengths)       # exact, valid at every real root
            cands, subsets = fd.grow(K, 6 - len(K.distances()))
            largest[10 + max(map(len, subsets))] += len(sol.embeddings)
            for s in subsets:
                if not s:
                    continue
                X = K.with_points(cands[j] for j in s)
                cell = (len(X.points), len(X.distances()))
                for root in sol.embeddings:
                    X.root = root
                    shapes[cell].add(shape_key(X, root, sol.field))
    print(f"all 282 cores extended in {time.time() - t0:.0f} s")
    print("largest number of points reached -> number of cores:", dict(sorted(largest.items())))
    for cell in sorted(shapes):
        print(f"{cell[0]} points, {cell[1]} distances: {len(shapes[cell])} sets up to similarity")
    assert {c: len(v) for c, v in shapes.items()} == EXPECTED
    assert sum(largest.values()) == 282 and max(largest) == 13
    print("as stated in the paper")


if __name__ == "__main__":
    main()
