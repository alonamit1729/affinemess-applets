"""Edge-coloured complete graphs ("patterns") up to vertex and colour permutations.

A pattern on n points with k colours is stored as a tuple col of length C(n,2) with entries in
range(k), indexed by EDGES[n] = list(itertools.combinations(range(n),2)), and normalised so that
colours appear in order 0,1,2,... along that edge order (restricted-growth form, 'rgf').

Canonical forms use nauty (pynauty) on the incidence encoding
    points 0..n-1 | edge-vertices n..n+m-1 | colour-vertices n+m..n+m+k-1
with three cells, so vertex permutations AND colour permutations are automorphisms of the encoding.
"""
import itertools
import pynauty

EDGES = {n: list(itertools.combinations(range(n), 2)) for n in range(2, 17)}
EIDX = {n: {e: i for i, e in enumerate(EDGES[n])} for n in EDGES}
MAXCOL = 6


def rgf(col):
    """Relabel colours by order of first appearance."""
    mp = {}
    out = []
    for c in col:
        if c not in mp:
            mp[c] = len(mp)
        out.append(mp[c])
    return tuple(out)


def ncolours(col):
    return max(col) + 1


def _graph(n, col):
    m = len(col)
    k = max(col) + 1
    adj = {}
    for e, (i, j) in enumerate(EDGES[n]):
        adj[n + e] = [i, j, n + m + col[e]]
    return pynauty.Graph(n + m + k, directed=False, adjacency_dict=adj,
                         vertex_coloring=[set(range(n)), set(range(n, n + m)), set(range(n + m, n + m + k))])


def certificate(n, col):
    """Bytes invariant under vertex and colour permutations (nauty certificate of the encoding)."""
    return pynauty.certificate(_graph(n, col))


def canon(n, col):
    """Canonical representative (rgf tuple) of the pattern col on n points."""
    lab = pynauty.canon_label(_graph(n, col))
    # positions 0..n-1 of the canonical labelling hold the point vertices
    pts = lab[:n]
    eidx = EIDX[n]
    out = []
    for i, j in EDGES[n]:
        a, b = pts[i], pts[j]
        out.append(col[eidx[(a, b) if a < b else (b, a)]])
    return rgf(out)


def autgroup_order(n, col):
    return pynauty.autgrp(_graph(n, col))[1]


def induced(n, col, verts):
    """Colouring of the complete graph induced on the sorted vertex list verts (rgf-normalised)."""
    eidx = EIDX[n]
    return rgf([col[eidx[(a, b)]] for a, b in itertools.combinations(verts, 2)])


def col_to_dict(n, col):
    return {f"{i},{j}": col[e] for e, (i, j) in enumerate(EDGES[n])}


def dict_to_col(n, d):
    return rgf([d[f"{i},{j}"] for i, j in EDGES[n]])


def all_rgf(m, maxk=MAXCOL):
    """All restricted-growth colourings of m edges with at most maxk colours."""
    def rec(i, cur, mx):
        if i == m:
            yield tuple(cur)
            return
        for b in range(min(mx + 1, maxk)):
            cur.append(b)
            yield from rec(i + 1, cur, max(mx, b + 1))
            cur.pop()
    yield from rec(0, [], 0)
