"""The candidates at each level, generated a second time by an independent program.

A pattern on n points with at most six colours is a CANDIDATE if each of its sub-patterns on
n-1 points is in the list KEPT at the level below.  This program lists all candidates on n
points, for n = 5..10, and checks that every one of them is accounted for: it is either kept
(L5.json ... L10.json) or has a certificate (level5.jsonl, level6.jsonl.gz, level7.jsonl.gz).

    python3 generate.py N [--jobs K] [--work DIR] [--certificates FILE]

It shares no code with the generator used first (gen.py, built on nauty).  Pure Python: the
only import from outside the standard library is fewdist.canonical, the canonical form of the
notebook.

How it works.  Take each kept pattern Q on n-1 points and add a point v.  The n-1 new pairs
are coloured one at a time, each with one of the colours present so far or with one new colour
(at most six in all).  After the pair (j, v) is coloured, all sub-patterns on v, j and some of
the points 0..j-1 are complete, and each must itself be kept at its own level, because the kept
lists are closed under taking sub-patterns.  So a partial colouring is dropped as soon as one
of them is not kept; nothing else is ever dropped.  For sub-patterns on 4, 5 or 6 points this
test is a table lookup: the tables hold every colouring of the complete graph on 4, 5, 6
labelled points that is a kept pattern up to renumbering.  When n - 1 > 6, the n-1
sub-patterns on n-1 points that contain v are tested at the end, by canonical form.

Every candidate P arises this way: remove any point v of P; what is left is a kept pattern,
hence a renumbering of some Q in the list; renumber P accordingly, and rename the colours of
the new pairs that do not occur in Q in the order in which they first appear.
"""
import gzip
import json
import os
import re
import sys
import time
from array import array
from bisect import bisect_left
from itertools import combinations, permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from fewdist import canonical, n_points, pairs                          # noqa: E402

MAX_COLOURS = 6
TABLE_SIZES = (4, 5, 6)


def code(colours):
    """Colours renamed in order of first appearance, then packed into one integer."""
    names, value = {}, 0
    for c in colours:
        value = value * MAX_COLOURS + names.setdefault(c, len(names))
    return value


def labelled_table(patterns, m):
    """Codes of all colourings of the complete graph on m labelled points that are one of
    `patterns` after renumbering the points; a sorted array."""
    P = pairs(m)
    where = {e: i for i, e in enumerate(P)}
    renumberings = [[where[tuple(sorted((s[a], s[b])))] for a, b in P]
                    for s in permutations(range(m))]
    codes = set()
    for p in patterns:
        for r in renumberings:
            codes.add(code([p[i] for i in r]))
    return array("q", sorted(codes))


def in_table(table, value):
    i = bisect_left(table, value)
    return i < len(table) and table[i] == value


def extensions(Q, tables, kept_below=None):
    """All candidates obtained from the kept pattern Q on m points by adding the point m.
    tables[s]: labelled table of the kept patterns on s points (s = 4, 5, 6);
    kept_below: canonical forms of the kept patterns on m points (needed only if m > 6)."""
    m = n_points(Q)
    old = [[None] * m for _ in range(m)]                # old[a][b]: colour of the old pair
    for (a, b), c in zip(pairs(m), Q):
        old[a][b] = old[b][a] = c
    new = [None] * m                                    # new[a]: colour of the pair (a, m)
    top = min(m, max(TABLE_SIZES))                      # largest size tested by table lookup
    found = []

    def sub_code(points):                 # the sub-pattern on these old points and the new one
        names, value = {}, 0
        for i, a in enumerate(points):
            for b in points[i + 1:]:
                value = value * MAX_COLOURS + names.setdefault(old[a][b], len(names))
            value = value * MAX_COLOURS + names.setdefault(new[a], len(names))
        return value

    def whole():                                        # the pattern on all m + 1 points
        return tuple(new[a] if b == m else old[a][b] for a, b in pairs(m + 1))

    def without(d):                                     # the sub-pattern without old point d
        rest = [v for v in range(m) if v != d] + [m]
        return tuple(new[a] if b == m else old[a][b] for a, b in combinations(rest, 2))

    def colour_pair(j, used):             # colour the pair (j, m); `used` colours so far
        if j == m:
            if m == top or all(canonical(without(d)) in kept_below for d in range(m)):
                found.append(whole())
            return
        size = min(j + 2, top)            # the sub-patterns completed now have this many points
        completed = [T + (j,) for T in combinations(range(j), size - 2)]
        if size < 4:                      # no list below four points: nothing to test yet
            completed = []
        for c in range(min(used + 1, MAX_COLOURS)):
            new[j] = c
            if all(in_table(tables[size], sub_code(points)) for points in completed):
                colour_pair(j + 1, max(used, c + 1))
        new[j] = None

    colour_pair(0, max(Q) + 1)
    return found


# ------------------------------------------------------------------ running a level
def load(name):
    return [tuple(p) for p in json.load(open(os.path.join(HERE, name)))["patterns"]]


def kept(n):
    """The patterns kept on n points."""
    if n == 4:      # every pattern on four points except the one-coloured one
        strings = [()]
        for _ in range(6):
            strings = [s + (c,) for s in strings for c in range(max(s, default=-1) + 2)]
        return sorted({canonical(s) for s in strings} - {(0,) * 6})
    return load(f"L{n}.json")


def build_tables(n, work=None):
    """The labelled tables needed for the candidates on n points."""
    tables = {}
    for s in TABLE_SIZES:
        if s <= n - 1:
            cache = work and os.path.join(work, f"table{s}.bin")
            if cache and os.path.exists(cache):
                tables[s] = array("q")
                with open(cache, "rb") as fh:
                    tables[s].frombytes(fh.read())
            else:
                tables[s] = labelled_table(kept(s), s)
                if cache:
                    with open(cache, "wb") as fh:
                        tables[s].tofile(fh)
    return tables


_STATE = {}


def _work(chunk):
    out = set()
    for Q in chunk:
        for P in extensions(Q, _STATE["tables"], _STATE["below"]):
            out.add(canonical(P))
    return sorted(out)


def candidates(n, jobs=1, work=None, log=print):
    """Canonical forms of all candidates on n points."""
    parents = kept(n - 1)
    t0 = time.time()
    _STATE["tables"] = build_tables(n, work)
    _STATE["below"] = {canonical(p) for p in parents} if n - 1 > max(TABLE_SIZES) else None
    sizes = [(s, len(t)) for s, t in _STATE["tables"].items()]
    log(f"{n} points: {len(parents)} kept patterns below; tables {sizes}"
        f" ({time.time() - t0:.0f} s)")
    size = max(1, min(200, len(parents) // (8 * jobs) + 1))
    chunks = [parents[i:i + size] for i in range(0, len(parents), size)]
    name_of = lambda i: work and os.path.join(work, f"gen{n}_{i:04d}.json")
    result, todo = set(), []
    for i, chunk in enumerate(chunks):
        if name_of(i) and os.path.exists(name_of(i)):
            result.update(tuple(p) for p in json.load(open(name_of(i))))
        else:
            todo.append((i, chunk))
    pool = None
    if jobs > 1 and todo:                 # worker processes inherit the tables (fork)
        import multiprocessing
        if "fork" in multiprocessing.get_all_start_methods():
            pool = multiprocessing.get_context("fork").Pool(jobs)
    stream = (pool.imap if pool else map)(_work, [chunk for _, chunk in todo])
    for done, ((i, _), found) in enumerate(zip(todo, stream)):
        result.update(tuple(p) for p in found)
        if name_of(i):
            with open(name_of(i) + ".tmp", "w") as fh:
                json.dump(found, fh, separators=(",", ":"))
            os.replace(name_of(i) + ".tmp", name_of(i))
        if (done + 1) % max(1, len(todo) // 10) == 0:
            log(f"   {done + 1}/{len(todo)} chunks, {len(result)} candidates so far"
                f" ({time.time() - t0:.0f} s)")
    if pool:
        pool.close()
    return result


def certified(name):
    """The patterns of the certificates in a file, in canonical form, by number of points."""
    name = str(name)
    opener = gzip.open if name.endswith(".gz") else open
    find = re.compile(r'"n":(\d+),"pattern":\[([0-9,]*)\]')
    out = {}
    with opener(name, "rt") as fh:
        for line in fh:
            n, pattern = find.search(line[:300]).groups()
            out.setdefault(int(n), set()).add(canonical(tuple(map(int, pattern.split(",")))))
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    n = int(args[0])
    option = lambda flag, default: args[args.index(flag) + 1] if flag in args else default
    jobs, work = int(option("--jobs", 1)), option("--work", None)
    if work:
        os.makedirs(work, exist_ok=True)
    t0 = time.time()
    found = candidates(n, jobs, work)
    kept_here = {canonical(p) for p in kept(n)}
    print(f"{n} points: {len(found)} candidates; {len(kept_here)} kept"
          f" ({time.time() - t0:.0f} s)", flush=True)
    assert kept_here <= found, "a kept pattern is not a candidate"
    rest = found - kept_here
    if n >= 8:
        print("every candidate is kept:", not rest)
    else:
        default = {5: "level5.jsonl", 6: "level6.jsonl.gz", 7: "level7.jsonl.gz"}[n]
        have = certified(option("--certificates", os.path.join(HERE, default))).get(n, set())
        print(f"candidates that are not kept: {len(rest)}; certificates for {n} points:",
              f"{len(have)}; every candidate is kept or has a certificate: {rest <= have};",
              f"exactly these: {rest == have}")
