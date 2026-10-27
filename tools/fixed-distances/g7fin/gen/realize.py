"""Complex / real realizability tests for a pattern (edge colouring of K_n) via msolve.

Conventions (same as scratchpad/certify5d.py and g7fin/pat_rigid.py):
  p0=(0,0), p1=(x1,0), pi=(xi,yi) for i>=2; the colour of edge (0,1) has value 1, the other
  colour values are variables s1..s_{k-1}; equations |pi-pj|^2 = s_{col(ij)}; Rabinowitsch
  variable t with t*prod_{a<b}(s_a-s_b)*prod_{c>=1} s_c - 1 = 0  (values distinct and nonzero).

Two equivalent formulations of complex realizability (both certified against each other):
  'coords' : the coordinate system above (2n-3 + k-1 + 1 variables).
  'gram'   : the s-variables only: all 3x3 minors of the (n-1)x(n-1) matrix
             C_ij = s_{0i} + s_{0j} - s_{ij} (C_ii = 2 s_{0i}), i.e. rank(Gram) <= 2, plus the same
             Rabinowitsch equation.  This is literally Szollosi-Ostergard's system (their Prop. 4.2
             with d=2), and it is equivalent to 'coords' over C because a complex symmetric matrix of
             rank <= 2 is V^T V with V having 2 rows.
A pattern is complex-realizable iff the ideal is proper (Groebner basis != {1}).  msolve computes the
basis modulo a prime p (the characteristic is written in the input file), so 'discard' means
"1 is in the ideal mod p"; use a second prime to confirm discards if wanted.

All polynomials are produced fully expanded with integer coefficients (msolve silently misparses
anything else).  Temp files carry the pid; never run two msolve jobs on the same file.
"""
import itertools, os, subprocess, re
from collections import defaultdict
from fractions import Fraction
from patlib import EDGES, EIDX

MSOLVE = '/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
PRIMES = [1073741827, 1073741789, 1073741783]   # msolve wants 32-bit primes < 2^31
TMPDIR = os.environ.get('GEN_TMPDIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tmp'))
os.makedirs(TMPDIR, exist_ok=True)

# ---------------------------------------------------------------- tiny sparse polynomial arithmetic
# polynomial = dict {monomial: coeff}, monomial = tuple of (var, exp) sorted by var (var = str)


def pvar(v):
    return {((v, 1),): 1}


def pconst(c):
    return {(): c} if c else {}


def padd(p, q, sign=1):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + sign * c
        if r[m] == 0:
            del r[m]
    return r


def mmul(m1, m2):
    d = dict(m1)
    for v, e in m2:
        d[v] = d.get(v, 0) + e
    return tuple(sorted(d.items()))


def pmul(p, q):
    r = defaultdict(int)
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            r[mmul(m1, m2)] += c1 * c2
    return {m: c for m, c in r.items() if c}


def pstr(p):
    """msolve syntax: expanded, integer coefficients, '^' powers, no spaces."""
    if not p:
        return '0'
    terms = []
    for m, c in sorted(p.items()):
        mono = '*'.join(v if e == 1 else f'{v}^{e}' for v, e in m)
        if mono == '':
            terms.append(str(c))
        elif c == 1:
            terms.append(mono)
        elif c == -1:
            terms.append('-' + mono)
        else:
            terms.append(f'{c}*{mono}')
    s = '+'.join(terms).replace('+-', '-')
    return s


def det3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    t1 = pmul(a, padd(pmul(e, i), pmul(f, h), -1))
    t2 = pmul(b, padd(pmul(d, i), pmul(f, g), -1))
    t3 = pmul(c, padd(pmul(d, h), pmul(e, g), -1))
    return padd(padd(t1, t2, -1), t3)


# ---------------------------------------------------------------- Rabinowitsch polynomial per k
_RAB = {}


def rabinowitsch(k):
    """t*prod_{a<b}(s_a-s_b)*prod_{c>=1}s_c - 1 with s_0 = 1, as a string, for k colours."""
    if k not in _RAB:
        S = [pconst(1)] + [pvar(f's{c}') for c in range(1, k)]
        f = pvar('t')
        for a, b in itertools.combinations(range(k), 2):
            f = pmul(f, padd(S[a], S[b], -1))
        for c in range(1, k):
            f = pmul(f, S[c])
        f = padd(f, pconst(1), -1)
        _RAB[k] = pstr(f)
    return _RAB[k]


# ---------------------------------------------------------------- precomputed Gram minors per n
_MINORS = {}


def gram_minors(n):
    """All 3x3 minors of C (size n-1) as polynomials in edge variables e_{ij} (index into EDGES[n])."""
    if n not in _MINORS:
        eidx = EIDX[n]

        def ev(i, j):
            return pvar(f'e{eidx[(i, j) if i < j else (j, i)]}')
        C = {}
        for i in range(1, n):
            for j in range(1, n):
                if i == j:
                    C[(i, j)] = pmul(pconst(2), ev(0, i))
                else:
                    C[(i, j)] = padd(padd(ev(0, i), ev(0, j)), ev(i, j), -1)
        mins = []
        rows = list(range(1, n))
        for I in itertools.combinations(rows, 3):
            for J in itertools.combinations(rows, 3):
                if J < I:
                    continue   # C symmetric: minor(I,J) = minor(J,I)
                M = [[C[(i, j)] for j in J] for i in I]
                mins.append(det3(M))
        _MINORS[n] = mins
    return _MINORS[n]


def colour_names(col):
    """colour c -> s-variable name (colour of edge (0,1) -> '1')."""
    k = max(col) + 1
    c01 = col[0]
    names = {c01: None}
    nxt = 1
    for c in range(k):
        if c != c01:
            names[c] = f's{nxt}'
            nxt += 1
    return names


_FLAT = {}


def gram_minors_flat(n):
    """Minors of gram_minors(n) flattened to lists of (coeff, (edge,edge,edge)) (every term is a cubic)."""
    if n not in _FLAT:
        flat = []
        for mn in gram_minors(n):
            terms = []
            for m, c in mn.items():
                es = []
                for v, e in m:
                    es.extend([int(v[1:])] * e)
                assert len(es) == 3
                terms.append((c, tuple(es)))
            flat.append(terms)
        _FLAT[n] = flat
    return _FLAT[n]


_MONOSTR = {}


def _monostr(key):
    """key: sorted tuple of colour indices >= 1 (s-variables) -> msolve monomial string."""
    s = _MONOSTR.get(key)
    if s is None:
        parts = []
        i = 0
        while i < len(key):
            j = i
            while j < len(key) and key[j] == key[i]:
                j += 1
            parts.append(f's{key[i]}' if j - i == 1 else f's{key[i]}^{j-i}')
            i = j
        s = '*'.join(parts)
        _MONOSTR[key] = s
    return s


def system_gram(n, col):
    """(variables, polynomial strings) for the Gram-minor formulation."""
    k = max(col) + 1
    names = colour_names(col)
    subs = [0 if names[c] is None else int(names[c][1:]) for c in col]   # edge -> colour index, 0 = value 1
    polys = set()
    for terms in gram_minors_flat(n):
        r = {}
        for c, es in terms:
            key = tuple(sorted(x for x in (subs[es[0]], subs[es[1]], subs[es[2]]) if x))
            r[key] = r.get(key, 0) + c
        items = sorted((kk, c) for kk, c in r.items() if c)
        if not items:
            continue
        out = []
        for kk, c in items:
            if kk == ():
                out.append(str(c))
            elif c == 1:
                out.append(_monostr(kk))
            elif c == -1:
                out.append('-' + _monostr(kk))
            else:
                out.append(f'{c}*{_monostr(kk)}')
        polys.add('+'.join(out).replace('+-', '-'))
    V = [f's{c}' for c in range(1, k)] + ['t']
    return V, sorted(polys) + [rabinowitsch(k)]


def system_gram_slow(n, col):
    """Reference implementation of system_gram (generic polynomial substitution)."""
    k = max(col) + 1
    names = colour_names(col)
    subs = [names[c] for c in col]
    polys = set()
    for mn in gram_minors(n):
        r = defaultdict(int)
        for m, c in mn.items():
            d = defaultdict(int)
            for v, e in m:
                nm = subs[int(v[1:])]
                if nm is not None:
                    d[nm] += e
            r[tuple(sorted(d.items()))] += c
        p = {m: c for m, c in r.items() if c}
        if p:
            polys.add(pstr(p))
    V = [f's{c}' for c in range(1, k)] + ['t']
    return V, sorted(polys) + [rabinowitsch(k)]


def system_coords(n, col):
    """(variables, polynomial strings) for the coordinate formulation."""
    k = max(col) + 1
    names = colour_names(col)
    XV = ['x1'] + [v for i in range(2, n) for v in (f'x{i}', f'y{i}')]
    V = XV + [f's{c}' for c in range(1, k)] + ['t']
    coords = {0: (None, None), 1: ('x1', None)}
    for i in range(2, n):
        coords[i] = (f'x{i}', f'y{i}')
    polys = []
    for e, (i, j) in enumerate(EDGES[n]):
        p = {}
        for a, b in zip(coords[i], coords[j]):
            if a is None and b is None:
                continue
            da = padd(pvar(a) if a else {}, pvar(b) if b else {}, -1)
            p = padd(p, pmul(da, da))
        nm = names[col[e]]
        p = padd(p, pvar(nm) if nm else pconst(1), -1)
        polys.append(pstr(p))
    polys.append(rabinowitsch(k))
    return V, polys


# ---------------------------------------------------------------- msolve drivers
def _files(tag=''):
    pid = os.getpid()   # evaluated per call: safe under multiprocessing/fork
    return (os.path.join(TMPDIR, f'ms_in_{pid}{tag}.txt'), os.path.join(TMPDIR, f'ms_out_{pid}{tag}.txt'))


def msolve_gb_trivial(V, polys, prime=PRIMES[0], timeout=600, tag=''):
    """True iff the reduced Groebner basis of <polys> over F_prime is {1}."""
    fi, fo = _files(tag)
    with open(fi, 'w') as f:
        f.write(','.join(V) + '\n' + str(prime) + '\n' + ',\n'.join(polys) + '\n')
    if os.path.exists(fo):
        os.remove(fo)
    r = subprocess.run([MSOLVE, '-f', fi, '-o', fo, '-g', '2'], capture_output=True, text=True, timeout=timeout)
    out = open(fo).read()
    body = out[out.index('['):].strip() if '[' in out else ''
    if r.returncode != 0 or body == '':
        raise RuntimeError(f'msolve failed: rc={r.returncode} stderr={r.stderr[:300]} out={out[:200]}')
    return body in ('[1]', '[1,]') or body.startswith('[1]')


def complex_realizable(n, col, method='gram', prime=PRIMES[0], **kw):
    V, polys = (system_gram if method == 'gram' else system_coords)(n, col)
    return not msolve_gb_trivial(V, polys, prime=prime, **kw)


def msolve_real(V, polys, precision=200, timeout=3600, tag=''):
    """msolve over Q with real root isolation; returns raw output text."""
    fi, fo = _files(tag)
    with open(fi, 'w') as f:
        f.write(','.join(V) + '\n0\n' + ',\n'.join(polys) + '\n')
    if os.path.exists(fo):
        os.remove(fo)
    subprocess.run([MSOLVE, '-f', fi, '-o', fo, '-p', str(precision)], capture_output=True, text=True, timeout=timeout)
    return open(fo).read()


def real_status(out):
    s = out.replace('\n', '').replace(' ', '')
    if s.startswith('[-1]'):
        return 'empty'
    if s.startswith('[1,'):
        return 'posdim'
    # zero-dimensional: count real solutions
    s2 = s.rstrip(':')
    # output is [0, [nvars, [sol, sol, ...]]] with each sol a list of intervals
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s2)

    def build(i):
        if tokens[i] == '[':
            lst = []
            i += 1
            while tokens[i] != ']':
                v, i = build(i)
                lst.append(v)
                if tokens[i] == ',':
                    i += 1
            return lst, i + 1
        tk = tokens[i]
        if '^' in tk:
            num, den = tk.split('/')
            base, ex = den.split('^')
            return Fraction(int(num), int(base) ** int(ex)), i + 1
        return Fraction(tk), i + 1
    tree, _ = build(0)
    sols = tree[1][1] if isinstance(tree[1], list) and len(tree[1]) > 1 else []
    return 'real' if len(sols) > 0 else 'zerodim_noreal'
