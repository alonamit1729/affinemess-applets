"""Exact classification of the cell (n, d) from the complete candidate-pattern list g7fin/gen/levels/L{n}.json, in the
Gram (value-space) formulation: the unknowns are the d-1 free squared-distance values s1..s_{d-1} (s0 = 1), the equations
the 3x3 minors of the Gram matrix (rank <= 2), the inequation "values pairwise distinct and nonzero" via a Rabinowitsch
variable t.  A full distance matrix determines the point set up to isometry, so the configuration space of a pattern
(mod similarity) is the set of real value vectors whose Gram matrix is positive semidefinite of rank <= 2 (rank 1 =
collinear); dimensions of families are dimensions of components of that variety.
Per pattern: msolve real roots of (minors + Rabinowitsch) -> 'empty' / zero-dimensional (real solutions -> PSD test ->
rigid configurations) / positive-dimensional -> Singular minAssGTZ of the minor ideal (no t): per prime dim, degeneracy
(an inequation factor vanishes on it), real PSD points by random slices for dim >= 1 (zero-dimensional primes solved
directly), singular locus (slocus) and its real PSD points (hidden rigid candidates).  Checkpointed per pattern.
Usage: python3 cells/cellcert.py n d   -> cells/cert_{n}_{d}.json, cells/cert_{n}_{d}.log"""
import os, sys, json, time, itertools, random, subprocess, re, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); SCR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(SCR, 'g7fin', 'gen')); sys.path.insert(0, os.path.join(SCR, 'g7count'))
os.environ.setdefault('GEN_TMPDIR', os.path.join(SCR, 'tmp'))
import realize
from realize import system_gram, colour_names, MSOLVE
from msolve_util import msolve_real_raw, parse_real
import mpmath as mp
mp.mp.dps = 60
SG = '/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/Singular'
n = int(sys.argv[1]); d = int(sys.argv[2])
PID = os.getpid()
L = json.load(open(os.path.join(SCR, 'g7fin', 'gen', 'levels', f'L{n}.json')))['patterns']
pats = [(i, tuple(p)) for i, p in enumerate(L) if max(p) + 1 == d]
OUT = os.path.join(HERE, f'cert_{n}_{d}.json'); LOG = open(os.path.join(HERE, f'cert_{n}_{d}.log'), 'a')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
E = list(itertools.combinations(range(n), 2))
TOL = mp.mpf('1e-25')

def gram_sets(col, vals):
    """vals: dict s_c -> mpf (c>=1), s0 = 1.  Returns (psd, rank, points) for the Gram matrix of the pattern."""
    names = colour_names(col)
    val = {c: (mp.mpf(1) if names[c] is None else vals[names[c]]) for c in range(d)}
    D = mp.zeros(n, n)
    for e, (i, j) in enumerate(E):
        D[i, j] = D[j, i] = val[col[e]]
    G = mp.zeros(n - 1, n - 1)
    for i in range(1, n):
        for j in range(1, n):
            G[i - 1, j - 1] = (D[0, i] + D[0, j] - D[i, j]) / 2
    ev, U = mp.eigsy(G)
    evs = sorted([ev[i] for i in range(n - 1)])
    psd = evs[0] > -TOL and (n - 1 < 3 or evs[-3] < TOL)
    rank = sum(1 for x in evs if x > TOL)
    if not psd: return False, rank, None
    idx = sorted(range(n - 1), key=lambda i: -ev[i])[:2]
    pts = [(mp.mpf(0), mp.mpf(0))] + [tuple(mp.sqrt(max(ev[k], mp.mpf(0))) * U[i, k] for k in idx) for i in range(n - 1)]
    return True, rank, [[float(p[0]), float(p[1])] for p in pts]

def sol_vals(V, sol):
    return {name: (mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator) / 2 for name, (a, b) in zip(V, sol)}

def congruent(P, Q, tol=1e-7):
    D1 = [[math.dist(a, b) for b in P] for a in P]; D2 = [[math.dist(a, b) for b in Q] for a in Q]
    m1 = min(x for r in D1 for x in r if x > 1e-12); m2 = min(x for r in D2 for x in r if x > 1e-12)
    D1 = [[x / m1 for x in r] for r in D1]; D2 = [[x / m2 for x in r] for r in D2]
    T = tol * max(max(r) for r in D2)
    r1 = [sorted(r) for r in D1]; r2 = [sorted(r) for r in D2]
    cand = [[j for j in range(n) if all(abs(a - b) < T for a, b in zip(r1[i], r2[j]))] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i])); perm = [-1] * n; used = [False] * n
    def bt(k):
        if k == n: return True
        i = order[k]
        for j in cand[i]:
            if used[j]: continue
            if all(abs(D1[i][order[t]] - D2[j][perm[order[t]]]) < T for t in range(k)):
                perm[i] = j; used[j] = True
                if bt(k + 1): return True
                used[j] = False; perm[i] = -1
        return False
    return bt(0)

def linform(SV):
    """random integer linear form, in msolve syntax (no parentheses, no leading '+')"""
    terms = []
    for v in SV:
        a = random.randint(-7, 7)
        if a: terms.append(f'{a}*{v}')
    b = random.randint(-9, 9)
    terms.append(str(b))
    return '+'.join(terms).replace('+-', '-')

def singular(script, timeout=600):
    fn = os.path.join(SCR, 'tmp', f'cc_{PID}.sing'); open(fn, 'w').write(script)
    r = subprocess.run([SG, '-q', fn], capture_output=True, text=True, timeout=timeout)
    return [l.strip() for l in r.stdout.splitlines() if l.strip() and not l.startswith('//')]

def _sing_primes(ring, gens, algo, tmo):
    """minimal associated primes of the ideal <gens> by `algo`; list of generator-lists (strings)"""
    script = (ring + '\nLIB "primdec.lib";\nideal i=' + ',\n'.join(gens) + ';\nlist L=' + algo + '(i);\nint k; for(k=1;k<=size(L);k++){ print("PRIME"); print(L[k]); }\nquit;\n')
    out = singular(script, timeout=tmo); res = []
    for l in out:
        if l == 'PRIME': res.append([])
        elif res: res[-1].append((l.split('=', 1)[1].strip() if '=' in l else l).rstrip(','))
    return res

def decompose(SV, polys, facs):
    """minimal primes: minAssGTZ on the whole ideal (20 s); failing that, split by facstd and decompose each piece with
    minAssGTZ (60 s) or else minAssChar (60 s) -- on the hard patterns one of the two always finishes in seconds on each
    piece (they stall on different pieces) -- then discard non-minimal primes."""
    ring = 'ring r=0,(' + ','.join(SV) + '),dp;'
    pre = ring + '\nLIB "primdec.lib"; LIB "sing.lib";\n'
    post = ('int k; ideal G; ideal S; int bad;\n'
            'for(k=1;k<=size(L);k++){ G=std(L[k]); print("PRIME"); print(dim(G)); bad=0; ' + ' '.join(f'if(reduce({f},G)==0){{bad=1;}}' for f in facs) +
            ' print(bad); S=std(slocus(L[k])); print(dim(S)); print("GENS"); print(L[k]); print("SLOCUS"); print(S); }\nquit;\n')
    try:
        lines = singular(pre + 'ideal i=' + ',\n'.join(polys) + ';\nlist L=minAssGTZ(i);\n' + post, timeout=20)
    except subprocess.TimeoutExpired:
        log('   minAssGTZ timeout (20s): splitting with facstd')
        rring = 'ring r=0,(' + ','.join(SV[::-1]) + '),dp;'   # reversed variable order: a different term order often makes the stalled step finish at once
        allp = None
        for rg in (ring, rring):
            try:
                sp = singular(rg + '\nLIB "primdec.lib";\nideal i=' + ',\n'.join(polys) + ';\nlist L=facstd(i); int k; for(k=1;k<=size(L);k++){ print("PIECE"); print(L[k]); }\nquit;\n', timeout=300)
                pieces = []
                for l in sp:
                    if l == 'PIECE': pieces.append([])
                    elif pieces: pieces[-1].append((l.split('=', 1)[1].strip() if '=' in l else l).rstrip(','))
                cur_p = []
                for pc in pieces:
                    for algo in ('minAssGTZ', 'minAssChar'):
                        try:
                            cur_p += _sing_primes(rg, pc, algo, 60); break
                        except subprocess.TimeoutExpired:
                            log(f'   piece: {algo} timeout ({"reversed" if rg is rring else "plain"} order)')
                    else:
                        raise subprocess.TimeoutExpired('piece', 0)
                allp = cur_p; break
            except subprocess.TimeoutExpired:
                log(f'   split/decompose failed in {"reversed" if rg is rring else "plain"} order')
        if allp is None: raise subprocess.TimeoutExpired('pieces', 0)
        decl = 'list P;\n' + ''.join(f'P[{k+1}]=ideal({",".join(g)});\n' for k, g in enumerate(allp))
        filt = ('list L; int a; int b; int keep; ideal Pa; ideal Pb;\n'
                'for(a=1;a<=size(P);a++){ keep=1; Pa=modStd(P[a]); for(b=1;b<=size(P);b++){ if(b!=a){ Pb=modStd(P[b]); '
                'if(size(reduce(Pb,Pa))==0){ if(size(reduce(Pa,Pb))!=0 || b<a){ keep=0; } } } } if(keep){ L[size(L)+1]=P[a]; } }\n')
        post2 = ('int k; ideal G; ideal S; ideal M; matrix Jm; int c; int bad;\n'
                 'for(k=1;k<=size(L);k++){ G=modStd(L[k]); print("PRIME"); print(dim(G)); bad=0; ' + ' '.join(f'if(reduce({f},G)==0){{bad=1;}}' for f in facs) +
                 ' print(bad); c=' + str(len(SV)) + '-dim(G); Jm=jacob(G); M=minor(Jm,c,G); S=modStd(G+M); print(dim(S)); print("GENS"); print(L[k]); print("SLOCUS"); print(S); }\nquit;\n')
        lines = singular(pre.replace('LIB "sing.lib";', 'LIB "sing.lib"; LIB "modstd.lib";') + decl + filt + post2, timeout=600)
    primes = []; cur = None; mode = None
    for line in lines:
        if line == 'PRIME': cur = {'dim': None, 'degenerate': None, 'sdim': None, 'gens': [], 'slocus': []}; primes.append(cur); mode = None; continue
        if cur is None: continue
        if line == 'GENS': mode = 'gens'; continue
        if line == 'SLOCUS': mode = 'slocus'; continue
        if mode is None:
            if cur['dim'] is None and re.fullmatch(r'-?\d+', line): cur['dim'] = int(line); continue
            if cur['degenerate'] is None and re.fullmatch(r'[01]', line): cur['degenerate'] = bool(int(line)); continue
            if cur['sdim'] is None and re.fullmatch(r'-?\d+', line): cur['sdim'] = int(line); continue
        else:
            g = line.split('=', 1)[1].strip() if '=' in line else line
            cur[mode].append(g.rstrip(','))
    # msolve needs integer coefficients: clear denominators of every generator
    import sympy as sp
    syms = {v: sp.Symbol(v) for v in SV}
    def integral(g):
        if '/' not in g: return g
        e = sp.sympify(g.replace('^', '**'), locals=syms)
        _, p = sp.Poly(e, *[syms[v] for v in SV]).clear_denoms()
        return str(sp.expand(p.as_expr())).replace('**', '^').replace(' ', '')
    for cur in primes:
        cur['gens'] = [integral(g) for g in cur['gens']]; cur['slocus'] = [integral(g) for g in cur['slocus']]
    return primes

results = json.load(open(OUT)) if os.path.exists(OUT) else {}
t0 = time.time()
for pi, (idx, col) in enumerate(pats):
    if str(idx) in results: continue
    V, polys = system_gram(n, col)
    SV = V[:-1]; rab = polys[-1]; minors = polys[:-1]
    k = d
    facs = [f's{a}-s{b}' for a, b in itertools.combinations(range(1, k), 2)] + [f's{c}-1' for c in range(1, k)] + [f's{c}' for c in range(1, k)]
    rec = {'idx': idx, 'col': list(col), 'rigid': [], 'nonpsd': 0, 'primes': []}
    out = msolve_real_raw(V, polys, precision=200, timeout=3600)
    sols = parse_real(out)
    if sols == 'empty':
        rec['status'] = 'empty'
    elif sols != 'posdim':
        rec['status'] = 'zerodim'; rec['nsol'] = len(sols)
        for sol in sols:
            psd, rank, P = gram_sets(col, sol_vals(V, sol))
            if psd:
                if not any(congruent(P, Q) for Q, _ in rec['rigid']): rec['rigid'].append((P, rank))
            else: rec['nonpsd'] += 1
    else:
        rec['status'] = 'posdim'
        try:
            primes = decompose(SV, minors, facs)
        except subprocess.TimeoutExpired:
            rec['status'] = 'posdim_unresolved'; primes = []
            log(f'{pi+1}/{len(pats)} pat {idx} decomposition TIMEOUT')
        random.seed(idx)
        for p in primes:
            q = {'dim': p['dim'], 'degenerate': p['degenerate'], 'slocus_dim': p['sdim']}
            if not p['degenerate']:
                if p['dim'] == 0:
                    s2 = parse_real(msolve_real_raw(V, p['gens'] + [rab], precision=200))
                    q['nsol'] = 0 if s2 in ('empty', 'posdim') else len(s2); q['psd'] = 0
                    for sol in (s2 if isinstance(s2, list) else []):
                        psd, rank, P = gram_sets(col, sol_vals(V, sol))
                        if psd:
                            q['psd'] += 1
                            if not any(congruent(P, Q) for Q, _ in rec['rigid']): rec['rigid'].append((P, rank))
                else:
                    hits = 0; ranks = set()
                    for _ in range(40):
                        lins = [linform(SV) for _ in range(p['dim'])]
                        s2 = parse_real(msolve_real_raw(V, p['gens'] + lins + [rab], precision=200))
                        for sol in (s2 if isinstance(s2, list) else []):
                            psd, rank, P = gram_sets(col, sol_vals(V, sol))
                            if psd:
                                hits += 1; ranks.add(rank)
                                if 'point' not in q: q['point'] = P
                                break
                        if hits >= 2: break
                    q['real_psd'] = hits > 0; q['slice_hits'] = hits; q['ranks'] = sorted(ranks)
                    # hidden rigid candidates: real PSD points of the singular locus
                    if p['sdim'] >= 0:
                        s3 = parse_real(msolve_real_raw(V, p['slocus'] + [rab], precision=200))
                        if s3 == 'posdim': q['singular_real'] = 'posdim'
                        elif s3 == 'empty': q['singular_real'] = 0
                        else:
                            cnt = 0
                            for sol in s3:
                                psd, rank, P = gram_sets(col, sol_vals(V, sol))
                                if psd: cnt += 1; q.setdefault('singular_points', []).append(P)
                            q['singular_real'] = cnt
                    else: q['singular_real'] = 0
            rec['primes'].append(q)
    results[str(idx)] = rec
    tmp = OUT + f'.tmp{PID}'; json.dump(results, open(tmp, 'w')); os.replace(tmp, OUT)
    fams = [(q['dim'], q.get('real_psd'), q.get('singular_real')) for q in rec['primes'] if not q['degenerate'] and q['dim'] > 0]
    log(f'{pi+1}/{len(pats)} pat {idx} {rec["status"]} rigid {len(rec["rigid"])} nonpsd {rec["nonpsd"]} fam {fams} {time.time()-t0:.0f}s')
log('done')
