"""Settle the curves of the (6,6) classification on which random slices found no real point.

cellcert.py 6 6 left 19 non-degenerate one-dimensional primes (one in each of 19 patterns) with no real PSD point
found on 40 random slices; a miss is not a proof.  Here each such pattern is decomposed again (the same Singular
routine as cellcert.py) and every non-degenerate positive-dimensional prime P is tested exactly:

  Critical-point test.  Let C be the real points of V(P), a closed set in R^5 (P is prime of dimension 1, so its
  complex curve has finitely many singular points).  Pick a rational point A and four random rational combinations
  g1..g4 of the generators of P, and let W = P + (det [grad g1; ...; grad g4; x - A]).  If C is not empty, the
  squared distance to A attains its minimum on C (C is closed and the distance is proper) at a point p.  Either p is
  a singular point of V(P) (gradients of P of rank <= 3, so the determinant vanishes), or the gradients of P span the
  normal space at p and x - A is normal there, so x - A lies in the span of all gradients of P; if grad g1..g4 span
  that space the determinant vanishes, and if they do not, it vanishes as well.  So p lies in W.  When W is
  zero-dimensional (checked), msolve isolates all its real points.

  Positivity.  C usually has real points; what matters is K = the points of C whose Gram matrix is positive
  semidefinite (the real configurations).  On C the Gram matrix has rank <= 2, so its eigenvalues are la, lb, 0, 0, 0
  and K = {la >= 0, lb >= 0} is closed in C.  If K is not empty, take a connected component G of C meeting K.  Either
  K contains G, and then the minimiser of the distance to A on G is a point of W in K; or K has a boundary point in
  G, where one of la, lb changes sign, so the Gram matrix there has rank <= 1: a real point of
  R1 = P + (2x2 minors of the Gram matrix) in K.  So: if no real point of W and no real point of R1 is positive
  semidefinite, then K is empty and the prime carries no configuration.

  Sign certificate.  When the positivity test is inconclusive (a positive semidefinite point exists, typically a
  degenerate one), look for an element of P whose coefficients all have the same sign (the constant term stands for
  the unit value s0 = 1): at a real configuration every value is positive, so such an element cannot vanish there.
  Candidates are the elements of reduced Groebner bases of P for several lexicographic orders.

Writes cells/settle66.json.   cd tools/fixed-distances && python3 cells/settle66.py   (needs Singular and msolve)
"""
import os, sys, json, random, time, itertools, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); SCR = os.path.dirname(HERE)
os.makedirs(os.path.join(SCR, 'tmp'), exist_ok=True)

# reuse cellcert.py's definitions (everything before its main loop), for n = 6, d = 6
src = open(os.path.join(HERE, 'cellcert.py')).read()
src = src[:src.index('results = json.load(open(OUT))')]
src = src.replace("L = json.load(open(os.path.join(SCR, 'g7fin', 'gen', 'levels', f'L{n}.json')))['patterns']", "L = []")
src = src.replace("LOG = open(os.path.join(HERE, f'cert_{n}_{d}.log'), 'a')", "LOG = open(os.devnull, 'w')")
sys.argv = ['cellcert.py', '6', '6']
G = {'__file__': os.path.join(HERE, 'cellcert.py'), '__name__': 'cellcert'}
exec(compile(src, 'cellcert.py', 'exec'), G)
system_gram, decompose, singular = G['system_gram'], G['decompose'], G['singular']
msolve_real_raw, parse_real, gram_sets, sol_vals = G['msolve_real_raw'], G['parse_real'], G['gram_sets'], G['sol_vals']

n = d = 6


def critical_ideal(SV, gens, rng):
    """generators of W = P + det[grad g1..g4; x - A] (integer coefficients) and dim(W), via Singular"""
    k = len(SV)
    A = [rng.randint(-23, 23) for _ in SV]
    combos = [[rng.randint(-9, 9) for _ in gens] for _ in range(k - 1)]
    lines = ['ring r=0,(' + ','.join(SV) + '),dp;', 'LIB "modstd.lib";',
             'ideal P=' + ',\n'.join(gens) + ';', f'matrix M[{k}][{k}];', 'poly g; int j;']
    for i, c in enumerate(combos):
        lines.append('g=' + '+'.join(f'({a})*P[{t + 1}]' for t, a in enumerate(c) if a) + ';')
        lines.append(f'for(j=1;j<={k};j++){{ M[{i + 1},j]=diff(g,var(j)); }}')
    for j, v in enumerate(SV):
        lines.append(f'M[{k},{j + 1}]={v}-({A[j]});')
    lines += ['ideal W=P+ideal(det(M));', 'ideal S=modStd(W);', 'print("DIM"); print(dim(S));',
              'print("GENS"); print(S);', 'quit;']
    out = singular('\n'.join(lines) + '\n', timeout=1800)
    dim = int(out[out.index('DIM') + 1]); gl = out[out.index('GENS') + 1:]
    W = [(l.split('=', 1)[1].strip() if '=' in l else l).rstrip(',') for l in gl]
    import sympy as sp
    syms = {v: sp.Symbol(v) for v in SV}
    def integral(g):
        if '/' not in g: return g
        _, p = sp.Poly(sp.sympify(g.replace('^', '**'), locals=syms), *syms.values()).clear_denoms()
        return str(sp.expand(p.as_expr())).replace('**', '^').replace(' ', '')
    return A, combos, dim, [integral(g) for g in W]


def rank1_minors(col):
    """the 2x2 minors of the Gram matrix C_ij = (D_0i + D_0j - D_ij)/2, i, j = 1..5, as polynomial strings"""
    import sympy as sp
    names = G['colour_names'](col)
    E = list(itertools.combinations(range(n), 2))
    val = lambda c: sp.Integer(1) if names[c] is None else sp.Symbol(names[c])
    D = [[sp.Integer(0)] * n for _ in range(n)]
    for e, (i, j) in enumerate(E):
        D[i][j] = D[j][i] = val(col[e])
    Cm = sp.Matrix(n - 1, n - 1, lambda i, j: (D[0][i + 1] + D[0][j + 1] - D[i + 1][j + 1]))   # 2 x Gram
    out = set()
    for r in itertools.combinations(range(n - 1), 2):
        for c in itertools.combinations(range(n - 1), 2):
            m = sp.expand(Cm[r[0], c[0]] * Cm[r[1], c[1]] - Cm[r[0], c[1]] * Cm[r[1], c[0]])
            if m != 0: out.add(str(m).replace('**', '^').replace(' ', ''))
    return sorted(out)


def sign_certificate(SV, gens):
    """an element of <gens> with all coefficients of one sign, from reduced lex Groebner bases, or None"""
    import sympy as sp
    syms = [sp.Symbol(v) for v in SV]
    for r in range(len(SV)):
        order = SV[r:] + SV[:r]
        for rev in (False, True):
            o = order[::-1] if rev else order
            try:
                out = singular('ring r=0,(' + ','.join(o) + '),lp;\nideal P=' + ',\n'.join(gens) + ';\noption(redSB);\n'
                               'ideal S=std(P);\nprint("GENS"); print(S);\nquit;\n', timeout=60)
            except subprocess.TimeoutExpired:
                subprocess.run(['pkill', '-f', 'cc_%d.sing' % os.getpid()]); continue
            for l in out[out.index('GENS') + 1:]:
                g = (l.split('=', 1)[1].strip() if '=' in l else l).rstrip(',')
                coeffs = sp.Poly(sp.sympify(g.replace('^', '**'), locals={v: x for v, x in zip(SV, syms)}), *syms).coeffs()
                if coeffs and (all(c > 0 for c in coeffs) or all(c < 0 for c in coeffs)):
                    return g
    return None


def main():
    C = json.load(open(os.path.join(HERE, 'cert_6_6.json')))
    todo = [k for k, r in C.items() if any(q.get('dim', 0) >= 1 and not q.get('degenerate') and not q.get('real_psd')
                                           for q in r.get('primes', []))]
    OUT = os.path.join(HERE, 'settle66.json')
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    print(len(todo), 'patterns', flush=True)
    for key, out in res.items():
        # only where the exact tests prove fewer primes empty than cellcert.py found without real points
        bad = sum(1 for dim, deg, real in out['cellcert_primes'] if dim >= 1 and not deg and real is False)
        proved = sum(1 for q in out['primes'] if q.get('no_configuration'))
        if proved >= bad: continue
        changed = False
        for q in sorted((q for q in out['primes'] if q.get('no_configuration') is False), key=lambda q: q.get('psd_points', 0)):
            if 'sign_certificate' not in q:
                q['sign_certificate'] = sign_certificate(['s%d' % i for i in range(1, d)], q['gens']); changed = True
                if q['sign_certificate']: break
        if changed:
            print(key, 'sign certificates', [q.get('sign_certificate') for q in out['primes']], flush=True)
            json.dump(res, open(OUT + '.tmp', 'w')); os.replace(OUT + '.tmp', OUT)
    for key in todo:
        if key in res: continue
        t0 = time.time(); rec = C[key]; col = rec['col']
        V, polys = system_gram(n, col); SV = V[:-1]; minors = polys[:-1]
        facs = [f's{a}-s{b}' for a, b in itertools.combinations(range(1, d), 2)] + [f's{c}-1' for c in range(1, d)] + [f's{c}' for c in range(1, d)]
        primes = decompose(SV, minors, facs)
        out = {'col': col, 'cellcert_primes': [(q['dim'], q['degenerate'], q.get('real_psd')) for q in rec['primes']],
               'primes': []}
        for p in primes:
            q = {'dim': p['dim'], 'degenerate': p['degenerate']}
            if p['dim'] >= 1 and not p['degenerate']:
                q['gens'] = p['gens']
                rng = random.Random(int(key) * 7919 + len(out['primes']))
                for attempt in range(4):
                    A, combos, wdim, W = critical_ideal(SV, p['gens'], rng)
                    if wdim == 0: break
                q.update({'A': A, 'combos': combos, 'W_dim': wdim})
                if wdim == 0:
                    sols = parse_real(msolve_real_raw(SV, W, precision=128, timeout=1800))
                    if sols == 'empty' or sols == []:
                        q['real_points'] = 0
                    elif sols == 'posdim':
                        q['real_points'] = 'posdim?'
                    else:
                        q['real_points'] = len(sols)
                        q['psd_points'] = sum(1 for s in sols if gram_sets(col, sol_vals(SV, s))[0])
                    r1 = parse_real(msolve_real_raw(SV, p['gens'] + rank1_minors(col), precision=128, timeout=1800))
                    if r1 == 'posdim': q['R1_real'] = 'posdim'
                    else:
                        r1 = [] if r1 == 'empty' else r1
                        q['R1_real'] = len(r1)
                        q['R1_psd'] = sum(1 for s in r1 if gram_sets(col, sol_vals(SV, s))[0])
                    q['no_configuration'] = (q.get('psd_points') in (0, None) and q.get('real_points') != 'posdim?'
                                             and q.get('R1_real') != 'posdim' and q.get('R1_psd', 0) == 0)
                    if not q['no_configuration']:
                        q['sign_certificate'] = sign_certificate(SV, p['gens'])
            out['primes'].append(q)
        out['seconds'] = round(time.time() - t0, 1)
        res[key] = out
        json.dump(res, open(OUT + '.tmp', 'w')); os.replace(OUT + '.tmp', OUT)
        print(key, [(q['dim'], q['degenerate'], q.get('W_dim'), q.get('real_points'), q.get('psd_points'), q.get('R1_real'), q.get('R1_psd'), q.get('no_configuration'), q.get('sign_certificate')) for q in out['primes']],
              out['cellcert_primes'], out['seconds'], 's', flush=True)


if __name__ == '__main__':
    main()
