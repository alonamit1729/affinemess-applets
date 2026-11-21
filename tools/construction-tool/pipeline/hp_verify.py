"""python3 hp_verify.py <recipes dir> [cells...]
Replay every stored construction without equalities (lattice, straightedge and compass, regular polygon; families
with free points at their stored positions) in 50-digit arithmetic and check that its points have exactly d
distinct distances, equal ones agreeing to 1e-35 and different ones differing by more than 1e-20 (relative).
The page and the search work in doubles, where a near miss (a construction landing within 1e-9 of the set
without being exact) cannot be told from the real thing. Prints the misses; --json writes them as a list."""
import json, os, sys, glob
import mpmath as mp
mp.mp.dps = 50

def P(V, i):
    v = V[i]; return v if isinstance(v, tuple) and len(v) == 2 else None

def latvec(kind, a, b):
    ux, uy = b[0] - a[0], b[1] - a[1]
    if kind == 'tri':
        c, s = mp.mpf(1) / 2, mp.sqrt(3) / 2
        return ux, uy, c * ux - s * uy, s * ux + c * uy
    return ux, uy, -uy, ux

def curve(n, V):
    t = n[0]
    if t in ('l', 'bis', 'par', 'perp'):
        if t in ('l', 'bis'): a, b = P(V, n[1]), P(V, n[2])
        else: a, b = P(V, n[1]), P(V, n[2])
        dx, dy = b[0] - a[0], b[1] - a[1]; L = mp.sqrt(dx * dx + dy * dy); dx, dy = dx / L, dy / L
        if t == 'l': return ('L', a[0], a[1], dx, dy)
        if t == 'bis': return ('L', (a[0] + b[0]) / 2, (a[1] + b[1]) / 2, -dy, dx)
        p = P(V, n[3])
        return ('L', p[0], p[1], dx, dy) if t == 'par' else ('L', p[0], p[1], -dy, dx)
    if t == 'c':
        c, p = P(V, n[1]), P(V, n[2]); return ('C', c[0], c[1], mp.sqrt((p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2))
    if t == 'k':
        p, q, c = P(V, n[1]), P(V, n[2]), P(V, n[3]); return ('C', c[0], c[1], mp.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2))
    return None

def meet(A, B):
    """the same branch order as the page's meet(): line-line one point; line-circle and circle-circle two, ordered as there"""
    if A[0] == 'L' and B[0] == 'L':
        _, px, py, dx, dy = A; _, qx, qy, ex, ey = B
        den = dx * ey - dy * ex
        t = ((qx - px) * ey - (qy - py) * ex) / den
        return [(px + t * dx, py + t * dy)]
    if A[0] == 'C' and B[0] == 'L': A, B = B, A
    if A[0] == 'L':
        _, px, py, dx, dy = A; _, cx, cy, r = B
        fx, fy = px - cx, py - cy
        b = fx * dx + fy * dy; c = fx * fx + fy * fy - r * r; h = b * b - c
        if h < 0: h = mp.mpf(0)
        s = mp.sqrt(h)
        return [(px + (-b - s) * dx, py + (-b - s) * dy), (px + (-b + s) * dx, py + (-b + s) * dy)]
    _, ax, ay, r1 = A; _, bx, by, r2 = B
    dx, dy = bx - ax, by - ay; d = mp.sqrt(dx * dx + dy * dy)
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d); h2 = r1 * r1 - a * a
    if h2 < 0: h2 = mp.mpf(0)
    h = mp.sqrt(h2); mx, my = ax + a * dx / d, ay + a * dy / d
    return [(mx - h * dy / d, my + h * dx / d), (mx + h * dy / d, my - h * dx / d)]   # left of c1->c2 first, as in the page

def evaluate(nodes):
    V = []
    for n in nodes:
        t = n[0]
        if t in ('b', 'f'): V.append((mp.mpf(n[1]), mp.mpf(n[2])))
        elif t in ('t', 's'):
            a, b = P(V, n[1]), P(V, n[2]); ux, uy, vx, vy = latvec('tri' if t == 't' else 'sq', a, b)
            V.append((a[0] + n[3] * ux + n[4] * vx, a[1] + n[3] * uy + n[4] * vy))
        elif t == 'poly':
            a, b, m = P(V, n[1]), P(V, n[2]), n[3]
            w = mp.expjpi(mp.mpf(2) / m); A, B = mp.mpc(*a), mp.mpc(*b)
            c = (B - w * A) / (1 - w)
            V.append(('poly', c, A - c, w))
        elif t == 'pv':
            _, c, r, w = V[n[1]]; z = c + r * w ** n[2]; V.append((z.real, z.imag))
        elif t == 'pc':
            c = V[n[1]][1]; V.append((c.real, c.imag))
        elif t == 'g':
            C = V[n[1]]; u = mp.mpf(n[2])
            V.append((C[1] + u * C[3], C[2] + u * C[4]) if C[0] == 'L' else (C[1] + C[3] * mp.cos(u), C[2] + C[3] * mp.sin(u)))
        elif t == 'x':
            m = meet(V[n[1]], V[n[2]]); V.append(m[min(n[3], len(m) - 1)])
        else: V.append(curve(n, V))
    return V

def check(rec, d):
    V = evaluate(rec['n']); pts = [V[i] for i in rec['m']]
    D = sorted(mp.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) for i, p in enumerate(pts) for q in pts[i + 1:])
    classes = 1; worst_eq = mp.mpf(0); best_gap = mp.inf
    for a, b in zip(D, D[1:]):
        rel = (b - a) / b
        if rel < mp.mpf(10) ** -20: worst_eq = max(worst_eq, rel)
        else: classes += 1; best_gap = min(best_gap, rel)
    return classes == d and worst_eq < mp.mpf(10) ** -35, classes, worst_eq

if __name__ == '__main__':
    rec_dir = sys.argv[1]; cells = [a for a in sys.argv[2:] if not a.startswith('--')]
    bad = []; k = 0
    for cdir in sorted(glob.glob(os.path.join(rec_dir, '*_*'))):
        cell = os.path.basename(cdir)
        if cells and cell not in cells: continue
        d = int(cell.split('_')[1])
        for f in glob.glob(os.path.join(cdir, '[0-9]*.json')):
            for num, it in json.load(open(f))['items'].items():
                if it[2] in (1, 2, 3):
                    k += 1
                    ok, cl, eq = check(it[6], d)
                    if not ok: bad.append([cell.replace('_', ','), int(num), it[2], cl, mp.nstr(eq, 3)])
    print(k, 'constructions replayed at 50 digits;', len(bad), 'not exact')
    for b in bad[:30]: print(b)
    if '--json' in sys.argv: json.dump(bad, open('/tmp/hp_bad.json', 'w'))
