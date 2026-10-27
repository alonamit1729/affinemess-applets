"""msolve real-root driver for 10-point patterns (coordinate system of g7fin/gen/realize.py)."""
import os, sys, re, subprocess, json, time
from fractions import Fraction
SCR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SCR, 'g7fin', 'gen'))
os.environ.setdefault('GEN_TMPDIR', os.path.join(SCR, 'tmp'))
import realize
from realize import system_coords, MSOLVE

def msolve_real_raw(V, polys, precision=400, timeout=7200, tag=''):
    pid = os.getpid()
    fi = os.path.join(SCR, 'tmp', f'g7c_in_{pid}{tag}.txt')
    fo = os.path.join(SCR, 'tmp', f'g7c_out_{pid}{tag}.txt')
    with open(fi, 'w') as f:
        f.write(','.join(V) + '\n0\n' + ',\n'.join(polys) + '\n')
    if os.path.exists(fo): os.remove(fo)
    r = subprocess.run([MSOLVE, '-f', fi, '-o', fo, '-p', str(precision)], capture_output=True, text=True, timeout=timeout)
    out = open(fo).read()
    os.remove(fi); os.remove(fo)
    return out

def parse_real(out):
    """-> 'empty' | 'posdim' | list of solutions, each a list of (lo,hi) Fractions per variable."""
    s = out.replace('\n', '').replace(' ', '').rstrip(':')
    if s.startswith('[-1'): return 'empty'
    if s.startswith('[1,'): return 'posdim'
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s)
    def build(i):
        if tokens[i] == '[':
            lst = []; i += 1
            while tokens[i] != ']':
                v, i = build(i); lst.append(v)
                if tokens[i] == ',': i += 1
            return lst, i + 1
        tk = tokens[i]
        if '^' in tk:
            num, den = tk.split('/'); base, ex = den.split('^')
            return Fraction(int(num), int(base) ** int(ex)), i + 1
        return Fraction(tk), i + 1
    tree, _ = build(0)
    assert tree[0] == 0, tree[0]
    sols = tree[1][1]
    return [[(a, b) for a, b in sol] for sol in sols]

def realize_pattern(n, col, precision=400, timeout=7200):
    V, P = system_coords(n, col)
    t0 = time.time()
    out = msolve_real_raw(V, P, precision=precision, timeout=timeout)
    sols = parse_real(out)
    return V, sols, time.time() - t0
