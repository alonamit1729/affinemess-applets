"""Step 2: random points on every positive-dimensional family of 5-point d-distance sets (saturated ideals in
sat5d.json), via random linear slices solved with msolve exactly as strat.py does; symmetry at each point.
Usage: python3 step2_sample.py d [nslices]   -> writes symtest/samples_5d.json (incremental, resumable)"""
import json, sys, os, re, math, random, subprocess, itertools, time, collections
from os import getpid
from fractions import Fraction
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
PID = str(getpid())
MS = '/usr/local/lib/python3.13/dist-packages/sage_wheels/bin/msolve'
TMP = os.path.join(HERE, 'tmp'); os.makedirs(TMP, exist_ok=True)

def run_msolve(eqs, vars_):
    fi, fo = os.path.join(TMP, f'ms_in_{PID}.txt'), os.path.join(TMP, f'ms_out_{PID}.txt')
    open(fi, 'w').write(','.join(vars_) + '\n0\n' + ',\n'.join(eqs) + '\n')
    subprocess.run([MS, '-f', fi, '-o', fo, '-p', '200'], capture_output=True, text=True, timeout=900)
    return open(fo).read()

def parse_real(out):
    s = out.replace('\n', '').replace(' ', '')
    if s.startswith('[-1]'): return 'empty', []
    if s.startswith('[1,'): return 'posdim', []
    tokens = re.findall(r'-?\d+/2\^\d+|-?\d+/\d+|-?\d+|[\[\],]', s.rstrip(':'))
    if not tokens or tokens[0] != '[': return 'error', []
    def build(i):
        if tokens[i] == '[':
            lst = []; i += 1
            while tokens[i] != ']':
                v, i = build(i); lst.append(v)
                if tokens[i] == ',': i += 1
            return lst, i + 1
        tk = tokens[i]
        if '^' in tk:
            num, den = tk.split('/'); base, ex = den.split('^'); return Fraction(int(num), int(base) ** int(ex)), i + 1
        return Fraction(tk), i + 1
    tree, _ = build(0)
    sols = tree[1][1] if isinstance(tree[1], list) and len(tree[1]) > 1 else []
    return 'zerodim', [[float((lo + hi) / 2) for lo, hi in sol] for sol in sols]

