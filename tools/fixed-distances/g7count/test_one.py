import sys, json, time
from msolve_util import *
L = json.load(open(os.path.join(SCR,'g7fin','gen','levels','L10.json')))['patterns']
idx = int(sys.argv[1]); prec = int(sys.argv[2])
V, sols, dt = realize_pattern(10, tuple(L[idx]), precision=prec)
print(idx, 'nsol', len(sols) if isinstance(sols, list) else sols, 'time', round(dt,1))
if isinstance(sols, list):
    w = max(float(b-a) for sol in sols for a,b in sol)
    print('max interval width', w)
