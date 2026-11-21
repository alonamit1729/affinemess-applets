"""python3 galois_probe.py <jobs dir> <out dir> <recipes dir> [--levels=3,4]
For every rigid card built at one of the given levels (default: regular polygon or equalizer) that verdicts.json
does not cover yet: refine the set to 1300 digits, find the minimal polynomial of each squared-distance ratio and
probe its Galois group by factoring modulo many primes (analysis/fields.py). Appends to <out dir>/<n>_<d>.json
({num: result}); resumable. classify.py turns the probes into verdicts."""
import json, os, sys, glob, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analysis'))
import fields
from multiprocessing import Pool

def work(a):
    cell, num, flat, cls, n = a
    o = dict(cell=cell, i=num, label='', frame='', miss=0, coords=[(flat[k], flat[k + 1]) for k in range(0, len(flat), 2)], cls=cls, n=n)
    r = fields.work(o)
    return cell, num, {k: r[k] for k in ('res', 'ratios', 'err') if k in r}

if __name__ == '__main__':
    jobs_dir, out_dir, rec_dir = sys.argv[1:4]
    levels = {3, 4}
    for a in sys.argv[4:]:
        if a.startswith('--levels='): levels = {int(x) for x in a[9:].split(',')}
    os.makedirs(out_dir, exist_ok=True)
    known = {(c, n) for c, n, v in json.load(open(os.path.join(HERE, '..', 'verdicts.json'))) if v}
    todo = []
    for f in sorted(glob.glob(os.path.join(jobs_dir, '*.json'))):
        cell = os.path.basename(f)[:-5]
        lev = {}
        for g in glob.glob(os.path.join(rec_dir, cell, '*.json')):
            for num, it in json.load(open(g))['items'].items(): lev[int(num)] = it[2]
        done = json.load(open(os.path.join(out_dir, cell + '.json'))) if os.path.exists(os.path.join(out_dir, cell + '.json')) else {}
        for num, kind, dim, flat, cls, label, frame in json.load(open(f)):
            if kind == 'rigid' and lev.get(num) in levels and (cell.replace('_', ','), num) not in known and str(num) not in done:
                todo.append((cell, num, flat, cls, len(flat) // 2))
    # small cells first
    todo.sort(key=lambda a: (sum(map(int, a[0].split('_'))), a[1]))
    print(len(todo), 'sets to probe', flush=True)
    out = {}; t0 = time.time(); k = 0
    def save():
        for cell, D in out.items():
            fn = os.path.join(out_dir, cell + '.json')
            old = json.load(open(fn)) if os.path.exists(fn) else {}
            old.update(D); json.dump(old, open(fn, 'w'))
        out.clear()
    with Pool(2) as P:
        for cell, num, r in P.imap_unordered(work, todo, chunksize=4):
            out.setdefault(cell, {})[str(num)] = r; k += 1
            if k % 200 == 0:
                save(); print(k, '/', len(todo), '%.0fs' % (time.time() - t0), flush=True)
    save(); print('done', k, '%.0fs' % (time.time() - t0), flush=True)
