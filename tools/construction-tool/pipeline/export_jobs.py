"""python3 export_jobs.py <jobs dir>
   Export every catalog card (rigid + families) as construction jobs: jobs/<n>_<d>.json
   items: [num, kind, dim, flat coords (refined for rigid), cls, label]"""
import sys, json, os
import os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..')
sys.path.insert(0,R+'/fixed-distances'); import datautil
sys.path.insert(0,R+'/construction-tool'); import refine
from multiprocessing import Pool
OUT=sys.argv[1] if len(sys.argv)>1 else 'jobs'
os.makedirs(OUT,exist_ok=True)
D=datautil.load(R+'/fixed-distances/data/applet_data.json')
def job(a):
    cell,num,it=a
    coords=it['coords']
    if True:   # refine every set (families too) onto its distance pattern
        try:
            z,res=refine.refine(coords,it['cls'])
            if res>1e-9: z=[complex(float(x),float(y)) for x,y in coords]
        except Exception: z=[complex(float(x),float(y)) for x,y in coords]
    else: z=[complex(float(x),float(y)) for x,y in coords]
    return cell,[num,it['kind'],it.get('dim',0),[round(v,13) for q in z for v in (q.real,q.imag)],it['cls'],it['label'],it.get('frame') or '']
jobs=[(cell,k+1,it) for cell,c in D['cats'].items() for k,it in enumerate(c['items'])]
out={}
with Pool(2) as P:
    for cell,r in P.imap(job,jobs,chunksize=200): out.setdefault(cell,[]).append(r)
for cell,L in out.items():
    json.dump(L,open(OUT+'/'+cell.replace(',','_')+'.json','w'),separators=(',',':'))
print(len(jobs), sum(len(v) for v in out.values()))
