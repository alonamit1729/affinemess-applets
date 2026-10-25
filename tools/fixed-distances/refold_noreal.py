import json
for d in (5,6,7):
    C=json.load(open(f'cert5{d}.json')); byidx={r['idx']:r for r in C}
    for rec in json.load(open(f'noreal2_{d}.json')):
        if rec['idx'] not in byidx: continue
        p=byidx[rec['idx']]['primes'][rec['prime']]
        if rec['verdict']=='real': p['real']=True; p['point']=rec['point']; p['certified_by']='hyperplanes through a numerically found point'
        else: p['real']=False; p['empty_certificate']=rec['certificate']
    json.dump(C,open(f'cert5{d}.json','w'))
print('refolded')
