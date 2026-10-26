"""Saturated ideals (closure of the non-degenerate locus) for the positive-dimensional types of (5,d), with x1 = 1.
Usage: python3 sat5.py d [idx ...]; writes sat5{d}.json incrementally."""
import sys, json, os, subprocess, time
src=open('certify_nd.py').read().split('import os')[0]
exec(src)   # defines classes, system, facs, V, SG, PID, expand, ...
n=int(sys.argv[1]); d=int(sys.argv[2]); 
out={}
fn=f'sat{n}{d}.json'
if os.path.exists(fn): out=json.load(open(fn))
R=json.load(open(f'cert{n}{d}.json'))
want=[r['idx'] for r in R if r.get('unresolved') or any((q.get('real') or q.get('real_slices')) and q['dim']>0 for q in r['primes'])]
if len(sys.argv)>3: want=[int(x) for x in sys.argv[3:]]
t0=time.time()
for ci in want:
    if str(ci) in out: continue
    eqs=system(classes[ci])+['x1-1']
    script=('ring r=0,('+','.join(V)+'),dp;\nLIB "elim.lib";\nideal i='+',\n'.join(eqs)+';\nideal s=i;\n'
            +''.join(f'poly f{k}={f};\ns=sat(s,f{k});\n' for k,f in enumerate(facs))
            +'ideal m=minbase(s);\nprint("DIM"); print(dim(std(m))); print("GENS"); print(m); quit;\n')
    open(f'tmp/sg_in_{PID}.sing','w').write(script)
    try: r=subprocess.run([SG,'-q',f'tmp/sg_in_{PID}.sing'],capture_output=True,text=True,timeout=600)
    except subprocess.TimeoutExpired:
        out[str(ci)]={'timeout':True}; json.dump(out,open(fn,'w')); print(ci,'TIMEOUT',flush=True); continue
    lines=r.stdout.splitlines(); gens=[]; dim=None; mode=None
    for ln in lines:
        ln=ln.strip()
        if ln=='DIM': mode='d'; continue
        if ln=='GENS': mode='g'; continue
        if mode=='d' and ln: dim=int(ln); mode=None
        elif mode=='g' and ln and not ln.startswith('//'):
            g=ln.split('=',1)[1].strip().rstrip(',') if '=' in ln else ln.rstrip(',')
            gens.append(g)
    out[str(ci)]={'dim':dim,'gens':gens}
    json.dump(out,open(fn,'w')); print(ci,'dim',dim,'gens',len(gens),'%.0fs'%(time.time()-t0),flush=True)
print('done')
