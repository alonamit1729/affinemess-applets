"""H1: independent re-verification of the unrealizable verdicts.
d=4: types msolve called empty -> Singular: std(I + Rabinowitsch) == 1 ?
d=5,6,7: types with no real component -> if all primes are degenerate: Singular std(I+Rab)==1 AND msolve says empty;
         otherwise the components without real points carry a positive-definite-generator certificate (noreal2)."""
import sys, json, subprocess
d=int(sys.argv[1])
if d==4:
    src=open('certify54.py').read().split('results=[]')[0]; exec(src)
    C=json.load(open('cert54.json')); V=VARS
    todo=[(r['idx'],r['col']) for r in C if r['status']=='empty']
    allsys=lambda col: system(col)   # includes ineq
else:
    src=open('certify5d.py').read().split('import os')[0]; sys.argv=[sys.argv[0],str(d)]; exec(src)
    C=json.load(open(f'cert5{d}.json')); V=V+['t']
    todo=[]
    for r in C:
        if r.get('unresolved') or r['rigid']: continue
        if any(p.get('real') for p in r['primes']): continue
        if all(p['degenerate'] for p in r['primes']): todo.append((r['idx'],r['col']))
        else:
            assert all(p.get('empty_certificate') for p in r['primes'] if not p['degenerate']), r['idx']
    allsys=lambda col: system(col)+[expand(ineq)]
print(d,'types to re-verify:',len(todo),flush=True)
bad=[]
for i,(idx,col) in enumerate(todo):
    eqs=allsys(col)
    ring='ring r=0,('+','.join(V)+'),dp;'
    script=ring+'\nideal i='+',\n'.join(eqs)+';\nideal G=std(i);\nprint("RES"); print(size(G)==1 && deg(G[1])==0);\nquit;\n'
    open('rv.sing','w').write(script)
    r=subprocess.run([SG,'-q','rv.sing'],capture_output=True,text=True,timeout=600)
    out=[l.strip() for l in r.stdout.splitlines() if l.strip()]
    ok_sg = out[out.index('RES')+1]=='1'
    ok_ms=True
    if d!=4:
        o=run_msolve(eqs,V); ok_ms=o.replace(' ','').startswith('[-1]')
    if not(ok_sg and ok_ms): bad.append((idx,ok_sg,ok_ms))
    if i%20==0: print(i,'...',flush=True)
print('d',d,'re-verified',len(todo),'unrealizable types; disagreements:',bad)
json.dump({'d':d,'n':len(todo),'disagreements':bad},open(f'reverify_{d}.json','w'))
