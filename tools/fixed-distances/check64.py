"""(6,4) families: for each positive-dimensional non-degenerate prime of the four flexible types, the singular locus
(Singular slocus) and its real points with the inequation (msolve): no real singular point => the real locus is a smooth
curve, so there are no real-isolated (hidden rigid) points and the family has real dimension 1 everywhere.
Also: Jacobian rank at the certified slice points.  Usage: python3 check64.py  (reads cert64L.json)"""
import sys, json, subprocess, re, itertools
sys.argv=['certify_nd.py','6','4','L']
src=open('certify_nd.py').read().split('import os\nresults=')[0]
exec(src)   # classes, system, decompose, run_msolve, parse_real, expand, ineq, V, XV, SV, facs, SG, MS, PID
import sympy as sp
R=json.load(open('cert64L.json'))
flex=[r for r in R if any(q['real'] and q['dim']>0 for q in r['primes'])]
out=[]
for r in flex:
    col=classes[r['idx']]; eqs=system(col); primes=decompose(eqs)
    for p in primes:
        if p['degenerate'] or p['dim']==0: continue
        gens=p['gens']
        # singular locus via Singular
        ring='ring r=0,('+','.join(V)+'),dp;'
        script=ring+'\nLIB "sing.lib";\nideal P='+',\n'.join(gens)+';\nideal S=slocus(P);\nS=std(S);\nprint(dim(S));\nprint(S);\nquit;\n'
        open(f'tmp/sl_{PID}.sing','w').write(script)
        rr=subprocess.run([SG,'-q',f'tmp/sl_{PID}.sing'],capture_output=True,text=True,timeout=600)
        lines=[l.strip() for l in rr.stdout.splitlines() if l.strip() and not l.startswith('//')]
        sdim=int(lines[0]); sgens=[l.split('=',1)[1].strip().rstrip(',') if '=' in l else l.rstrip(',') for l in lines[1:]]
        # real points of the singular locus satisfying the inequation
        st,md=parse_real(run_msolve(sgens+[expand(ineq)],V+['t'])) if sdim>=0 else ('empty',[])
        # Jacobian rank at the slice point
        pt=p.get('point') if 'point' in p else None
        J=sp.Matrix([[sp.diff(sp.sympify(g.replace('^','**'),locals=SYM),SYM[v]) for v in V] for g in gens])
        rank=None
        # evaluate at the certified real point from cert64L.json for this prime (match by dim; take q['point'])
        q=[q for q in r['primes'] if q['real'] and q['dim']==p['dim']]
        if q and 'point' in q[0]:
            P=q[0]['point']; vals={ 'x1':P[1][0] }
            for k in range(2,6): vals[f'x{k}']=P[k][0]; vals[f'y{k}']=P[k][1]
            dd=[[ (P[i][0]-P[j][0])**2+(P[i][1]-P[j][1])**2 for j in range(6)] for i in range(6)]
            E6=list(itertools.combinations(range(6),2))
            for c in range(1,4):
                e=[e for e,(i,j) in enumerate(E6) if col[e]==c][0]; vals[f's{c}']=dd[E6[e][0]][E6[e][1]]
            Jn=J.subs({SYM[v]:vals[v] for v in V}); rank=Jn.rank(iszerofunc=lambda x: abs(x)<1e-7)
        rec={'type':r['idx'],'dim':p['dim'],'slocus_dim':sdim,'real_singular_points':len(md) if st=='zerodim' else st,'jacobian_rank_at_point':rank,'nvars':len(V)}
        out.append(rec); print(rec,flush=True)
json.dump(out,open('check64.json','w'))
