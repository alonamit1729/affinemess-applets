import json,sys,os,random
n=int(sys.argv[1]); d=int(sys.argv[2]); idx=sys.argv[3]
sys.argv=['cellcert.py',str(n),str(d)]
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,'cellcert.py')).read(); exec(src[:src.index('results = json.load')])
R=json.load(open(os.path.join(HERE,f'cert_{n}_{d}.json'))); r=R[idx]; col=tuple(r['col'])
V,polys=system_gram(n,col); SV=V[:-1]; rab=polys[-1]
facs=[f's{a}-s{b}' for a,b in itertools.combinations(range(1,d),2)]+[f's{c}-1' for c in range(1,d)]+[f's{c}' for c in range(1,d)]
primes=decompose(SV,polys[:-1],facs)
for p in primes:
    if p['degenerate'] or p['dim']==0: continue
    print('prime dim',p['dim'],'gens',p['gens'])
    tot=0; psdc=0; ranks=set()
    for _ in range(30):
        lins=['+'.join(f'({random.randint(-7,7)})*{v}' for v in SV)+f'-({random.randint(-9,9)})' for _ in range(p['dim'])]
        s2=parse_real(msolve_real_raw(V,p['gens']+lins+[rab],precision=200))
        if isinstance(s2,list):
            for sol in s2:
                tot+=1; vals=sol_vals(V,sol); psd,rank,P=gram_sets(col,vals)
                if psd: psdc+=1; ranks.add(rank)
                elif tot<=4: print('  nonpsd point vals',{k:float(v) for k,v in vals.items()})
    print('  real points on slices',tot,'psd',psdc,ranks)
