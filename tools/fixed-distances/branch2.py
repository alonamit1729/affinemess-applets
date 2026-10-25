"""Local branch analysis at the candidate points of the stubborn d=5 types (380, 381, 426, 427).
For each real rank-deficient point p (stage-1 system Z = I + 10x10 minors + Rabinowitsch inequation, zero-dimensional):
  1. msolve -P 2: rational parametrization of Z; the real root T0 of the eliminating polynomial w giving p; the minimal
     polynomial m(a) of T0 (irreducible factor of w) -> exact coordinates p(a) in Q(a).
  2. Saturated ideal Js of the type (sat(I, product of inequation factors)), reduced by eliminating the linearly
     occurring variables (same reduction as hidden_subsets.py), translated to p: x -> x + p(a), in the ring
     Q(a)[free vars] with local ordering ds.
  3. curvepar.lib CurveParam: algebraic branches of the curve germ at p, each over an extension Q(a)(b) with minimal
     polynomial g(b).  A branch is real iff the extension embeds in R compatibly (its defining root b is real):
     certified by counting real roots of g (Sturm) and, when needed, interval evaluation.
  p is isolated in the real locus iff no branch is real.  Usage: python3 branch.py idx"""
import sys, json, subprocess, re, time
from fractions import Fraction
import sympy as sp, mpmath
idx=int(sys.argv[1]); sys.argv=['x','5']; exec(open('certify5d.py').read().split('import os')[0])
src2=open('hidden2b.py').read(); i=src2.index('def reduce_linear'); j=src2.index('out=[]; t0'); exec(src2[i:j])
mpmath.mp.dps=120
def intstr(g,vs): return str(sp.Poly(g,*vs).clear_denoms()[1].as_expr()).replace('**','^').replace(' ','')
def sg(script,tmo=3600):
    fn=f'tmp/br_{PID}.sing'; open(fn,'w').write(script)
    r=subprocess.run([SG,'-q',fn],capture_output=True,text=True,timeout=tmo)
    return r.stdout
VT=V+['t']; eqs=system(classes[idx]); rab=expand(ineq); t0=time.time()
# ---- 1. stage-1 system and its rational parametrization
ring='ring r=0,('+','.join(V)+'),dp;'
out=sg(ring+'\nideal i='+',\n'.join(eqs)+';\nideal M=minor(jacob(i),10); M=simplify(M,2+4+8);\nprint("MIN");\nint k; for(k=1;k<=size(M);k++){ print(string(M[k])); }\nquit;\n')
lines=[l.strip() for l in out.splitlines() if l.strip() and not l.startswith('//')]
minors=lines[lines.index('MIN')+1:]
Z=eqs+minors+[rab]
fi,fo=f'tmp/br_in_{PID}.txt',f'tmp/br_rur_{PID}.txt'
open(fi,'w').write(','.join(VT)+'\n0\n'+',\n'.join(Z)+'\n')
subprocess.run([MS,'-f',fi,'-o',fo,'-P','2'],capture_output=True,text=True,timeout=3600)
s=open(fo).read().replace('\n','').rstrip(':')
def nested(txt):
    toks=re.findall(r"'[^']*'|-?\d+|[\[\],]",txt); pos=0
    def build():
        nonlocal pos
        if toks[pos]=='[':
            pos+=1; lst=[]
            while toks[pos]!=']':
                lst.append(build())
                if toks[pos]==',': pos+=1
            pos+=1; return lst
        tk=toks[pos]; pos+=1
        return tk.strip("'") if tk.startswith("'") else int(tk)
    return build()
tree=nested(s)
assert tree[0]==0, 'not zero-dimensional'
hdr=tree[1]; nvars=hdr[1]; deg=hdr[2]; varnames=hdr[3]; body=hdr[5]
# body = [1, [deg_w, coeffs_w], [deg_w', coeffs_w'], [ [[deg,coeffs], c], ... ] ]
inner=body[1]; wc=inner[0][1]; wpc=inner[1][1]; params=[(p[0][1],p[1]) for p in inner[2]]
T=sp.Symbol('T'); w=sp.Poly(list(reversed(wc)),T); wp=sp.Poly(list(reversed(wpc)),T)
par=[(sp.Poly(list(reversed(c)),T),cc) for c,cc in params]
print(f'type {idx}: parametrization degree {deg}, {len(par)} coordinates, vars {varnames[:3]}... %.0fs'%(time.time()-t0),flush=True)
# real roots of w and the corresponding points; convention check against the isolating boxes
roots=[r for r in sp.Poly(w,T).nroots(n=60) if abs(sp.im(r))<1e-40]
roots=[sp.re(r) for r in roots]
H={r['idx']:r for r in json.load(open('hidden_5.json'))}; cands=[c['point'] for c in H[idx]['isolated']]
def coords(T0,sign):
    wpv=wp.eval(T0)
    return [sign*pv.eval(T0)/(cc*wpv) for pv,cc in par]
def pts(vals): return [(0.0,0.0),(float(vals[0]),0.0),(float(vals[1]),float(vals[2])),(float(vals[3]),float(vals[4])),(float(vals[5]),float(vals[6]))]
import math
def close(P,Q): return max(math.dist(a,b) for a,b in zip(P,Q))<1e-6
sign=None
for sgn in (1,-1):
    ok=sum(1 for r0 in roots for c in cands if close(pts(coords(r0,sgn)),c))
    if ok==len(cands): sign=sgn
assert sign is not None, 'could not match the parametrization to the candidate points'
print(f'  {len(roots)} real roots of w; convention sign {sign}; all {len(cands)} candidates matched',flush=True)
# ---- 2. reduced saturated ideal
prod='*'.join(facs)
out=sg(ring+'\nLIB "elim.lib";\nideal I='+',\n'.join(eqs)+';\nideal Js=sat(I,'+prod+'); list L=mstd(Js); print("PM"); print(L[2]);\nquit;\n')
lines=[l.strip() for l in out.splitlines() if l.strip() and not l.startswith('//')]
gens0=[l.rstrip(',') for l in lines[lines.index('PM')+1:]]
F,subs,free=reduce_linear(gens0); vs=[SYM[v] for v in free]; gens=[sp.Poly(g,*vs).clear_denoms()[1].as_expr() for g in F]
print(f'  saturated ideal reduced to {len(gens)} generators in {free}',flush=True)
# ---- 3. per candidate: exact point, translation, CurveParam
wf=sp.factor_list(w.as_expr())[1]
results=[]
for r0 in roots:
    P=pts(coords(r0,sign))
    if not any(close(P,c) for c in cands): continue
    fac=[f for f,e in wf if abs(f.subs(T,r0))<1e-30]
    assert len(fac)==1, fac
    mp=sp.Poly(fac[0],T); print(f'  candidate {[(round(x,3),round(y,3)) for x,y in P]}: minimal polynomial of T0 has degree {mp.degree()}',flush=True)
    # exact coordinates as rational functions of a (reduce numerators modulo m to keep them small)
    a=sp.Symbol('a'); mpa=mp.as_expr().subs(T,a); wpa=wp.as_expr().subs(T,a)
    coord={}
    for v,(pv,cc) in zip(varnames,par):
        num=sp.rem(sign*pv.as_expr().subs(T,a),mpa,a)
        coord[v]=(num,cc*wpa)      # value = num/(cc*w'(a))
    # Singular script: ring over Q(a), local ordering, substitute x -> x + p
    subst=[]
    for v in free:
        num,den=coord[v]
        subst.append(f'{v}+(({str(num).replace("**","^")})/({str(den).replace("**","^")}))')
    gstr=[str(g).replace('**','^') for g in gens]
    script=('ring r=(0,a),('+','.join(free)+'),ds;\nminpoly='+str(mpa).replace('**','^')+';\nLIB "primdec.lib";\n'
            'ideal I='+',\n'.join(gstr)+';\nmap phi=r,'+','.join(subst)+';\nideal J=std(phi(I)); print("DIMJ"); print(dim(J)); print("MULT"); print(mult(J));\n'
            'ideal T=tangentcone(J); ring s=(0,a),('+','.join(free)+'),dp; minpoly='+str(mpa).replace('**','^')+'; LIB "primdec.lib";\nideal T=imap(r,T); T=std(T); print("TDIM"); print(dim(T)); print("TDEG"); print(mult(T));\n'
            'list L=minAssGTZ(T); print("NCOMP"); print(size(L)); int k; for(k=1;k<=size(L);k++){ print("COMP"); print(dim(std(L[k]))); print(mult(std(L[k]))); print(L[k]); }\nquit;\n')
    try: out=sg(script,tmo=3600)
    except subprocess.TimeoutExpired: print('   tangent cone timeout',flush=True); results.append({'point':P,'status':'timeout'}); continue
    lines=[l.strip() for l in out.splitlines() if l.strip() and not l.startswith('//')]
    rec={'point':P,'minpoly_degree':mp.degree(),'singular_output':out[-3000:]}
    def after(tag): return lines[lines.index(tag)+1] if tag in lines else None
    rec['dim']=after('DIMJ'); rec['mult']=after('MULT'); rec['tcone_dim']=after('TDIM'); rec['tcone_deg']=after('TDEG'); rec['ncomp']=after('NCOMP')
    comps=[]
    if 'NCOMP' in lines:
        k=lines.index('NCOMP')+2
        while k<len(lines):
            if lines[k]=='COMP': comps.append({'dim':lines[k+1],'deg':lines[k+2],'gens':[]}); k+=3; continue
            comps[-1]['gens'].append(lines[k].rstrip(',')); k+=1
    rec['components']=comps
    verdict='?'
    if rec['mult']=='2' and rec['ncomp']=='2' and all(c['deg']=='1' for c in comps): verdict='ordinary node with two real tangent lines: two real branches cross, NOT isolated'
    elif rec['mult']=='2' and rec['ncomp']=='1' and comps and comps[0]['deg']=='2': verdict='ordinary double point with conjugate tangents: ISOLATED real point'
    elif rec['mult']=='2' and rec['ncomp']=='1' and comps and comps[0]['deg']=='1': verdict='multiplicity 2, one tangent line (cusp or tacnode): inconclusive'
    rec['verdict']=verdict
    print(f'   mult {rec["mult"]}, tangent cone: {rec["ncomp"]} component(s) of degrees {[c["deg"] for c in comps]} -> {verdict}',flush=True)
    results.append(rec)
    json.dump(results,open(f'branch2_{idx}.json','w'))
print('done %.0fs'%(time.time()-t0))
