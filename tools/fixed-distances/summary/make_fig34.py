import json, math, itertools, collections, numpy as np, networkx as nx
A=json.load(open('applet_data.json')); G=json.load(open('strat_graph_4_3.json'))
V={v['id']:v for v in G['vertices']}
items2=A['cats']['4,2']['items']
def glyph(v):
    if v['item'] is None:   # fold: collinear equally spaced
        P=np.array([[0,0],[1,0],[2,0],[3,0]],float)
    else: P=np.array(items2[v['item']]['draw'],float)
    P-=P.mean(0); P/=max(np.abs(P).max(),1e-9)*2.2; return P
Gx=nx.Graph(); par=collections.Counter(); half=collections.defaultdict(list)
for i in V: Gx.add_node(i)
for e in G['edges']:
    a,b=e['ends']
    if a[0]=='v' and b[0]=='v': Gx.add_edge(a[1],b[1]); par[tuple(sorted((a[1],b[1])))]+=1
    else:
        v=a if a[0]=='v' else b; o=b if a[0]=='v' else a; half[v[1]].append('degenerate' if o[1] in ('degenerate','inf') else o[1])
pos=nx.circular_layout(Gx); pos=nx.spring_layout(Gx,pos=pos,iterations=200,seed=5,k=0.9)
W,H=700,520; pad=80
xs=np.array([p[0] for p in pos.values()]); ys=np.array([p[1] for p in pos.values()])
def tr(p): return (pad+(p[0]-xs.min())/(xs.max()-xs.min()+1e-9)*(W-2*pad), pad+(p[1]-ys.min())/(ys.max()-ys.min()+1e-9)*(H-2*pad))
s=[f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Incidence graph of the nine one-parameter families of four-point three-distance sets">','<rect width="100%" height="100%" fill="var(--fig-bg)"/>']
for i,ends in half.items():
    x0,y0=tr(pos[i]); n=len(ends)
    for k,kind in enumerate(ends):
        ang=2*math.pi*(k+0.5)/n+0.9; x=x0+24*math.cos(ang); y=y0+24*math.sin(ang); x2=x0+46*math.cos(ang); y2=y0+46*math.sin(ang)
        s.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="var(--end-deg)" stroke-width="1.3" stroke-dasharray="3 2"/>')
# per-type colours for the 9 families
cols=['#1f6f78','#b5542b','#8a3fa0','#d9a441','#3a7d44','#2f5fa8','#a83a5c','#6b7a2a','#555f68']
for e in G['edges']:
    a,b=e['ends']
    if not (a[0]=='v' and b[0]=='v'): continue
    (x1,y1),(x2,y2)=tr(pos[a[1]]),tr(pos[b[1]])
    key=tuple(sorted((a[1],b[1]))); cnt=par[key]; j=e.setdefault('_j',0)
    # distribute parallel arcs
    done=sum(1 for f in G['edges'] if f is not e and f.get('_drawn') and tuple(sorted((f['ends'][0][1] if f['ends'][0][0]=='v' else -1,f['ends'][1][1] if f['ends'][1][0]=='v' else -1)))==key)
    e['_drawn']=True; off=(done-(cnt-1)/2)*18
    dx,dy=x2-x1,y2-y1; L=math.hypot(dx,dy) or 1; nx_,ny_=-dy/L,dx/L; cx,cy=(x1+x2)/2+nx_*off,(y1+y2)/2+ny_*off
    s.append(f'<path d="M{x1:.1f},{y1:.1f} Q{cx:.1f},{cy:.1f} {x2:.1f},{y2:.1f}" fill="none" stroke="{cols[e["type"]]}" stroke-width="1.8" opacity="0.85"/>')
for i,v in V.items():
    x,y=tr(pos[i]); P=glyph(v)*34
    s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="24" fill="var(--node-bg)" stroke="var(--node-edge)" stroke-width="1"/>')
    for a,b in itertools.combinations(range(4),2): s.append(f'<line x1="{x+P[a,0]:.1f}" y1="{y-P[a,1]:.1f}" x2="{x+P[b,0]:.1f}" y2="{y-P[b,1]:.1f}" stroke="var(--glyph-line)" stroke-width="0.7"/>')
    for a in range(4): s.append(f'<circle cx="{x+P[a,0]:.1f}" cy="{y-P[a,1]:.1f}" r="2.1" fill="var(--glyph-dot)"/>')
    lab='fold' if v['item'] is None else str(v['item']+1)
    s.append(f'<text x="{x+26:.1f}" y="{y-18:.1f}" font-size="11" font-family="var(--mono)" fill="var(--fg-soft)">{lab}</text>')
s.append('</svg>'); fig4=''.join(s)
E4=list(itertools.combinations(range(4),2))
perms=[[E4.index(tuple(sorted((p[a],p[b])))) for a,b in E4] for p in itertools.permutations(range(4))]
def ccol(rgs):
    best=None
    for pe in perms:
        c=[rgs[pe[i]] for i in range(6)]; mp={}; out=[]
        for x in c:
            if x not in mp: mp[x]=len(mp)
            out.append(mp[x])
        t=tuple(out)
        if best is None or t<best: best=t
    return best
R=json.load(open('cert43.json')); bycol={ccol(tuple(r['col'])):r['idx'] for r in R}
fam=A['cats']['4,3']['items']
leg=[]
for k,it in enumerate(fam):
    t=bycol[ccol(tuple(it['cls']))]; leg.append((cols[t],f"#{k+1} "+it['label']))

leg_v=[(str(k+1),it['label']) for k,it in enumerate(items2)]+[('fold','four equally spaced collinear points — the family passes through its own mirror image')]
# ---- n = 3: the shape region of triangles ----
W3,H3=560,440
def tri(cx,cy,pts,sc=16):
    P=np.array(pts,float); P-=P.mean(0); P/=max(np.abs(P).max(),1e-9); out=[]
    for a,b in ((0,1),(1,2),(0,2)): out.append(f'<line x1="{cx+sc*P[a,0]:.1f}" y1="{cy-sc*P[a,1]:.1f}" x2="{cx+sc*P[b,0]:.1f}" y2="{cy-sc*P[b,1]:.1f}" stroke="var(--glyph-line)" stroke-width="0.9"/>')
    for a in range(3): out.append(f'<circle cx="{cx+sc*P[a,0]:.1f}" cy="{cy-sc*P[a,1]:.1f}" r="2" fill="var(--glyph-dot)"/>')
    return ''.join(out)
E=(280,70); C=(100,380); X=(460,380)
t=[f'<svg viewBox="0 0 {W3} {H3}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="The shape space of triangles with its strata">','<rect width="100%" height="100%" fill="var(--fig-bg)"/>']
t.append(f'<path d="M{E[0]},{E[1]} Q{170},{230} {C[0]},{C[1]} Q{280},{330} {X[0]},{X[1]} Q{370},{230} {E[0]},{E[1]} Z" fill="var(--region)" stroke="none"/>')
t.append(f'<path d="M{E[0]},{E[1]} Q{170},{230} {C[0]},{C[1]}" fill="none" stroke="var(--arc)" stroke-width="2.2"/>')
t.append(f'<path d="M{E[0]},{E[1]} Q{390},{230} {X[0]},{X[1]}" fill="none" stroke="var(--arc)" stroke-width="2.2"/>')
t.append(f'<path d="M{C[0]},{C[1]} Q{280},{330} {X[0]},{X[1]}" fill="none" stroke="var(--end-inf)" stroke-width="2.2"/>')
# glyphs along the isosceles arcs and the collinear edge
def isoc(theta):  # apex angle theta, legs 1
    return [(-math.sin(theta/2),0),(math.sin(theta/2),0),(0,math.cos(theta/2))]
for frac,theta in ((0.0,math.pi/3),(0.35,math.radians(100)),(0.7,math.radians(150)),(1.0,math.pi-0.001)):
    # quadratic bezier E->C
    p=lambda f:((1-f)**2*E[0]+2*(1-f)*f*170+f*f*C[0],(1-f)**2*E[1]+2*(1-f)*f*230+f*f*C[1])
    cx,cy=p(frac); t.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="19" fill="var(--node-bg)" stroke="var(--node-edge)"/>'+tri(cx,cy,isoc(theta)))
for frac,theta in ((0.35,math.radians(40)),(0.7,math.radians(18)),(1.0,math.radians(1))):
    p=lambda f:((1-f)**2*E[0]+2*(1-f)*f*390+f*f*X[0],(1-f)**2*E[1]+2*(1-f)*f*230+f*f*X[1])
    cx,cy=p(frac); t.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="19" fill="var(--node-bg)" stroke="var(--node-edge)"/>'+tri(cx,cy,isoc(theta)))
for frac,a in ((0.33,0.6),(0.66,0.25)):
    p=lambda f:((1-f)**2*C[0]+2*(1-f)*f*280+f*f*X[0],(1-f)**2*C[1]+2*(1-f)*f*330+f*f*X[1])
    cx,cy=p(frac); t.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="19" fill="var(--node-bg)" stroke="var(--node-edge)"/>'+tri(cx,cy,[(0,0),(a,0),(1,0)]))
t.append(f'<circle cx="{280}" cy="{225}" r="19" fill="var(--node-bg)" stroke="var(--node-edge)"/>'+tri(280,225,[(0,0),(1,0),(0.3,0.75)]))
t.append(f'<text x="{E[0]}" y="{E[1]-28}" text-anchor="middle" font-size="12" font-family="var(--body)" fill="var(--fg)">equilateral (1 distance)</text>')
t.append(f'<text x="{C[0]+10}" y="{C[1]+34}" text-anchor="middle" font-size="12" font-family="var(--body)" fill="var(--fg)">collinear, apex at the midpoint</text><text x="{C[0]+10}" y="{C[1]+48}" text-anchor="middle" font-size="11" font-family="var(--body)" fill="var(--fg-soft)">(fold of the isosceles family)</text>')
t.append(f'<text x="{X[0]}" y="{X[1]+34}" text-anchor="middle" font-size="12" font-family="var(--body)" fill="var(--fg)">two points coincide</text><text x="{X[0]}" y="{X[1]+48}" text-anchor="middle" font-size="11" font-family="var(--body)" fill="var(--fg-soft)">(not in the space)</text>')
t.append(f'<text x="196" y="118" text-anchor="end" font-size="11.5" font-family="var(--body)" fill="var(--arc)">isosceles, apex angle 60° → 180°</text>')
t.append(f'<text x="364" y="118" text-anchor="start" font-size="11.5" font-family="var(--body)" fill="var(--arc)">isosceles, apex angle 60° → 0°</text>')
t.append(f'<text x="280" y="285" text-anchor="middle" font-size="11.5" font-family="var(--body)" fill="var(--fg-soft)">scalene (3 distances): the interior</text>')
t.append(f'<text x="280" y="322" text-anchor="middle" font-size="11.5" font-family="var(--body)" fill="var(--end-inf)">collinear scalene triples: a mirror edge</text>')
t.append('</svg>'); fig3=''.join(t)
json.dump({'fig4':fig4,'legend4':leg,'legend4v':leg_v,'fig3':fig3},open('summary/figs34.json','w')); print('ok')
