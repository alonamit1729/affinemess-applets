"""Figures for the stratification summary page: the (5,4) incidence graph with configuration glyphs, the (5,5) valence
histogram, the degeneration-poset diagram.  Writes summary/figs.json with SVG strings."""
import json, math, collections, itertools, random
import networkx as nx, numpy as np
A=json.load(open('applet_data.json')); G4=json.load(open('strat_graph_4.json')); G5=json.load(open('strat_graph_5.json'))
def esc(s): return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
# ---------------- (5,4) graph ----------------
V={v['id']:v for v in G4['vertices']}
items3=A['cats']['5,3']['items']; pent=A['cats']['5,2']['items'][0]
def glyph_pts(v):
    it=pent if v['cell']==2 else items3[v['item']]
    P=np.array(it['draw'],float); P-=P.mean(0); P/=max(np.abs(P).max(),1e-9)*2.2
    return P
G=nx.Graph()
for i in V: G.add_node(i)
par=collections.Counter(); half=collections.defaultdict(list)
for e in G4['edges']:
    a,b=e['ends']
    if a[0]=='v' and b[0]=='v':
        G.add_edge(a[1],b[1],weight=G.get_edge_data(a[1],b[1],{}).get('weight',0)+1); par[tuple(sorted((a[1],b[1])))]+=1
    else:
        v=a if a[0]=='v' else b; o=b if a[0]=='v' else a; half[v[1]].append(o[1])
random.seed(3); pos=nx.kamada_kawai_layout(G,weight=None)
pos=nx.spring_layout(G,pos=pos,iterations=300,seed=3,k=0.32)
W,H=900,760; pad=70
xs=np.array([p[0] for p in pos.values()]); ys=np.array([p[1] for p in pos.values()])
def tr(p): return (pad+(p[0]-xs.min())/(xs.max()-xs.min())*(W-2*pad), pad+(p[1]-ys.min())/(ys.max()-ys.min())*(H-2*pad))
svg=[f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Incidence graph of the 25 one-parameter families of five-point four-distance sets">']
svg.append('<rect width="100%" height="100%" fill="var(--fig-bg)"/>')
# half-open ends as stubs
for i,ends in half.items():
    x,y=tr(pos[i]); c=collections.Counter(ends)
    k=0; n=len(ends)
    for kind,cnt in c.items():
        for j in range(cnt):
            x,y=tr(pos[i]); ang=2*math.pi*(k+0.5)/n+0.3; k+=1
            L=36 if kind=='degenerate' else 44
            x,y=x+19*math.cos(ang),y+19*math.sin(ang); x2,y2=x+(L-19)*math.cos(ang),y+(L-19)*math.sin(ang)
            if kind=='degenerate': svg.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="var(--end-deg)" stroke-width="1.1" stroke-dasharray="3 2"/>')
            else: svg.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="var(--end-inf)" stroke-width="1.4"/><circle cx="{x2:.1f}" cy="{y2:.1f}" r="2.2" fill="var(--end-inf)"/>')
# arcs: parallel arcs curved
for (a,b),cnt in par.items():
    (x1,y1),(x2,y2)=tr(pos[a]),tr(pos[b]); dx,dy=x2-x1,y2-y1; L=math.hypot(dx,dy) or 1; nx_,ny_=-dy/L,dx/L
    for j in range(cnt):
        off=(j-(cnt-1)/2)*14
        cx,cy=(x1+x2)/2+nx_*off,(y1+y2)/2+ny_*off
        svg.append(f'<path d="M{x1:.1f},{y1:.1f} Q{cx:.1f},{cy:.1f} {x2:.1f},{y2:.1f}" fill="none" stroke="var(--arc)" stroke-width="1.3" opacity="0.75"/>')
# glyph nodes
for i,v in V.items():
    x,y=tr(pos[i]); P=glyph_pts(v)*30
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="19" fill="var(--node-bg)" stroke="var(--node-edge)" stroke-width="1"/>')
    for a,b in itertools.combinations(range(5),2):
        svg.append(f'<line x1="{x+P[a,0]:.1f}" y1="{y-P[a,1]:.1f}" x2="{x+P[b,0]:.1f}" y2="{y-P[b,1]:.1f}" stroke="var(--glyph-line)" stroke-width="0.55"/>')
    for a in range(5): svg.append(f'<circle cx="{x+P[a,0]:.1f}" cy="{y-P[a,1]:.1f}" r="1.7" fill="var(--glyph-dot)"/>')
    lab='R₅' if v['cell']==2 else str(v['item']+1)
    svg.append(f'<text x="{x+21:.1f}" y="{y-14:.1f}" font-size="10.5" font-family="var(--mono)" fill="var(--fg-soft)">{lab}</text>')
svg.append('</svg>')
fig1=''.join(svg)
legend=[(('R₅' if v['cell']==2 else str(v['item']+1)), (v['label'] or '').replace(';',',')) for i,v in sorted(V.items(),key=lambda kv:(kv[1]['cell'],kv[1]['item'] if kv[1]['item'] is not None else -1))]
# ---------------- (5,5) valence histogram ----------------
V5={v['id']:v for v in G5['vertices']}
val=collections.Counter()
for e in G5['edges']:
    for x in e['ends']:
        if x[0]=='v': val[x[1]]+=1
def kind(v): return 'junction' if v['item'] is None else f"d{v['cell']}"
series={'d4':[],'d3':[],'junction':[],'d2':[]}
for i,v in V5.items(): series[kind(v)].append(val[i])
bins=[0,5,10,15,20,25,30,40,50,70,100,150,220]
W2,H2=720,300; padl,padb=46,46
rows=[]
for k,name,col in (('d4','rigid four-distance sets (468)','var(--bar-a)'),('d3','three-distance sets (34)','var(--bar-b)'),('junction','junction points on four-distance families (24)','var(--bar-c)'),('d2','the pentagon','var(--bar-d)')):
    h=[sum(1 for x in series[k] if bins[j]<=x<bins[j+1]) for j in range(len(bins)-1)]; rows.append((name,col,h))
nb=len(bins)-1; bw=(W2-padl-20)/nb; mx=max(max(r[2]) for r in rows)
s2=[f'<svg viewBox="0 0 {W2} {H2}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Valence distribution of the vertices of the five-distance incidence graph">','<rect width="100%" height="100%" fill="var(--fig-bg)"/>']
for g in (0,50,100,150,200):
    y=H2-padb-g/mx*(H2-padb-30); s2.append(f'<line x1="{padl}" y1="{y:.1f}" x2="{W2-20}" y2="{y:.1f}" stroke="var(--grid)" stroke-width="0.6"/><text x="{padl-6}" y="{y+3.5:.1f}" font-size="10" text-anchor="end" font-family="var(--mono)" fill="var(--fg-soft)">{g}</text>')
for j in range(nb):
    x0=padl+j*bw; y=H2-padb; 
    for name,col,h in rows:
        hh=h[j]/mx*(H2-padb-30)
        if hh>0: s2.append(f'<rect x="{x0+2:.1f}" y="{y-hh:.1f}" width="{bw-4:.1f}" height="{hh:.1f}" fill="{col}"/>'); y-=hh
    s2.append(f'<text x="{x0+bw/2:.1f}" y="{H2-padb+14}" font-size="10" text-anchor="middle" font-family="var(--mono)" fill="var(--fg-soft)">{bins[j]}–{bins[j+1]-1}</text>')
s2.append(f'<text x="{(W2+padl)/2:.0f}" y="{H2-8}" font-size="11" text-anchor="middle" font-family="var(--body)" fill="var(--fg-soft)">arcs at the vertex (valence)</text>')
s2.append(f'<text x="14" y="{(H2-padb)/2:.0f}" font-size="11" text-anchor="middle" font-family="var(--body)" fill="var(--fg-soft)" transform="rotate(-90 14 {(H2-padb)/2:.0f})">vertices</text>')
s2.append('</svg>'); fig2=''.join(s2); leg2=[(r[0],r[1]) for r in rows]
json.dump({'fig1':fig1,'legend1':legend,'fig2':fig2,'legend2':leg2},open('summary/figs.json','w'))
print('ok',len(fig1),len(fig2))
