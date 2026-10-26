"""Rebuild the Quora reader pages (quora/<slug>.html), quora/index.json and quora/search.json.

Inputs are PRIVATE and stay off the repo:
  QUORA_EXPORT  Alon's Quora data export folder ("2026 Download"; on his Mac under Alon/Quora/)
  QUORA_WORK    a folder with match.json, answers_raw.json and links1.txt (the qr.ae link list
                matched to the export; regenerate by resolving the qr.ae links, see project-knowledge.md)

usage: QUORA_EXPORT=... QUORA_WORK=... python3 src/build_quora.py
Then run python3 src/build_site.py so the Quora rows on the site are regenerated.
"""
import json, re, os, sys, shutil, html as H
from urllib.parse import urlparse, parse_qs, unquote
from datetime import datetime
from bs4 import BeautifulSoup, NavigableString
REPO=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q=os.environ.get('QUORA_WORK','/home/claude/quora')
UP=os.environ.get('QUORA_EXPORT','/mnt/user-data/uploads/2026 Download')
m=json.load(open(Q+'/match.json')); A=json.load(open(Q+'/answers_raw.json'))
titles=dict(l.rstrip('\n').split('|',1) for l in open(Q+'/links1.txt') if '|' in l)
os.makedirs(REPO+'/quora/img',exist_ok=True)
def ext(b):
    if b[:8]==b'\x89PNG\r\n\x1a\n': return 'png'
    if b[:3]==b'\xff\xd8\xff': return 'jpg'
    if b[:3]==b'GIF': return 'gif'
    if b[:4]==b'RIFF' and b[8:12]==b'WEBP': return 'webp'
    return 'png'
def slugify(t):
    s=re.sub(r'[^a-z0-9]+','-',t.lower()).strip('-')
    return s[:60].rstrip('-') or 'answer'
def parse_date(d):
    return datetime.strptime(re.sub(r'\s+\d+:\d+ [AP]M.*','',d),'%b %d, %Y').strftime('%Y-%m-%d')
CSS=open(REPO+'/hcn-explorer.html').read().split('/* am-kit v1 BEGIN */')[1].split('/* am-kit v1 END */')[0]
PAGE='''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>%(t)s · affine·mess</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght,SOFT@0,9..144,300..600,0..100;1,9..144,300..500,0..100&family=JetBrains+Mono:wght@400;500&display=swap">
<link rel="stylesheet" href="katex/katex.min.css">
<style>
/* am-kit v1 BEGIN */%(css)s/* am-kit v1 END */
html,body{height:auto}
.q-wrap{max-width:720px;margin:0 auto;padding:28px 22px 60px}
.q-wrap h1{font-style:italic;font-weight:400;font-size:28px;line-height:1.2;margin:6px 0 8px;font-variation-settings:"opsz" 144,"SOFT" 80;letter-spacing:-.01em}
.q-meta{font-family:var(--mono);font-size:10.5px;letter-spacing:.12em;text-transform:uppercase;color:var(--ink-faint);margin-bottom:22px}
.q-body{font-size:17.5px;line-height:1.65}
.q-body p{margin:0 0 1em}
.q-body a{color:var(--accent-deep)}
.q-body img{max-width:100%%;height:auto;display:block;margin:1em auto;background:#fff;border-radius:2px}
.q-body blockquote{margin:1em 0;padding:0 0 0 16px;border-left:3px solid var(--line);color:var(--ink-soft)}
.q-body h3{font-size:20px;font-style:italic;font-weight:400;margin:1.4em 0 .5em}
.q-body pre,.q-body code{font-family:var(--mono);font-size:14px;background:var(--paper-card);border-radius:3px}
.q-body pre{padding:10px 12px;overflow-x:auto}
.q-body .tex{white-space:nowrap}
.q-body .katex-display,.q-body .tex-block{overflow-x:auto;overflow-y:hidden}
.q-body p.citation{font-size:14px;color:var(--ink-soft);margin:0 0 .4em}
.q-foot{margin-top:34px;padding-top:14px;border-top:1px solid var(--line-soft);font-size:15px;color:var(--ink-soft)}
</style></head><body>
<article class="q-wrap">
<div class="q-meta">Alon Amit on Quora &middot; %(d)s</div>
<h1>%(th)s</h1>
<div class="q-body">%(body)s</div>
<div class="q-foot">Originally posted on Quora: <a href="%(url)s" target="_blank" rel="noopener">%(url)s</a></div>
</article>
<script src="katex/katex.min.js"></script>
<script>
document.querySelectorAll('.tex').forEach(function(e){try{katex.render(e.getAttribute('data-tex'),e,{throwOnError:false,displayMode:e.classList.contains('tex-block')||/\\\\tag/.test(e.getAttribute('data-tex'))});}catch(x){}});
if(window.top===window)document.documentElement.classList.add('am-standalone');
</script></body></html>'''
index=[]; used=set(); imgmap={}; missing=[]
for code,idxs in m.items():
    a=A[idxs[0]]; d=a['file'].split('/')[0]
    title=a['q'] if a['q'] else titles[code]
    raw=re.sub(r'\s+',' ',title).strip()
    parts=re.split(r'\[math\](.*?)\[/math\]',raw,flags=re.S)
    th=''.join(H.escape(p,quote=False) if k%2==0 else '<span class="tex" data-tex="%s"></span>'%H.escape(p.replace('\\displaystyle','').replace('_\\limits','\\limits_'),quote=True) for k,p in enumerate(parts))
    title=re.sub(r'\[/?math\]','',raw).replace('\\displaystyle','')
    s=slugify(title); base=s; n=2
    while s in used: s=f'{base}-{n}'; n+=1
    used.add(s)
    soup=BeautifulSoup(a['html'],'html.parser')
    # math
    for sp in soup.select('span.render_latex'):
        t=sp.get_text()
        mm=re.fullmatch(r'\s*\[math\](.*)\[/math\]\s*',t,re.S)
        tex=mm.group(1) if mm else t
        tex=tex.replace('_\\limits','\\limits_').replace('^\\limits','\\limits^')
        tex=re.sub(r'\\(begin|end)\{(align\*?|eqnarray\*?)\}',r'\\\1{aligned}',tex)
        new=soup.new_tag('span'); new['class']='tex'; new['data-tex']=tex; new.string=tex
        sp.replace_with(new)
    for el in list(soup.find_all(string=re.compile(r'\[math\]'))):
        pass
    # links
    for aa in soup.find_all('a'):
        h=aa.get('href','')
        if 'quora.com/_/redirect' in h:
            u=parse_qs(urlparse(H.unescape(h)).query).get('url',[''])[0]
            h=unquote(u) or h
        for k in list(aa.attrs):
            if k not in('href',): del aa.attrs[k]
        aa['href']=h
        if h.startswith('http'): aa['target']='_blank'; aa['rel']='noopener'
    for sp in soup.select('span.qlink_container'): sp.unwrap()
    # images
    for g in soup.select('div.ui_qtext_gif_noclick_wrapper'):
        pass
    for img in soup.find_all('img'):
        h=img['src'].split('/')[-1]
        src=f'{UP}/{d}/images/{h}'
        if not os.path.exists(src): missing.append(src); img.decompose(); continue
        b=open(src,'rb').read(16); e=ext(b); fn=f'{h[5:13]}{h[-6:]}.{e}'
        shutil.copy(src,f'{REPO}/quora/img/{fn}')
        for k in list(img.attrs): del img.attrs[k]
        img['src']=f'img/{fn}'; img['alt']=''; img['loading']='lazy'
    for w in soup.select('div.ui_qtext_gif_noclick_wrapper'):
        w.unwrap()
    for w in soup.select('div.qtext_image_wrapper,div.ui_qtext_image_outer'): w.unwrap()
    for t in soup.find_all(True):
        if t.name in('a','img','span'): 
            if t.name=='span' and 'tex' in (t.get('class') or []): pass
            continue
        for k in list(t.attrs):
            if not(t.name=='p' and t.get('class')==['citation'] and k in('class','id')): del t.attrs[k]
    for t in soup.find_all('a'):
        if t.get('href','').startswith('#cite'): pass
    # unwrap stray divs/spans
    for t in soup.find_all(['div','span']):
        if t.name=='span' and 'tex' in (t.get('class') or []): continue
        t.unwrap()
    # bare text/br chunks into paragraphs when no <p>
    if not soup.find('p'):
        parts=re.split(r'(?:<br\s*/?>\s*){2,}',str(soup))
        body=''.join(f'<p>{p.strip()}</p>' for p in parts if p.strip())
    else: body=str(soup)
    # plain text for excerpt
    pt=BeautifulSoup(body,'html.parser')
    for sp in pt.select('.tex'):
        tx=sp['data-tex'].strip()
        sp.replace_with(' … ' if ('\\' in tx or '{' in tx or len(tx)>14) else ' '+tx+' ')
    txt=re.sub(r'\s+',' ',pt.get_text(' ')).strip()
    txt=re.sub(r' ([,.;:!?)])',r'\1',txt)
    ex=txt[:200].rsplit(' ',1)[0]+('…' if len(txt)>200 else '')
    date=parse_date(a['date']); url=f'https://qr.ae/{code}'
    open(f'{REPO}/quora/{s}.html','w').write(PAGE%dict(t=H.escape(title),th=th,css=CSS,d=date,body=body,url=url))
    index.append(dict(slug=s,title=title,th=th,date=date,excerpt=ex,url=url,words=len(txt.split())))
gen_order=list(index)  # search.json keeps generation order
index.sort(key=lambda x:x['date'],reverse=True)
json.dump(index,open(REPO+'/quora/index.json','w'),indent=1,ensure_ascii=False)
# full-text search index: slug -> title + plain text of the answer (math kept as its TeX source)
def _plain(h):
    s=BeautifulSoup(h,'html.parser')
    for t in s.select('.tex'): t.replace_with(' '+t['data-tex']+' ')
    return s.get_text(' ')
search={}
for it in gen_order:
    b=BeautifulSoup(open(f"{REPO}/quora/{it['slug']}.html").read(),'html.parser').select_one('.q-body')
    search[it['slug']]=it['title']+' — '+re.sub(r'\s+',' ',_plain(str(b))).strip()
json.dump(search,open(REPO+'/quora/search.json','w'),ensure_ascii=False,separators=(',',':'))
print(len(index),'missing',missing)
