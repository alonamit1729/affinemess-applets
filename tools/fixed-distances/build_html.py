"""Inject data/applet_data.json into the template.  Cells with more than LIMIT items are *slimmed*: their items are written
to data/catalogs/{n}_{d}.json (in their original order, so card numbers are stable; fetched by the page on demand, path
relative to the site root) and removed from applet_data.json, which keeps only the counts and a `remote` pointer.
`datautil.load()` merges them back for scripts that need every item.  Idempotent: an already slim applet_data.json
is left as it is, and its catalogs are not rewritten."""
import json, os, re
HERE=os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(HERE,'data'))            # the data and its catalogs/ live in data/
LIMIT=1000
data=json.load(open('applet_data.json'))
os.makedirs('catalogs',exist_ok=True)
slimmed=False
for key,cat in data['cats'].items():
    if cat['items'] and len(cat['items'])>LIMIT:
        n,d=key.split(','); fn=f'catalogs/{n}_{d}.json'
        new=json.dumps(cat['items'],ensure_ascii=False,separators=(',',':'))
        if not os.path.exists(fn) or open(fn).read()!=new:      # rewrite only when the content changed (keeps git history small)
            open(fn,'w').write(new)
        cat['remote']='tools/fixed-distances/data/'+fn; cat['remoteCount']=len(cat['items']); cat['items']=[]; slimmed=True
        print('remote catalog',key,cat['remoteCount'],'items ->',fn,os.path.getsize(fn),'bytes')
if slimmed:
    json.dump(data,open('applet_data.json','w'),ensure_ascii=False,separators=(',',':'))
# the Construction Tool's geometry, cut verbatim from its template (the catalog draws cards from the tool's recipes)
CT_TPL=os.path.join(HERE,'..','construction-tool','construction-tool.template.html')
CT_CONSTS=['EPS']
CT_FUNCS=['latVec','polyGeom','curveVal','meet','evaluate','decodeRecipe','animAt']
def cut_function(src,name):
    i=src.index('function '+name+'('); j=src.index('{',i); depth=0
    for k in range(j,len(src)):
        if src[k]=='{': depth+=1
        elif src[k]=='}':
            depth-=1
            if not depth: return src[i:k+1]
    raise ValueError(name)
ct=open(CT_TPL,encoding='utf-8').read()
core='\n'.join(re.search(r'^const '+c+r' = .*$',ct,re.M).group(0) for c in CT_CONSTS)+'\n'+'\n'.join(cut_function(ct,f) for f in CT_FUNCS)
core+='\nreturn { evaluate, decodeRecipe, animAt };'
tpl=open(os.path.join(HERE,'fixed-distances.template.html'),encoding='utf-8').read()
assert tpl.count('/*CT-CORE*/')==1
js=json.dumps(data,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')
out=os.path.join(HERE,'..','..','fixed-distances.html')        # the applet, at the repository root
open(out,'w',encoding='utf-8').write(tpl.replace('/*CT-CORE*/',core.replace('</','<\\/')).replace('__DATA__',js))
print('bytes',os.path.getsize(out))
