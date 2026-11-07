"""Inject data/applet_data.json into the template.  Cells with more than LIMIT items are *slimmed*: their items are written
to data/catalogs/{n}_{d}.json (in their original order, so card numbers are stable; fetched by the page on demand, path
relative to the site root) and removed from applet_data.json, which keeps only the counts and a `remote` pointer.
`datautil.load()` merges them back for scripts that need every item.  Idempotent: an already slim applet_data.json
is left as it is, and its catalogs are not rewritten."""
import json, os
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
tpl=open('fixed-distances.template.html').read()
js=json.dumps(data,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')
open('fixed-distances.html','w').write(tpl.replace('__DATA__',js))
print('bytes',os.path.getsize('fixed-distances.html'))
