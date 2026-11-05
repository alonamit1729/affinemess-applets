"""Inject data/applet_data.json into the template.  Catalogs with more than LIMIT items are written to
data/catalogs/{n}_{d}.json (fetched by the page on demand, path relative to the site root) and left empty inline."""
import json, os, copy
LIMIT=1000
data=json.load(open('applet_data.json'))
os.makedirs('catalogs',exist_ok=True)
inline=copy.deepcopy(data)
for key,cat in inline['cats'].items():
    if len(cat['items'])>LIMIT:
        n,d=key.split(','); fn=f'catalogs/{n}_{d}.json'
        json.dump(cat['items'],open(fn,'w'),ensure_ascii=False,separators=(',',':'))
        cat['remote']='tools/fixed-distances/data/'+fn; cat['remoteCount']=len(cat['items']); cat['items']=[]
        print('remote catalog',key,cat['remoteCount'],'items ->',fn,os.path.getsize(fn),'bytes')
tpl=open('fixed-distances.template.html').read()
js=json.dumps(inline,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')
open('fixed-distances.html','w').write(tpl.replace('__DATA__',js))
print('bytes',os.path.getsize('fixed-distances.html'))
