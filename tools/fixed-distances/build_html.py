import json
data=json.load(open('applet_data.json'))
tpl=open('fixed-distances.template.html').read()
js=json.dumps(data,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')
open('fixed-distances.html','w').write(tpl.replace('__DATA__',js))
import os; print('bytes',os.path.getsize('fixed-distances.html'))
