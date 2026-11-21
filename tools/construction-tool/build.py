"""Build construction-tool.html from the template and the catalog fingerprints.
   python3 canon.py <fixed-distances applet_data.json> catalog.json   (only when the catalog changes)
   python3 build.py   -> ../../construction-tool.html"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(here, 'construction-tool.template.html'), encoding='utf-8').read()
cat = open(os.path.join(here, 'catalog.json'), encoding='utf-8').read()
a, b = '/*CATALOG*/null/*END-CATALOG*/', '/*CATALOG*/' + cat + '/*END-CATALOG*/'
assert t.count(a) == 1
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..', '..', 'construction-tool.html')
open(out, 'w', encoding='utf-8').write(t.replace(a, b))
print('wrote', out, os.path.getsize(out))
