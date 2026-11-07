"""Load data/applet_data.json with the items of the big cells restored.

The committed applet_data.json is *slim*: cells with more than 1000 items keep only their counts and a
`remote` pointer, and their items live in data/catalogs/<n>_<d>.json (the files the live page fetches), so the
same 30 MB are not stored twice.  Scripts that need every item call

    import datautil; D = datautil.load('data/applet_data.json')

and get exactly the full structure `build_data2.py` produces (items in their original order).
"""
import json, os

def load(path):
    D = json.load(open(path))
    here = os.path.dirname(os.path.abspath(path))
    for key, cat in D['cats'].items():
        if cat.get('remote') and not cat['items']:
            fn = None
            for d in (os.path.join(here, 'catalogs'), os.path.join(here, 'data', 'catalogs')):
                f = os.path.join(d, os.path.basename(cat['remote']))
                if os.path.exists(f): fn = f; break
            if fn is None: raise FileNotFoundError('catalog for cell %s not found next to %s' % (key, path))
            cat['items'] = json.load(open(fn))
            assert len(cat['items']) == cat['remoteCount'], key
    return D
