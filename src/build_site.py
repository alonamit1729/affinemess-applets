#!/usr/bin/env python3
"""Build the website files GitHub Pages serves from the readable master.

The master (src/master.html) is the ONLY place to edit the site. site/site.html,
site/site.js and site/preview.html are generated; never edit them by hand.

usage (from anywhere; paths are relative to the repo):
  python3 src/build_site.py            rebuild site/ from src/master.html
  python3 src/build_site.py --check    exit 1 if site/ differs from what the master builds
  python3 src/build_site.py --loader   also print the one-time Squarespace loader (already pasted)

  site/site.html     style + markup (everything except the script); Quora rows generated from quora/index.json
  site/site.js       the script (runs after site.html is injected)
  site/preview.html  standalone page containing the loader (the whole site without Squarespace)
"""
import html as _h, json, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = os.path.join(REPO, 'src', 'master.html')
BASE = 'https://apps.affinemess.com/site/'


def build(master=MASTER, base=BASE):
    src = open(master, encoding='utf-8').read()
    m = re.search(r'<script>\n(.*)</script>\s*$', src, re.S)
    assert m, 'master must end with one <script> block'
    js = m.group(1)
    html = src[:m.start()].rstrip() + '\n'
    assert '<script' not in html, 'unexpected extra <script> in the markup part'

    # --- Quora rows, generated from quora/index.json (written by build_quora.py) ---
    qbase = base.rstrip('/').rsplit('/', 1)[0] + '/quora/'
    qidx = os.path.join(REPO, 'quora', 'index.json')
    if os.path.exists(qidx):
        items = json.load(open(qidx, encoding='utf-8'))
        e = lambda t: _h.escape(t, quote=True)
        out_rows = []
        for it in items:
            find = (it['title'] + ' ' + it['excerpt']).lower()
            out_rows.append(
                '          <a class="qrow" data-slug="%s" data-name="%s" data-excerpt="%s" data-find="%s" href="%s%s.html" target="_blank" rel="noopener">'
                '<span class="yr">%s</span><span class="meta"><span class="title">%s</span><span class="sub">%s</span></span></a>'
                % (e(it['slug']), e(it['title']), e(it['excerpt']), e(find), qbase, it['slug'], it['date'][:4], it['th'], e(it['excerpt'])))
        rows = ('        <a class="qotd" href="#quora" data-slug="x"><div class="label">Answer of the day</div><div class="title"></div><div class="desc"></div><div class="open">read <span>→</span></div></a>\n'
                '        <div class="qtools"><input type="search" placeholder="Search the text of %d answers…" aria-label="Search answers"><button type="button" data-q="random">Random</button><span class="qcount"></span></div>\n'
                '        <div class="qlist">\n' % len(items) + '\n'.join(out_rows) + '\n          <div class="qnone" hidden>No answers match.</div>\n        </div>\n')
        html = re.sub(r'(<!--QUORA-ROWS BEGIN-->\n).*?(        <!--QUORA-ROWS END-->)', lambda m_: m_.group(1) + rows + m_.group(2), html, flags=re.S)
    html = html.replace('__QUORA_BASE__', qbase).replace('__SITE_BASE__', base)

    version = (re.search(r'VERSION (\S+)', src) or [None, 'unknown'])[1]

    loader = '''<!-- affinemess.com loader %(v)s: the site itself lives on GitHub Pages (%(b)s). Paste once; it never needs replacing. -->
<div id="am-mount"><noscript>affine&middot;mess needs JavaScript to show its pages.</noscript></div>
<script>
(function () {
  var B = '%(b)s', m = document.getElementById('am-mount');
  function fail() {
    m.innerHTML = '<p style="font:17px Georgia,serif;padding:48px 20px;text-align:center;color:#1E1810">affine&middot;mess is taking a moment to load. <a href="' + B + 'preview.html">Open it directly</a>.</p>';
  }
  function get(f) {
    return fetch(B + f, { cache: 'no-cache' }).then(function (r) { if (!r.ok) throw new Error(f); return r.text(); });
  }
  Promise.all([get('site.html'), get('site.js')]).then(function (a) {
    m.innerHTML = a[0];
    var s = document.createElement('script');
    s.textContent = a[1];
    document.body.appendChild(s);
  }).catch(fail);
})();
</script>
''' % {'v': version, 'b': base}

    preview = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
               '<meta name="viewport" content="width=device-width,initial-scale=1">'
               '<title>Affine Mess</title><style>body{margin:0}</style></head>\n<body>\n'
               + loader + '</body></html>\n')
    return version, {'site.html': html, 'site.js': js, 'preview.html': preview}, loader


def main():
    args = sys.argv[1:]
    version, files, loader = build()
    site = os.path.join(REPO, 'site')
    if '--check' in args:
        bad = []
        for name, text in files.items():
            p = os.path.join(site, name)
            if not os.path.exists(p) or open(p, encoding='utf-8').read() != text:
                bad.append('site/' + name)
        if bad:
            print('OUT OF DATE: %s does not match what src/master.html builds (v%s).\n'
                  'Edit src/master.html and run python3 src/build_site.py; never edit site/ by hand.' % (', '.join(bad), version))
            sys.exit(1)
        print('site/ matches src/master.html (v%s)' % version)
        return
    os.makedirs(site, exist_ok=True)
    for name, text in files.items():
        with open(os.path.join(site, name), 'w', encoding='utf-8') as f:
            f.write(text)
    print({'version': version, **{k: len(v) for k, v in files.items()}})
    if '--loader' in args:
        print(loader)


if __name__ == '__main__':
    main()
