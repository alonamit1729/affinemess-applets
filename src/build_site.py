#!/usr/bin/env python3
"""Build the website GitHub Pages serves from the readable master.

The master (src/master.html) is the ONLY place to edit the site. Everything below is
generated; never edit it by hand.

usage (from anywhere; paths are relative to the repo):
  python3 src/build_site.py            rebuild the generated files from src/master.html
  python3 src/build_site.py --check    exit 1 if a generated file differs from what the master builds
  python3 src/build_site.py --loader   also print the Squarespace loader (legacy)

  index.html         the site itself: a complete static page (served at affinemess.com/)
                     with relative links, so it works on whatever domain the repo is served from
  site/site.html     legacy, for the Squarespace loader: style + markup with absolute
  site/site.js         links to apps.affinemess.com; site.js is the script
  site/preview.html  legacy: the loader on a bare page

Placeholders in the master: __APPS__ (repo root: applets, papers), __SITE_BASE__ (site/
assets such as aom-logo.jpg), __QUORA_BASE__ (quora/). Quora rows are generated from
quora/index.json between <!--QUORA-ROWS BEGIN/END-->.
"""
import html as _h, json, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = os.path.join(REPO, 'src', 'master.html')
APPS_ABS = 'https://apps.affinemess.com/'

HEAD = '''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>affine·mess</title>
<meta name="description" content="Alon Amit's math arcade: interactive applets, Desmos sketches, talks, papers and Quora answers.">
<link rel="canonical" href="https://affinemess.com/">
<meta property="og:type" content="website">
<meta property="og:url" content="https://affinemess.com/">
<meta property="og:title" content="affine·mess">
<meta property="og:description" content="Alon Amit's math arcade: interactive applets, Desmos sketches, talks, papers and Quora answers.">
<style>html,body{margin:0;background:#F2EBD9}</style>
</head>
<body>
'''


def split_master(master=MASTER):
    src = open(master, encoding='utf-8').read()
    m = re.search(r'<script>\n(.*)</script>\s*$', src, re.S)
    assert m, 'master must end with one <script> block'
    js = m.group(1)
    html = src[:m.start()].rstrip() + '\n'
    assert '<script' not in html, 'unexpected extra <script> in the markup part'
    version = (re.search(r'VERSION (\S+)', src) or [None, 'unknown'])[1]
    return html, js, version


def markup(html, root):
    """Fill the placeholders and the Quora rows; root is '' (relative) or an absolute URL ending in /."""
    site_base, qbase = root + 'site/', root + 'quora/'
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
    return html.replace('__QUORA_BASE__', qbase).replace('__SITE_BASE__', site_base).replace('__APPS__', root)


def loader_for(version, base=APPS_ABS + 'site/'):
    return '''<!-- affinemess.com loader %(v)s: the site itself lives on GitHub Pages (%(b)s). Paste once; it never needs replacing. -->
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


def build(master=MASTER):
    html, js, version = split_master(master)
    loader = loader_for(version)
    preview = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
               '<meta name="viewport" content="width=device-width,initial-scale=1">'
               '<title>Affine Mess</title><style>body{margin:0}</style></head>\n<body>\n'
               + loader + '</body></html>\n')
    index = HEAD + markup(html, '') + '<script>\n' + js + '</script>\n</body>\n</html>\n'
    files = {
        'index.html': index,
        'site/site.html': markup(html, APPS_ABS),
        'site/site.js': js,
        'site/preview.html': preview,
    }
    return version, files, loader


def main():
    args = sys.argv[1:]
    version, files, loader = build()
    if '--check' in args:
        bad = [n for n, t in files.items()
               if not os.path.exists(os.path.join(REPO, n)) or open(os.path.join(REPO, n), encoding='utf-8').read() != t]
        if bad:
            print('OUT OF DATE: %s does not match what src/master.html builds (v%s).\n'
                  'Edit src/master.html and run python3 src/build_site.py; never edit generated files by hand.' % (', '.join(bad), version))
            sys.exit(1)
        print('generated files match src/master.html (v%s)' % version)
        return
    for name, text in files.items():
        p = os.path.join(REPO, name)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'w', encoding='utf-8') as f:
            f.write(text)
    print({'version': version, **{k: len(v) for k, v in files.items()}})
    if '--loader' in args:
        print(loader)


if __name__ == '__main__':
    main()
