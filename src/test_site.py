#!/usr/bin/env python3
"""Smoke-test the built site in headless Chromium before pushing.

usage: python3 src/test_site.py [--shots DIR]
Serves the repo on a free local port, loads index.html (the site) and site/preview.html
(the legacy Squarespace loader, base URL pointed at that server), and checks: version stamp, every tab, every applet tile, the in-place
viewer (#applets/<slug>), Quora search, and that no page errors are thrown.
Needs: pip install playwright && playwright install chromium (or an existing Chromium).
"""
import functools, glob, http.server, os, re, shutil, socket, sys, tempfile, threading

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    from playwright.sync_api import sync_playwright
    shots = sys.argv[sys.argv.index('--shots') + 1] if '--shots' in sys.argv else None
    tmp = tempfile.mkdtemp()
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0)); port = s.getsockname()[1]
    base = 'http://127.0.0.1:%d/' % port
    # a copy of preview.html whose loader fetches from the local server
    prev = open(os.path.join(REPO, 'site', 'preview.html'), encoding='utf-8').read()
    open(os.path.join(tmp, 'preview.html'), 'w', encoding='utf-8').write(
        prev.replace('https://apps.affinemess.com/site/', base + 'site/'))
    site_html = open(os.path.join(REPO, 'site', 'site.html'), encoding='utf-8').read()

    class H(http.server.SimpleHTTPRequestHandler):
        def translate_path(self, path):
            if path.split('?')[0] == '/_test/preview.html':
                return os.path.join(tmp, 'preview.html')
            return super().translate_path(path)
        def log_message(self, *a): pass
    srv = http.server.ThreadingHTTPServer(('127.0.0.1', port), functools.partial(H, directory=REPO))
    threading.Thread(target=srv.serve_forever, daemon=True).start()

    fails = []
    def check(ok, what):
        print(('  ok  ' if ok else '  FAIL') + ' ' + what)
        if not ok: fails.append(what)

    version = re.search(r'data-version="([^"]+)"', site_html).group(1)
    slugs = re.findall(r'<a class="tile" data-slug="([^"]+)" href="(?:https://apps\.affinemess\.com/)?([a-z0-9-]+)\.html"', site_html)
    exe = None
    for pat in ('/opt/pw-browsers/chromium-*/chrome-linux*/chrome',):
        hits = sorted(glob.glob(pat)); exe = hits[-1] if hits else None
    with sync_playwright() as p:
        b = p.chromium.launch(**({'executable_path': exe} if exe else {}))
        for entry, w, h in (('index.html', 1280, 900), ('index.html', 390, 844), ('_test/preview.html', 1280, 900)):
            pg = b.new_page(viewport={'width': w, 'height': h})
            errs = []; pg.on('pageerror', lambda e: errs.append(str(e)))
            url = base + entry
            pg.goto(url); pg.wait_for_selector('.am-root .tile')
            print('%s at %dx%d' % (entry, w, h))
            check(pg.get_attribute('.am-root', 'data-version') == version, 'version stamp v' + version)
            tiles = pg.eval_on_selector_all('#am-panel-applets .tile', 'els=>els.map(e=>e.dataset.slug)')
            check(len(tiles) >= 1, '%d applet tiles' % len(tiles))
            for slug, f in slugs:
                if slug in tiles:
                    check(os.path.exists(os.path.join(REPO, f + '.html')), 'applet file %s.html exists' % f)
            tabs = pg.eval_on_selector_all('.am-root nav button, .am-root [role=tab]', 'els=>els.map(e=>e.textContent.trim())')
            check(len(tabs) >= 6, 'tabs: ' + ', '.join(tabs))
            first = tiles[0] if tiles else None
            if first:
                pg.goto(url + '#applets/' + first); pg.wait_for_timeout(800)
                srcs = pg.eval_on_selector_all('iframe', 'els=>els.map(e=>e.getAttribute("src"))')
                check(any(first in (s or '') for s in srcs), 'viewer opens #applets/' + first)
                if entry == 'index.html':
                    fr = pg.frame_locator('iframe').locator('body')
                    check(fr.count() == 1, 'applet loads from this site (same origin): ' + ', '.join(srcs))
                    pg.goto(url + '#papers'); pg.wait_for_timeout(300)
                    pdf = pg.eval_on_selector_all('a.paper', 'els=>els.map(e=>e.getAttribute("href")).filter(h=>h.indexOf("papers/")>=0)')
                    check(all(os.path.exists(os.path.join(REPO, h)) for h in pdf) and pdf, 'thesis link resolves locally: ' + ', '.join(pdf))
                    pg.goto(url + '#media'); pg.wait_for_timeout(1200)
                    logo = pg.eval_on_selector_all('.media.podcast img', 'els=>els.map(e=>e.naturalWidth)')
                    check(logo and all(logo), 'podcast logos load')
            pg.goto(url + '#quora'); pg.wait_for_timeout(400)
            n = pg.eval_on_selector_all('#am-panel-quora .qrow', 'els=>els.length')
            check(n > 0, '%d Quora rows' % n)
            if shots:
                os.makedirs(shots, exist_ok=True)
                pg.goto(url); pg.wait_for_selector('.am-root .tile')
                pg.screenshot(path=os.path.join(shots, '%s-%d.png' % (entry.replace('/', '_'), w)))
            check(not errs, 'no page errors' + ('' if not errs else ': ' + '; '.join(errs[:3])))
            pg.close()
        b.close()
    srv.shutdown(); shutil.rmtree(tmp, ignore_errors=True)
    print('PASS' if not fails else 'FAILED: %d check(s)' % len(fails))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
