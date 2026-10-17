  (function() {
    const root = document.querySelector('.am-root');
    if (!root) return;
    try { console.info('affinemess.com site version ' + root.dataset.version); } catch (_) {}
    const tabs   = root.querySelectorAll('.tabs button[data-tab]');
    const panels = root.querySelectorAll('.panel');

    /* ---------- viewer ----------
       Any .tile[data-slug] opens in place at #<tab>/<slug>, where <tab> is the panel it sits in.
       The iframe loads data-src if present (YouTube/Desmos embed URLs), otherwise the tile's href.
       "New tab" always uses href, so changing a tile's URL needs no script edit. */
    const viewer   = document.getElementById('am-viewer');
    const stage    = document.getElementById('am-viewer-stage');
    let   frame    = null;      // created fresh on every open, removed on close (see mountFrame)
    const vTitle   = document.getElementById('am-viewer-title');
    const vNote    = document.getElementById('am-viewer-note');
    const vOpen    = document.getElementById('am-viewer-open');
    const vBack    = viewer.querySelector('[data-viewer="back"]');
    const vCopy    = viewer.querySelector('[data-viewer="copy"]');
    const overlayMQ = window.matchMedia('(max-width: 700px), (max-height: 520px)');
    const siteTitle = document.title;
    let currentKey = null;      // "<tab>/<slug>" of the open item
    let currentTab = null;
    let viewerPushed = false;   // true if WE pushed a history entry for the open viewer
    let lastTile = null;

    function tabOf(tile) {
      const p = tile.closest('.panel');
      return p ? p.id.replace('am-panel-', '') : null;
    }
    function tabLabel(tab) {
      const b = Array.from(tabs).find(x => x.dataset.tab === tab);
      return b ? b.textContent.trim() : tab;
    }
    function tileFor(tab, slug) {
      const panel = document.getElementById('am-panel-' + tab);
      if (!panel) return null;
      return Array.from(panel.querySelectorAll('.tile[data-slug]'))
                  .find(t => t.dataset.slug === slug) || null;
    }

    function openViewer(tile, push) {
      const tab = tabOf(tile), key = tab + '/' + tile.dataset.slug;
      if (push) {
        history.pushState({ amViewer: 1 }, '', '#' + key);
        viewerPushed = true;
      }
      if (currentKey === key && !viewer.hidden) return;
      const name = (tile.querySelector('.title') || tile).textContent.trim();
      currentKey = key;
      currentTab = tab;
      lastTile = tile;
      vTitle.textContent = name;
      vBack.textContent = '← ' + tabLabel(tab);
      vOpen.href = tile.getAttribute('href');
      vNote.textContent = '';
      if (tile.dataset.note) vNote.appendChild(document.createElement('span')).textContent = tile.dataset.note;
      vNote.hidden = !tile.dataset.note;
      vNote.classList.remove('open');
      viewer.classList.toggle('is-video', tab === 'youtube');
      mountFrame(tile.dataset.src || tile.getAttribute('href'), name);
      viewer.hidden = false;
      root.classList.add('am-viewing');
      if (overlayMQ.matches) document.documentElement.classList.add('am-lock');
      document.title = name + ' · ' + siteTitle;
    }

    // UI-only teardown (no history changes)
    function closeViewerUI() {
      if (viewer.hidden && !currentKey) return;
      viewer.hidden = true;
      root.classList.remove('am-viewing');
      document.documentElement.classList.remove('am-lock');
      unmountFrame();                        // destroys the applet / stops the video
      vNote.hidden = true;
      currentKey = null;
      currentTab = null;
      viewerPushed = false;
      document.title = siteTitle;
    }

    // User-initiated close (Back button, Esc): unwind our history entry if we made one
    function closeViewer() {
      const tile = lastTile, tab = currentTab || 'applets';
      if (viewerPushed) {
        history.back();                      // popstate -> route() -> closeViewerUI()
      } else {
        history.replaceState(null, '', '#' + tab);
        closeViewerUI();
      }
      if (tile) setTimeout(() => tile.focus({ preventScroll: true }), 0);
    }

    /* A brand-new iframe per open: re-pointing an existing iframe adds its own session-history
       entries, which would make the browser's Back button navigate *inside* the frame. */
    function mountFrame(url, name) {
      unmountFrame();
      const f = document.createElement('iframe');
      f.title = name;
      f.setAttribute('allow', 'autoplay; fullscreen; encrypted-media; picture-in-picture');
      f.setAttribute('allowfullscreen', '');
      f.addEventListener('load', () => f.classList.add('loaded'));
      f.src = url;
      stage.appendChild(f);
      frame = f;
    }
    function unmountFrame() {
      if (frame) { frame.remove(); frame = null; }
    }

    root.querySelectorAll('.tile[data-slug]').forEach(tile => {
      tile.addEventListener('click', e => {
        // let ctrl/cmd/shift/middle-click do their normal new-tab thing
        if (e.metaKey || e.ctrlKey || e.shiftKey || e.altKey || e.button !== 0) return;
        e.preventDefault();
        openViewer(tile, true);
      });
    });

    vNote.addEventListener('click', () => vNote.classList.toggle('open'));   // phone: tap to expand

    viewer.addEventListener('click', e => {
      const b = e.target.closest('[data-viewer]');
      if (!b) return;
      if (b.dataset.viewer === 'back') closeViewer();
      if (b.dataset.viewer === 'copy') {
        const url = location.origin + location.pathname + '#' + currentKey;
        const done = () => {
          vCopy.textContent = 'Copied ✓';
          setTimeout(() => { vCopy.textContent = 'Copy link'; }, 1600);
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(url).then(done, () => window.prompt('Copy this link:', url));
        } else {
          window.prompt('Copy this link:', url);
        }
      }
    });

    /* ---------- tabs ---------- */
    function activate(name, pushHash) {
      let found = false;
      tabs.forEach(btn => {
        const isMatch = btn.dataset.tab === name;
        btn.setAttribute('aria-selected', isMatch ? 'true' : 'false');
        if (isMatch) found = true;
      });
      if (!found) return false;
      if (currentTab && name !== currentTab) closeViewerUI();
      panels.forEach(p => p.classList.toggle('active', p.id === 'am-panel-' + name));
      if (pushHash) history.replaceState(null, '', '#' + name);
      return true;
    }

    tabs.forEach(btn => {
      btn.addEventListener('click', e => {
        e.preventDefault();
        closeViewerUI();                     // clicking a tab (even the current one) leaves the viewer
        activate(btn.dataset.tab, true);
      });
    });

    document.addEventListener('keydown', e => {
      if (e.target.matches('input, textarea')) return;
      if (e.key === 'Escape' && !viewer.hidden) { closeViewer(); return; }
      if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return;
      const order = Array.from(tabs).map(b => b.dataset.tab);
      const current = Array.from(tabs).findIndex(b => b.getAttribute('aria-selected') === 'true');
      const next = e.key === 'ArrowRight' ? (current + 1) % order.length : (current - 1 + order.length) % order.length;
      closeViewerUI();
      activate(order[next], true);
    });

    /* ---------- routing: #applets, #applets/<slug>, #youtube/<slug>, ... ---------- */
    function route() {
      const h = (location.hash || '').replace(/^#\/?/, '');
      const m = h.match(/^([a-z]+)\/(.+?)\/*$/);
      if (m) {
        let slug = m[2];
        try { slug = decodeURIComponent(slug); } catch (_) {}
        const tile = tileFor(m[1], slug);
        if (tile) { activate(m[1], false); openViewer(tile, false); return; }
        closeViewerUI();                     // unknown item: just show its tab
        if (activate(m[1], false)) history.replaceState(null, '', '#' + m[1]);
        else activate('applets', false);
        return;
      }
      closeViewerUI();
      if (!activate(h || 'applets', false)) activate('applets', false);
    }
    window.addEventListener('popstate', route);
    window.addEventListener('hashchange', route);
    route();

    // Keep this LAST: proves the whole script arrived and ran.
    const ver = root.querySelector('.am-ver');
    if (ver) { ver.textContent = 'v' + root.dataset.version; ver.classList.remove('am-ver-bad'); }
  })();
