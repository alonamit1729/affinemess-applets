# affinemess-applets: rules for every Claude session

This repo is the **master copy of everything on affinemess.com**. GitHub Pages
serves it at apps.affinemess.com: the applets (`<slug>.html` in the root, plus
each applet's own tools folder), the website (`site/`, generated from `src/`),
Quora answers (`quora/`) and papers (`papers/`). Squarespace only holds a tiny
loader that fetches `site/`. There is no separate master on iCloud any more
(retired 2026-10-25); several chats and Claude Projects push here, and these
rules keep them from overwriting each other.

## The website: edit `src/`, never `site/`

- `src/master.html` is the readable source of the whole site (style, markup,
  script). `python3 src/build_site.py` regenerates `index.html` (the site
  itself, a complete static page with relative links) and the legacy
  Squarespace-loader files `site/site.html`, `site/site.js`,
  `site/preview.html`. A hand edit to any of them is lost at the next build,
  and the GitHub check (`site matches master`) fails.
- `python3 src/build_site.py --check` tells you whether `site/` is in sync.
- `python3 src/test_site.py` smoke-tests the built site in headless Chromium.
- Quora rows in the site are generated from `quora/index.json`; never
  hand-edit between `<!--QUORA-ROWS BEGIN/END-->`. `src/build_quora.py`
  rebuilds `quora/` from Alon's private Quora export (not in the repo).
- **Homepage changes are made by the "Affine Mess Website" Project.** Any
  other session that wants a card for its applet pushes the applet only and
  hands the card details (slug, title, label, one-line description, SVG
  glyph) to Alon to pass to that Project.

## Applets

- An applet's own files belong to the session working on that applet.
- Commit only the paths you meant to change: `git add <slug>.html` (and the
  applet's own folder), never `git add -A` or `git add .`.

## Before every push

1. `git fetch origin main` and rebase onto it; never force-push.
2. `git diff --stat origin/main` must list only your own paths. If it shows
   `index.html`, `site/`, `src/` or another applet and you are not the Website Project, stop.
3. Website changes: edit `src/master.html`, bump the version in both places
   (the `VERSION` header comment and `data-version` on `.am-root`; the footer
   shows it), run `build_site.py`, run `test_site.py`, then read every removed
   line in `git diff origin/main -- site/ src/` and make sure each removal is
   intended. Commit `src/`, `index.html` and `site/` together.
4. Site versions (`v2026-MM-DDx` in the footer) are assigned only by the
   Website Project: take the latest version on `origin/main` and bump the
   letter, so no two deploys share a number.
