# affinemess-applets: rules for every Claude session

This repo is the **master copy of everything on affinemess.com**, and GitHub
Pages serves it directly at **https://affinemess.com** (since 2026-10-25; the
domain is in `CNAME`, DNS at Squarespace Domains): the site itself
(`index.html`, generated from `src/`), the applets (`<slug>.html` in the root,
plus each applet's own tools folder), Quora answers (`quora/`) and papers
(`papers/`). The old address apps.affinemess.com is served by the separate
repo `alonamit1729/apps-redirect`, which forwards every path to the same path
on affinemess.com. Squarespace no longer hosts anything. Several chats and
Claude Projects push here; these rules keep them from overwriting each other.

## The website: edit `src/`, never `index.html` or `site/`

- `src/master.html` is the readable source of the whole site (style, markup,
  script). `python3 src/build_site.py` regenerates `index.html` (the site
  itself, a complete static page with relative links) and the legacy
  Squarespace-loader files `site/site.html`, `site/site.js`,
  `site/preview.html` (to be retired once Squarespace is cancelled;
  `site/aom-logo.jpg` is a real asset and stays). A hand edit to any of them is lost at the next build,
  and the GitHub check (`site matches master`) fails.
- `python3 src/build_site.py --check` tells you whether the generated files are in sync.
- Links inside the site are relative; applets link home to
  `https://affinemess.com/#applets/<slug>`.
- `python3 src/test_site.py` smoke-tests the built site in headless Chromium.
- Quora rows in the site are generated from `quora/index.json`; never
  hand-edit between `<!--QUORA-ROWS BEGIN/END-->`. `src/build_quora.py`
  rebuilds `quora/` from Alon's private Quora export (not in the repo).
- **Any session may change the website** (`src/master.html`: homepage cards,
  header, head tags) when Alon asks, provided it (a) keeps the site's look and
  feel: the `affinemess-applet-style` tokens, fonts and components, no new
  colours or fonts; and (b) follows the concurrency rules below, so it never
  overrides a session working at the same time.
- Keep website commits small and focused: one change per commit, touching
  only the lines it needs, so concurrent sessions rarely edit the same lines.
- A homepage card needs: slug, card title, label, one-line description, and a
  small SVG glyph (64x64 viewBox, stroke-only, site colours).

## Applets

- An applet's own files belong to the session working on that applet.
- Commit only the paths you meant to change: `git add <slug>.html` (and the
  applet's own folder), never `git add -A` or `git add .`.

## Large files

GitHub Pages publishes everything on `main`, and every committed version of a file stays in history for every
clone. So:

- Don't commit generated or intermediate files over **1 MB** under `tools/`. Exceptions: `tools/fixed-distances/data/catalogs/*`
  (fetched live by the applet) and `tools/fixed-distances/data/applet_data.json` (slim: the big cells' items live in
  the catalogs; `tools/fixed-distances/datautil.py` merges them back). Anything else that big goes behind a
  `fetch_data.py` + `DATA.md` pointer, as in `tools/fixed-distances/DATA.md` (files stay in history at a pinned
  commit; never rewrite history, the pointer depends on it).
- Rebuild generated data only when its content changed, and commit it only then (the builders write a file only if
  the bytes differ). A 30 MB catalog rewritten on every rebuild adds megabytes to every future clone.
- The GitHub check `size-check` flags (red X and email) a push that adds or changes a file over 5 MB outside the exceptions; it flags, it does not block.

## Before every push

1. `git fetch origin main` and rebase onto it; never force-push.
2. `git diff --stat origin/main` must list only the paths your change
   touched (for a website change: `src/`, `index.html`, `site/`, plus any new
   assets). Anything else means you picked up stray changes: stop.
3. Website changes: edit `src/master.html`, bump the version in both places
   (the `VERSION` header comment and `data-version` on `.am-root`; the footer
   shows it), run `build_site.py`, run `test_site.py`, then read every removed
   line in `git diff origin/main -- site/ src/` and make sure each removal is
   intended. Commit `src/`, `index.html` and `site/` together.
4. Site versions (`v2026-MM-DDx` in the footer) are assigned at push time by
   whichever session pushes a website change: after your final rebase, take
   the latest version on `origin/main` and bump the letter (a new date starts
   at `a`), so no two deploys share a number.

## Concurrency: other sessions may be pushing at the same moment

- Fetch and rebase onto `origin/main` right before committing, not only at
  the start of the session.
- Generated files (`index.html`, `site/`) are never hand-merged: on a
  conflict, take `origin/main`'s version, rerun `python3 src/build_site.py`,
  and commit the result.
- A rebase conflict in `src/` or any hand-written file, where another session
  changed the same lines: stop and ask Alon. Never discard another session's
  lines to make your change apply.
- If the push is rejected (non-fast-forward), someone pushed meanwhile:
  fetch, rebase, rebuild, re-test, re-bump the version, retry. Never force-push.
- After pushing, confirm both GitHub checks (`site matches master`,
  `size-check`) went green.
