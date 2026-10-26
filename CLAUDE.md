# affinemess-applets: rules for every Claude session

This repo serves two things on GitHub Pages (apps.affinemess.com): the
applets (`<slug>.html` in the root, plus each applet's own tools folder)
and the website itself (`site/`). Several chats and Claude Projects push
here. These rules keep them from overwriting each other.

## Who may change what

- **`site/` belongs to the "Affine Mess Website" Project only.** Its files
  are generated from a master file by a build script; a hand edit here is
  lost at the next build. Any other session that wants a homepage change
  (a new card, a changed description) does not edit `site/`: it pushes its
  applet and hands the card details (slug, title, one-line description,
  label, SVG glyph) to Alon to pass to that Project.
- **An applet's own files belong to the session working on that applet.**
  Commit only the paths you meant to change: `git add <slug>.html` (and the
  applet's own folder), never `git add -A` or `git add .`.

## Before every push

1. `git fetch origin main` and rebase onto it; never force-push.
2. `git diff --stat origin/main` must list only your own paths. If it shows
   `site/` (or another applet) and you are not the Website Project, stop.
3. Website Project only: before editing the master, build it **unmodified**
   and compare with `origin/main:site/site.html`. Any difference means the
   live site has something the master lacks; port it into the master first.
   After building the change, read every removed line in
   `git diff origin/main -- site/` and make sure each removal is intended.
4. Site versions (`v2026-MM-DDx` in the footer) are assigned only by the
   Website Project: take the latest version on `origin/main` and bump the
   letter, so no two deploys share a number.
