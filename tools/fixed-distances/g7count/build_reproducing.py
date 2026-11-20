"""Render REPRODUCING.md as reproducing.html, a standalone page in the site's style.
Needs pandoc on the path.   cd tools/fixed-distances/g7count && python3 build_reproducing.py"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
body = subprocess.run(['pandoc', '-f', 'gfm', '-t', 'html5', os.path.join(HERE, 'REPRODUCING.md')],
                      check=True, capture_output=True, text=True).stdout

PAGE = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Reproducing the sixteen-point classification · affine·mess</title>
<meta name="description" content="Walkthrough: how to re-run the msolve computation that classifies the planar sets of ten to sixteen points with few distances.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght,SOFT@0,9..144,300..600,0..100;1,9..144,300..500,0..100&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
:root{
  --paper:#F2EBD9; --paper-card:#FBF6E8; --ink:#1E1810; --ink-soft:#5C5142; --ink-faint:#8A7D66;
  --accent:#B85820; --accent-deep:#8A3F12; --line:#B8AB8A; --line-soft:#D4C9AC;
  --serif:'Fraunces', Georgia, 'Times New Roman', serif;
  --mono:'JetBrains Mono', ui-monospace, Menlo, Consolas, monospace;
  color-scheme: light;
}
*, *::before, *::after { box-sizing:border-box; }
html, body { margin:0; }
body { background:var(--paper); color:var(--ink); font-family:var(--serif); font-weight:380;
  font-variation-settings:"opsz" 14, "SOFT" 50; -webkit-font-smoothing:antialiased; }
.bar { display:flex; align-items:center; gap:14px; padding:10px 14px; background:var(--paper-card);
  border-bottom:1px solid var(--line); flex-wrap:wrap; }
.home { text-decoration:none; color:var(--ink); font-weight:500; font-size:17px;
  font-variation-settings:"opsz" 144, "SOFT" 30; padding-right:12px; border-right:1px solid var(--line-soft); }
.home .sep { color:var(--accent); font-style:italic; font-weight:400; margin:0 .12em; }
.bar a.lnk { font-family:var(--mono); font-size:10.5px; letter-spacing:.1em; text-transform:uppercase;
  color:var(--ink-soft); text-decoration:none; border:1px solid var(--line); border-radius:3px; padding:7px 10px; }
.bar a.lnk:hover { color:var(--accent); border-color:var(--accent); }
.page { max-width:820px; margin:0 auto; padding:26px 22px 96px; font-size:17px; line-height:1.62; }
.page h1 { font-style:italic; font-weight:400; font-size:40px; line-height:1.08; letter-spacing:-.015em;
  font-variation-settings:"opsz" 144, "SOFT" 80; margin:6px 0 18px; }
.page h2 { font-style:italic; font-weight:400; font-size:28px; line-height:1.15;
  font-variation-settings:"opsz" 144, "SOFT" 80; margin:48px 0 12px; padding-top:14px; border-top:1px solid var(--line-soft); }
.page p { margin:0 0 14px; }
.page strong { font-weight:560; }
.page a { color:var(--accent-deep); }
.page ul, .page ol { margin:0 0 14px; padding-left:24px; }
.page li { margin:0 0 6px; }
.page code { font-family:var(--mono); font-size:.8em; background:var(--paper-card); border:1px solid var(--line-soft);
  border-radius:3px; padding:1px 5px; font-variant-ligatures:none; overflow-wrap:anywhere; }
.page pre { background:var(--paper-card); border:1px solid var(--line-soft); border-radius:3px; padding:12px 14px;
  overflow-x:auto; font-size:14px; line-height:1.5; }
.page pre code { border:0; padding:0; background:none; font-size:.86em; overflow-wrap:normal; }
.page blockquote { margin:18px 0; padding:12px 18px; background:var(--paper-card); border:1px solid var(--line-soft);
  border-left:3px solid var(--accent); border-radius:3px; }
.page blockquote p:last-child { margin-bottom:0; }
.tablewrap { overflow-x:auto; margin:0 0 16px; }
.page table { border-collapse:collapse; font-size:14.5px; line-height:1.45; }
.page th, .page td { border-bottom:1px solid var(--line-soft); padding:6px 10px; text-align:left; vertical-align:top; }
.page th { font-family:var(--mono); font-size:11px; letter-spacing:.08em; text-transform:uppercase; color:var(--ink-faint); font-weight:500; }
@media (max-width:600px) { .page { font-size:16px; padding:18px 16px 72px; } .page h1 { font-size:31px; } .page h2 { font-size:24px; } }
</style>
</head><body>
<div class="bar">
  <a class="home" href="https://affinemess.com/#applets/fixed-distances" target="_top">affine<span class="sep">·</span>mess</a>
  <a class="lnk" href="https://affinemess.com/#applets/fixed-distances" target="_top">The table</a>
  <a class="lnk" href="https://affinemess.com/tools/fixed-distances/report/table-methods.pdf">How the table was computed (PDF)</a>
</div>
<main class="page">
__BODY__
</main>
</body></html>
"""
body = body.replace('<table>', '<div class="tablewrap"><table>').replace('</table>', '</table></div>')
open(os.path.join(HERE, 'reproducing.html'), 'w').write(PAGE.replace('__BODY__', body))
print('wrote reproducing.html', len(body), 'bytes of body')
