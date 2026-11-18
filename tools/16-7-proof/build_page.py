#!/usr/bin/env python3
"""Build the applet page ../../16-7-proof.html from the executed notebook two-shapes.ipynb,
and the download bundle 16-7-proof-notebook.zip.

    python3 tools/16-7-proof/build_page.py

Needs mistune (3.x) and pygments (both come with Jupyter).  The page is one self-contained
file in the site's applet frame (am-kit v1): a slim bar, the notebook as a scrolling document
in the stage, and an About drawer.  Nothing here changes the mathematics: the notebook's text,
code and saved outputs are only re-dressed (site fonts, site colours in the figures, KaTeX for
the formulas with a plain-text fallback).  Files are rewritten only when their bytes change.

To update the applet: edit and re-run the notebook (Run All, save), then run this script.
"""
import html
import json
import re
import zipfile
from pathlib import Path

import mistune
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import PythonLexer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SLUG = "16-7-proof"
NAME = "(16, 7) proof"
BUNDLE = "16-7-proof-notebook"
BUNDLE_FILES = ["README.md", "requirements.txt", "two-shapes.ipynb", "fewdist.py", "figures.py",
                "data/L10.json", "data/solutions.json", "paper/seven-distances.pdf",
                "fifteen/fifteen.py", "fifteen/solutions9.json", "fifteen/result.json"] + ["certificates/" + name for name in (
                    "README.md", "check.py", "generate.py", "summary.json", "SHA256SUMS", "level5.jsonl",
                    "level6_sample.jsonl.gz", "level7_sample.jsonl.gz",
                    "L5.json", "L6.json", "L7.json", "L8.json", "L9.json", "L10.json")]

# ---------------------------------------------------------------- figures: notebook colours -> site colours
RECOLOUR = {
    "#fcfcfb": "#FBF6E8",   # panel surface      -> paper-card
    "#e4e2dc": "#D4C9AC",   # panel frame        -> line-soft
    "#0b0b0b": "#1E1810",   # ink
    "#52514e": "#5C5142",   # muted text         -> ink-soft
    "#c9c7c0": "#B8AB8A",   # faint marks        -> line
    "#2a78d6": "#35598F", "#eb6834": "#B85820", "#1baf7a": "#2F7F86", "#eda100": "#A8842A",
    "#e87ba4": "#5B4FA0", "#008300": "#4E7D52", "#6250d6": "#9A5577", "#e34948": "#A03A2E",
}
SVG_FONT_OLD = "font-family=\"system-ui, -apple-system, 'Segoe UI', sans-serif\""
SVG_FONT_NEW = "font-family=\"'JetBrains Mono', ui-monospace, Menlo, Consolas, monospace\""


def site_svg(svg):
    svg = re.sub(r"<\?xml[^>]*\?>", "", svg).strip()
    svg = svg.replace(";background:#fcfcfb", "")                 # let the paper show between panels
    for old, new in RECOLOUR.items():
        svg = svg.replace(old, new)
    svg = svg.replace(SVG_FONT_OLD, SVG_FONT_NEW)
    svg = svg.replace('font-size="12.5"', 'font-size="11"').replace('font-size="12"', 'font-size="11"')
    width = re.search(r'viewBox="0 0 ([\d.]+) ', svg).group(1)          # natural width in px: never enlarge
    svg = re.sub(r'(<svg[^>]*?) width="[\d.]+"', r"\1", svg, count=1)
    svg = re.sub(r'(<svg[^>]*?) style="[^"]*"', rf'\1 style="width:100%;max-width:{width}px;height:auto"', svg, count=1)
    return svg


# ---------------------------------------------------------------- formulas
SUP = str.maketrans("0123456789+-=()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ")
SUB = str.maketrans("0123456789+-=()ijknmx", "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ᵢⱼₖₙₘₓ")
WORDS = {"displaystyle": "", "left": "", "right": "", "bigl": "", "bigr": "", "qquad": "  ", "quad": " ",
         "le": "≤", "ge": "≥", "ne": "≠", "times": "×", "cdot": "·", "pm": "±", "approx": "≈", "in": "∈",
         "subseteq": "⊆", "setminus": "∖", "cup": "∪", "lfloor": "⌊", "rfloor": "⌋", "lceil": "⌈",
         "rceil": "⌉", "langle": "⟨", "rangle": "⟩", "lvert": "|", "rvert": "|", "dots": "…",
         "Delta": "Δ", "theta": "θ", "pi": "π", "omega": "ω", "sum": "Σ", "blacksquare": "∎", "det": "det",
         "sin": "sin", "cos": "cos", "tfrac": "\\tfrac", "frac": "\\frac", "sqrt": "√", "binom": "\\binom",
         "text": "", "mathbb": ""}


def plain(tex):
    """A readable plain-text stand-in for a formula, shown if KaTeX cannot be loaded."""
    s = tex
    s = re.sub(r"\\mathbb\s*\{?Q\}?", "ℚ", s)
    s = re.sub(r"\\tfrac\s*(\d)\s*(\d)", r"\1/\2", s)
    s = re.sub(r"\\t?frac\{([^{}]*)\}\{([^{}]*)\}", r"(\1)/(\2)", s)
    s = re.sub(r"\\binom\s*(\w)\s*(\w)", r"C(\1,\2)", s)
    s = re.sub(r"\\text\{([^{}]*)\}", r"\1", s)
    s = s.replace("\\{", "\x01").replace("\\}", "\x02").replace("\\#", "#").replace("\\,", " ").replace("\\;", " ")
    s = re.sub(r"\\([A-Za-z]+)", lambda m: WORDS.get(m.group(1), m.group(1)), s)
    s = s.replace("\\ ", " ").replace("\\\\", "; ").replace("\\", "")
    s = re.sub(r"\^\{([0-9+\-=()n]+)\}|\^([0-9n])", lambda m: (m.group(1) or m.group(2)).translate(SUP), s)
    s = re.sub(r"_\{([0-9+\-=()ijknmx]+)\}|_([0-9ijknmx])", lambda m: (m.group(1) or m.group(2)).translate(SUB), s)
    s = s.replace("{", "").replace("}", "") if "begin" not in s else s
    s = s.replace("\x01", "{").replace("\x02", "}")                       # braces that are part of the formula
    return re.sub(r"\s+", " ", s).strip()


def tex(source, display=False, fallback=None):
    tag, cls = ("div", "tex formula") if display else ("span", "tex")
    text = fallback if fallback is not None else plain(source)
    return f'<{tag} class="{cls}" data-tex="{html.escape(source, quote=True)}">{html.escape(text)}</{tag}>'


# ---------------------------------------------------------------- markdown
MD = mistune.create_markdown(escape=False, plugins=["table", "strikethrough"])
MATH = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)


def markdown(text):
    """Markdown -> HTML, with $...$ and $$...$$ taken out first so that Markdown cannot touch them."""
    held = []

    def hold(m):
        held.append(tex(m.group(1).strip(), True) if m.group(1) is not None else tex(m.group(2).strip()))
        return f"@@MATH{len(held) - 1}@@"

    out = MD(MATH.sub(hold, text))
    out = re.sub(r"<p>\s*(@@MATH\d+@@)\s*</p>",
                 lambda m: m.group(1) if 'class="tex formula"' in held[int(m.group(1)[6:-2])] else m.group(0), out)
    out = re.sub(r"@@MATH(\d+)@@", lambda m: held[int(m.group(1))], out)
    out = re.sub(r'<a href="(https?://[^"]+)"', r'<a href="\1" target="_blank" rel="noopener"', out)
    return out


# ---------------------------------------------------------------- code
FORMAT = HtmlFormatter(nowrap=True)


def code_block(source, cls="code"):
    return f'<pre class="{cls}">{highlight(source, PythonLexer(), FORMAT).rstrip()}</pre>'


def render_output(out):
    kind = out["output_type"]
    if kind == "stream":
        return f'<pre class="out">{html.escape("".join(out["text"]).rstrip())}</pre>'
    if kind not in ("execute_result", "display_data"):
        raise ValueError(f"unexpected output type {kind}")
    data = {k: "".join(v) if isinstance(v, list) else v for k, v in out["data"].items()}
    if "image/svg+xml" in data:
        return f'<figure class="fig">{site_svg(data["image/svg+xml"])}</figure>'
    if "text/markdown" in data:
        return f'<div class="out-md">{markdown(data["text/markdown"])}</div>'
    if "text/html" in data and 'class="highlight"' in data["text/html"]:      # a function shown in full
        plain_source = html.unescape(re.sub(r"<[^>]+>", "", re.search(
            r'<div class="highlight"><pre>(.*?)</pre></div>', data["text/html"], re.S).group(1)))
        plain_source = plain_source.lstrip("\n")
        return code_block(plain_source, "code source")
    if "text/latex" in data:
        source = re.sub(r"^\$\\displaystyle\s*|\$$", "", data["text/latex"].strip())
        return tex(source, True, fallback=data.get("text/plain", ""))
    return f'<pre class="out">{html.escape(data.get("text/plain", ""))}</pre>'


def render_notebook(nb):
    parts, sections = [], []
    for cell in nb["cells"]:
        source = "".join(cell["source"])
        if cell["cell_type"] == "markdown":
            if source.startswith("<small>*") and source.endswith("*</small>"):
                parts.append(f'<div class="cap">{markdown(source[8:-9])}</div>')
                continue
            block = markdown(source)

            def anchor(m):
                sections.append((f"s{len(sections) + 1}", re.sub(r"<[^>]+>", "", m.group(1))))
                return f'<h2 id="{sections[-1][0]}">{m.group(1)}</h2>'

            parts.append(re.sub(r"<h2>(.*?)</h2>", anchor, block))
        else:
            parts.append('<div class="cell">' + code_block(source)
                         + "".join(render_output(o) for o in cell.get("outputs", [])) + "</div>")
    return "\n".join(parts), sections


# ---------------------------------------------------------------- the page
KIT_CSS = r"""
/* am-kit v1 BEGIN */
:root{
  --paper:#F2EBD9; --paper-card:#FBF6E8; --paper-hover:#FFFAEC;
  --ink:#1E1810; --ink-soft:#5C5142; --ink-faint:#8A7D66;
  --accent:#B85820; --accent-deep:#8A3F12;
  --grid:#D8CCAE; --line:#B8AB8A; --line-soft:#D4C9AC;
  --serif:'Fraunces', Georgia, 'Times New Roman', serif;
  --mono:'JetBrains Mono', ui-monospace, Menlo, Consolas, monospace;
  color-scheme: light;
}
html, body { margin:0; height:100%; }
body {
  background:var(--paper); color:var(--ink);
  font-family:var(--serif); font-weight:380; line-height:1.5;
  font-variation-settings:"opsz" 14, "SOFT" 50;
  -webkit-font-smoothing:antialiased;
}
*, *::before, *::after { box-sizing:border-box; }
.am-app { height:100%; display:flex; flex-direction:column; overflow:hidden; }
.am-bar {
  display:flex; align-items:center; gap:10px 14px; flex-wrap:wrap; flex-shrink:0;
  padding:10px 14px; background:var(--paper-card); border-bottom:1px solid var(--line);
}
.am-spacer { flex:1 1 auto; }
.am-hint { font-size:14px; color:var(--ink-soft); line-height:1.3; }
.am-hint b { font-weight:500; font-style:italic; color:var(--ink); }
.am-label { font-family:var(--mono); font-size:10px; letter-spacing:.14em; text-transform:uppercase; color:var(--ink-faint); }
.am-stage { position:relative; flex:1 1 auto; min-height:0; }
.am-dots {
  background-color:var(--paper-card);
  background-image:radial-gradient(circle, var(--grid) 1px, transparent 1.4px);
  background-size:32px 32px;
}
.am-home {
  display:none; white-space:nowrap; text-decoration:none; color:var(--ink);
  font-family:var(--serif); font-weight:500; font-size:17px;
  font-variation-settings:"opsz" 144, "SOFT" 30;
  padding-right:12px; border-right:1px solid var(--line-soft);
}
.am-home .sep { color:var(--accent); font-style:italic; font-weight:400; margin:0 .12em; }
.am-home:hover { color:var(--accent-deep); }
html.am-standalone .am-home { display:inline-block; }
.am-btn {
  font-family:var(--mono); font-size:10.5px; letter-spacing:.1em; text-transform:uppercase;
  color:var(--ink-soft); background:transparent; border:1px solid var(--line); border-radius:3px;
  padding:8px 11px; cursor:pointer; white-space:nowrap; touch-action:manipulation;
  transition:color .2s, border-color .2s, background .2s;
}
.am-btn:hover, .am-btn:focus-visible { color:var(--accent); border-color:var(--accent); background:var(--paper-hover); outline:none; }
.am-btn:disabled { opacity:.4; cursor:default; color:var(--ink-soft); border-color:var(--line); background:transparent; }
.am-btn[aria-pressed="true"], .am-btn.is-on { color:var(--paper-card); background:var(--ink); border-color:var(--ink); }
.am-seg { display:inline-flex; border:1px solid var(--line); border-radius:3px; overflow:hidden; }
.am-seg .am-btn { border:0; border-radius:0; border-right:1px solid var(--line); }
.am-seg .am-btn:last-child { border-right:0; }
.am-input { font-family:var(--mono); font-size:13px; color:var(--ink); background:var(--paper-card); border:1px solid var(--line); border-radius:3px; padding:7px 9px; }
.am-input:focus { outline:none; border-color:var(--accent); }
input[type=range] { accent-color:var(--accent); }
.am-about {
  position:absolute; top:0; right:0; bottom:0; z-index:20;
  width:min(440px, 100%); overflow-y:auto;
  background:var(--paper-card); border-left:1px solid var(--line);
  box-shadow:-14px 0 30px -20px rgba(30,24,16,.45);
  padding:20px 24px 28px;
  transform:translateX(104%); visibility:hidden;
  transition:transform .35s cubic-bezier(.16,1,.3,1), visibility 0s linear .35s;
}
.am-about.open { transform:none; visibility:visible; transition:transform .35s cubic-bezier(.16,1,.3,1); }
.am-about h2 { font-family:var(--serif); font-style:italic; font-weight:400; font-size:26px; line-height:1.1; font-variation-settings:"opsz" 144, "SOFT" 80; letter-spacing:-.01em; margin:4px 0 14px; }
.am-about h3 { font-family:var(--mono); font-size:10px; letter-spacing:.14em; text-transform:uppercase; color:var(--ink-faint); font-weight:400; margin:22px 0 8px; }
.am-about p, .am-about li { font-size:15.5px; line-height:1.6; }
.am-about p { margin:0 0 12px; }
.am-about em { color:var(--accent-deep); }
.am-about a { color:var(--accent-deep); }
.am-about .am-close { float:right; margin:0 0 8px 12px; }
@media (max-width:600px) {
  .am-bar { padding:8px 10px; }
  .am-hint { font-size:12.5px; flex:1 1 100%; order:5; }
  .am-about { width:100%; border-left:0; padding:16px 18px 24px; }
}
/* am-kit v1 END */
"""

APPLET_CSS = r"""
/* ---- this applet: a notebook as a scrolling document ---- */
:root { --code-key:#8A3F12; --code-str:#35598F; }
a.am-btn { text-decoration:none; display:inline-block; }
#jump { max-width:min(300px, 52vw); font-size:12px; }
.am-about { width:min(520px, 100%); }
.am-about code { font-family:var(--mono); font-size:.86em; }
.doc { position:absolute; inset:0; overflow-y:auto; overflow-x:hidden; outline:none; scroll-behavior:auto; }
.page { max-width:800px; margin:0 auto; padding:26px 22px 96px; font-size:17px; line-height:1.62; }
.page h1 { font-style:italic; font-weight:400; font-size:44px; line-height:1.05; letter-spacing:-.015em;
  font-variation-settings:"opsz" 144, "SOFT" 80; margin:6px 0 6px; }
.page h1 + h3 { font-family:var(--serif); text-transform:none; letter-spacing:0; font-size:19px; font-weight:380;
  color:var(--ink-soft); margin:0 0 22px; }
.page h2 { font-style:italic; font-weight:400; font-size:30px; line-height:1.15; letter-spacing:-.01em;
  font-variation-settings:"opsz" 144, "SOFT" 80; margin:54px 0 12px; padding-top:14px; border-top:1px solid var(--line-soft); }
.page h3 { font-family:var(--mono); font-size:11.5px; letter-spacing:.13em; text-transform:uppercase; font-weight:500;
  color:var(--ink-faint); margin:32px 0 10px; }
.page p { margin:0 0 14px; }
.page strong { font-weight:560; }
.page em { color:var(--accent-deep); }
.page a { color:var(--accent-deep); }
.page ul, .page ol { margin:0 0 14px; padding-left:24px; }
.page li { margin:0 0 6px; }
.page code { font-family:var(--mono); font-size:.8em; background:var(--paper-card); border:1px solid var(--line-soft);
  border-radius:3px; padding:1px 5px; white-space:nowrap; font-variant-ligatures:none; }
.page blockquote { margin:18px 0; padding:12px 18px; background:var(--paper-card); border:1px solid var(--line-soft);
  border-left:3px solid var(--accent); border-radius:3px; }
.page blockquote p:last-child { margin-bottom:0; }
.page blockquote em { color:inherit; }
.page table { border-collapse:collapse; margin:16px 0 20px; font-size:15px; line-height:1.4; display:block; max-width:100%; overflow-x:auto; }
.page th { font-family:var(--mono); font-size:10px; letter-spacing:.12em; text-transform:uppercase; font-weight:500;
  color:var(--ink-faint); text-align:left; }
.page th, .page td { padding:7px 12px; border-bottom:1px solid var(--line-soft); vertical-align:top; }
.page td code { white-space:normal; }
.cell { margin:16px 0 20px; }
/* code lines are at most 96 characters; the font shrinks a little on narrower windows so that they still fit */
pre { margin:0; font-family:var(--mono); font-size:clamp(11px, calc((100vw - 100px) / 57.8), 12.4px); line-height:1.55; overflow-x:auto;
  font-variant-ligatures:none; font-feature-settings:"liga" 0, "calt" 0; }
pre.code { background:var(--paper-card); border:1px solid var(--line-soft); border-left:3px solid var(--line); border-radius:3px;
  padding:11px 14px; color:var(--ink); }
pre.code.source { border-left-color:var(--accent); margin-top:8px; }
pre.out { padding:8px 14px 2px; color:var(--ink-soft); white-space:pre-wrap; overflow-wrap:anywhere; }
.doc.no-code pre.code:not(.source) { display:none; }
.doc.no-code pre.out { padding-left:0; }
.code .k, .code .kn, .code .kc, .code .ow, .code .kd { color:var(--code-key); font-weight:500; }
.code .s, .code .s1, .code .s2, .code .sa, .code .sd, .code .si, .code .se { color:var(--code-str); }
.code .c1, .code .c { color:var(--ink-faint); font-style:italic; }
.code .sd { font-style:italic; }
.fig { margin:12px 0 4px; }
.fig svg { display:block; }
.out-md table { margin-top:8px; }
.cap { font-size:14.5px; line-height:1.5; color:var(--ink-soft); font-style:italic; margin:-8px 0 20px; }
.cap p { margin:0; }
.cap em { color:inherit; }
.tex.formula { display:block; margin:12px 0 16px; overflow-x:auto; overflow-y:hidden; text-align:center; }
.cell .tex.formula { text-align:left; padding-left:14px; }
.katex { font-size:1.04em; }
@media (max-width:600px) {
  .page { padding:18px 14px 80px; font-size:16px; }
  .page h1 { font-size:34px; }
  .page h2 { font-size:25px; margin-top:42px; }
  pre { font-size:11.5px; }
  .page code { white-space:normal; overflow-wrap:anywhere; }   /* a long command may break on a phone */
}
"""

KIT_JS = r"""
/* am-kit v1 BEGIN */
(function () {
  try { if (window.top === window) document.documentElement.classList.add('am-standalone'); }
  catch (_) { /* cross-origin parent: we are framed */ }
  var btn = document.querySelector('[data-am-about]');
  var pane = document.getElementById('am-about');
  if (!btn || !pane) return;
  function set(open) {
    pane.classList.toggle('open', open);
    btn.setAttribute('aria-expanded', open ? 'true' : 'false');
    pane.setAttribute('aria-hidden', open ? 'false' : 'true');
  }
  btn.addEventListener('click', function () { set(!pane.classList.contains('open')); });
  var close = pane.querySelector('.am-close');
  if (close) close.addEventListener('click', function () { set(false); });
  addEventListener('keydown', function (e) { if (e.key === 'Escape') set(false); });
  if (location.hash === '#about') set(true);
})();
/* am-kit v1 END */
"""

APPLET_JS = r"""
(function () {
  var doc = document.getElementById('doc'), jump = document.getElementById('jump'), codeBtn = document.getElementById('codeBtn');
  var heads = Array.prototype.slice.call(doc.querySelectorAll('h2[id]'));
  // section menu: scroll inside the document; the address (and the browser history) is never touched
  jump.addEventListener('change', function () {
    var el = jump.value === 'top' ? null : document.getElementById(jump.value);
    doc.scrollTop = el ? el.offsetTop - 6 : 0;
    doc.focus({ preventScroll: true });
  });
  var ticking = false;
  doc.addEventListener('scroll', function () {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(function () {
      ticking = false;
      var here = 'top', y = doc.scrollTop + 40;
      for (var i = 0; i < heads.length; i++) if (heads[i].offsetTop <= y) here = heads[i].id;
      if (jump.value !== here) jump.value = here;
    });
  }, { passive: true });
  // show / hide the code cells, keeping the reader's place
  codeBtn.addEventListener('click', function () {
    var on = codeBtn.getAttribute('aria-pressed') !== 'true';
    var y = doc.scrollTop + 40, anchor = null;
    var marks = doc.querySelectorAll('.page > h1, .page > h2, .page > h3, .page > p, .page > blockquote, .fig');
    for (var i = 0; i < marks.length; i++) { if (marks[i].offsetTop <= y) anchor = marks[i]; else break; }
    var before = anchor ? anchor.getBoundingClientRect().top : 0;
    codeBtn.setAttribute('aria-pressed', on ? 'true' : 'false');
    doc.classList.toggle('no-code', !on);
    if (anchor) doc.scrollTop += anchor.getBoundingClientRect().top - before;
  });
  // formulas: KaTeX when it can be loaded, the plain-text stand-ins otherwise
  function typeset() {
    if (!window.katex) return;
    Array.prototype.forEach.call(document.querySelectorAll('.tex'), function (el) {
      try { katex.render(el.dataset.tex, el, { throwOnError: false, displayMode: el.classList.contains('formula') }); } catch (e) {}
    });
  }
  var s = document.createElement('script');
  s.src = 'https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.11/katex.min.js';
  s.onload = typeset; s.onerror = function () {};
  document.head.appendChild(s);
})();
"""

ABOUT = """
      <button class="am-btn am-close">Close</button>
      <h2>(16, 7) proof</h2>
      <p>Sixteen points in the plane can be placed so that only seven different distances occur between them, and
      seventeen cannot. This page is a complete proof that, up to similarity, there are <em>exactly two</em> ways to
      place the sixteen. Both are pieces of the triangular lattice.</p>
      <h3>What you are looking at</h3>
      <p>A Jupyter notebook together with its results: the text is the argument, the code cells are the exact
      computations it rests on, and the pictures and numbers are what those computations produced.
      <em>Code</em> in the bar hides or shows the cells, and the menu jumps between sections.</p>
      <p>The proof is computer-assisted. Its last section lists precisely what it takes from the literature and
      what is computed outside the notebook; everything is exact, in integers, rational numbers and small number
      fields, with no floating point in any decision.</p>
      <h3>Run it yourself</h3>
      <p><em>Notebook</em> downloads the bundle: the notebook, the Python module behind it, the lists of patterns
      and the program that checks certificates. It needs Python 3.10 or later and SymPy, nothing else. Unzip,
      <code>pip install -r requirements.txt</code>, open <code>two-shapes.ipynb</code> in Jupyter and run all
      cells; that takes about three and a half minutes. Sections 3.5 and 8.3 say how to run the longer computations
      again.</p>
      <h3>Where it comes from</h3>
      <p>The bound of sixteen points is from A. Amit and GPT-6 Astra, <i>The maximum cardinality of a planar
      seven-distance set</i> (Hexagon 2610.00534). The way of listing the ten-point patterns, one point at a time, is
      due to F. Szöllősi and P. Östergård, <i>Constructions of maximum few-distance sets in Euclidean spaces</i>
      (2020); here every pattern that is discarded on the way has a certificate that can be checked. The two sets
      themselves were found by P. Erdős and P. Fishburn in 1996, who conjectured that they are the only ones.</p>
      <p>The files are in the site's repository, under
      <a href="https://github.com/alonamit1729/affinemess-applets/tree/main/tools/16-7-proof" target="_blank" rel="noopener">tools/16-7-proof</a>.</p>
"""


def page(body, sections):
    options = '<option value="top">Top</option>' + "".join(
        f'<option value="{sid}">{html.escape(html.unescape(title))}</option>' for sid, title in sections)
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{NAME} · affine·mess</title>
<meta name="description" content="A runnable, exact proof that exactly two sets of sixteen points in the plane determine only seven distinct distances.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght,SOFT@0,9..144,300..600,0..100;1,9..144,300..500,0..100&family=JetBrains+Mono:wght@400;500&display=swap">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.11/katex.min.css">
<style>{KIT_CSS}{APPLET_CSS}</style></head>
<body>
<!-- generated by tools/{SLUG}/build_page.py from tools/{SLUG}/two-shapes.ipynb; edit those, not this file -->
<div class="am-app">
  <header class="am-bar">
    <a class="am-home" href="https://affinemess.com/#applets/{SLUG}" target="_top">affine<span class="sep">·</span>mess</a>
    <select class="am-input" id="jump" aria-label="Jump to a section">{options}</select>
    <button class="am-btn" id="codeBtn" aria-pressed="true" title="Show or hide the code cells">Code</button>
    <span class="am-spacer"></span>
    <a class="am-btn" href="tools/{SLUG}/{BUNDLE}.zip" download title="Download the notebook, its module and data">Notebook</a>
    <button class="am-btn" data-am-about aria-expanded="false" aria-controls="am-about">About</button>
  </header>
  <div class="am-stage">
    <main class="doc" id="doc" tabindex="0"><div class="page">
{body}
    </div></main>
    <aside class="am-about" id="am-about" aria-hidden="true">{ABOUT}    </aside>
  </div>
</div>
<script>{KIT_JS}</script>
<script>{APPLET_JS}</script>
</body></html>
"""


def write_if_changed(path, data):
    data = data if isinstance(data, bytes) else data.encode("utf-8")
    if path.exists() and path.read_bytes() == data:
        print("unchanged:", path.relative_to(ROOT))
        return
    path.write_bytes(data)
    print(f"wrote {path.relative_to(ROOT)}  ({len(data):,} bytes)")


def bundle():
    """The download: a zip with a fixed timestamp and order, so equal contents give equal bytes."""
    import io
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name in BUNDLE_FILES:
            info = zipfile.ZipInfo(f"{BUNDLE}/{name}", date_time=(2026, 11, 16, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, (HERE / name).read_bytes())
    return buffer.getvalue()


if __name__ == "__main__":
    nb = json.loads((HERE / "two-shapes.ipynb").read_text(encoding="utf-8"))
    body, sections = render_notebook(nb)
    write_if_changed(ROOT / f"{SLUG}.html", page(body, sections))
    write_if_changed(HERE / f"{BUNDLE}.zip", bundle())
