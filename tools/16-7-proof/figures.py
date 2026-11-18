"""figures -- small dependency-free SVG helpers for the notebook.

Every drawing is a list of *panels*; a panel is a point set with optional coloured segments.
Segments are coloured by the rank of their length (shortest = first colour), and a legend
with the actual values accompanies each figure, so colour is never the only carrier.
"""
from __future__ import annotations

import math
from html import escape
from itertools import combinations

# categorical colours, in a fixed order checked for colour-vision separation
EDGE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, FAINT, SURFACE, FRAME = "#0b0b0b", "#52514e", "#c9c7c0", "#fcfcfb", "#e4e2dc"
ACCENT = "#eb6834"
FONT = "font-family=\"system-ui, -apple-system, 'Segoe UI', sans-serif\""


class SVG:
    """An SVG string that Jupyter displays as a picture."""

    def __init__(self, text):
        self.text = text

    def _repr_svg_(self):
        return self.text

    def save(self, path):
        with open(path, "w") as fh:
            fh.write(self.text)


def _fit(xy, box, margin, top=0, bottom=0):
    """Map float points into a square box of side `box`, keeping the aspect ratio and leaving
    `top` / `bottom` pixels free for a title / subtitle."""
    xs, ys = [p[0] for p in xy], [p[1] for p in xy]
    w, h = (max(xs) - min(xs)) or 1.0, (max(ys) - min(ys)) or 1.0
    s = min((box - 2 * margin) / w, (box - top - bottom - 2 * margin) / h)
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    mid_y = top + (box - top - bottom) / 2
    return [(box / 2 + s * (x - cx), mid_y - s * (y - cy)) for x, y in xy]


def ranks(values):
    """Rank equal-length segments: values are floats, equal within 1e-9 relative tolerance."""
    order = sorted(set(round(v, 9) for v in values))
    return [order.index(round(v, 9)) for v in values], order


def panel(xy, segments="all", title=None, box=220, margin=26, labels=None, marks=None,
          colour_of=None, radius=4.2, width=2.0, faint_segments=None, extra=None, subtitle=None,
          title_size=12.5, fill=INK):
    """One panel as an SVG fragment.

    xy            list of float points
    segments      "all", "shortest", None, or a list of (i, j) pairs to draw in colour
    colour_of     function (i, j) -> colour index; default: rank of the Euclidean length
    labels        optional list of labels for the points
    marks         {point index: fill colour} to highlight points
    faint_segments  list of (i, j) drawn thin and grey underneath
    extra         list of (x, y, kind) extra markers in the same coordinates; kind "ring" or "dot"
    """
    n = len(xy)
    all_pts = list(xy) + [(e[0], e[1]) for e in (extra or [])]
    fitted = _fit(all_pts, box, margin, top=(title_size + 3 if title else 0), bottom=(13 if subtitle else 0))
    pts, ext = fitted[:n], fitted[n:]
    dist = lambda i, j: math.dist(xy[i], xy[j])
    all_pairs = list(combinations(range(n), 2))
    rk, _ = ranks([dist(i, j) for i, j in all_pairs])
    rank_of = dict(zip(all_pairs, rk))
    if segments == "all":
        segs = all_pairs
    elif segments == "shortest":
        segs = [e for e in all_pairs if rank_of[e] == 0]
    else:
        segs = [tuple(sorted(e)) for e in (segments or [])]
    cidx = colour_of or (lambda i, j: rank_of[(i, j)])
    out = [f'<rect x="1" y="1" width="{box - 2}" height="{box - 2}" rx="10" fill="{SURFACE}" stroke="{FRAME}"/>']
    if title:
        out.append(f'<text x="{box / 2}" y="{title_size + 4.5:.1f}" text-anchor="middle" font-size="{title_size}" fill="{INK}" {FONT}>{escape(title)}</text>')
    if subtitle:
        out.append(f'<text x="{box / 2}" y="{box - 8}" text-anchor="middle" font-size="11" fill="{MUTED}" {FONT}>{escape(subtitle)}</text>')
    for i, j in (faint_segments or []):
        (a, b), (c, d) = pts[i], pts[j]
        out.append(f'<line x1="{a:.1f}" y1="{b:.1f}" x2="{c:.1f}" y2="{d:.1f}" stroke="{FAINT}" stroke-width="1"/>')
    for i, j in sorted(segs, key=lambda e: -cidx(*e)):          # long segments underneath
        (a, b), (c, d) = pts[i], pts[j]
        col = EDGE[cidx(i, j) % len(EDGE)]
        out.append(f'<line x1="{a:.1f}" y1="{b:.1f}" x2="{c:.1f}" y2="{d:.1f}" stroke="{col}" stroke-width="{width}" stroke-linecap="round"/>')
    for (x, y), e in zip(ext, extra or []):
        if e[2] == "ring":
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius + 1.5}" fill="{SURFACE}" stroke="{ACCENT}" stroke-width="1.8"/>')
        else:
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius}" fill="{ACCENT}" stroke="{SURFACE}" stroke-width="1.5"/>')
    for k, (x, y) in enumerate(pts):
        colour = (marks or {}).get(k, fill)
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius}" fill="{colour}" stroke="{SURFACE}" stroke-width="{min(1.5, radius / 2.5):.1f}"/>')
        if labels:
            out.append(f'<text x="{x + 7:.1f}" y="{y - 6:.1f}" font-size="12" fill="{MUTED}" {FONT}>{escape(str(labels[k]))}</text>')
    return "".join(out), box


def legend(entries, width, title=None):
    """A one-line legend: entries are (colour index or colour string, label)."""
    out, x = [], 6
    if title:
        out.append(f'<text x="{x}" y="15" font-size="12" fill="{MUTED}" {FONT}>{escape(title)}</text>')
        x += 7 * len(title) + 14
    for c, label in entries:
        col = EDGE[c % len(EDGE)] if isinstance(c, int) else c
        out.append(f'<line x1="{x}" y1="11" x2="{x + 18}" y2="11" stroke="{col}" stroke-width="3" stroke-linecap="round"/>')
        out.append(f'<text x="{x + 24}" y="15" font-size="12" fill="{INK}" {FONT}>{escape(str(label))}</text>')
        x += 24 + 7.2 * len(str(label)) + 18
    return "".join(out), 22


def figure(panels, columns=None, legend_entries=None, legend_title=None, gap=10):
    """Arrange panels (as returned by `panel`) in a grid, with an optional legend line."""
    columns = columns or len(panels)
    box = panels[0][1]
    rows = -(-len(panels) // columns)
    width = columns * box + (columns - 1) * gap
    body, y0 = [], 0
    if legend_entries:
        leg, h = legend(legend_entries, width, legend_title)
        body.append(leg)
        y0 = h + 4
    for k, (frag, _) in enumerate(panels):
        r, c = divmod(k, columns)
        body.append(f'<g transform="translate({c * (box + gap)},{y0 + r * (box + gap)})">{frag}</g>')
    height = y0 + rows * box + (rows - 1) * gap
    return SVG(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" '
               f'style="max-width:100%;height:auto;background:{SURFACE}" role="img">{"".join(body)}</svg>')


def polygon_xy(n, start=None):
    """n points evenly spaced on a circle: positions for drawing an abstract pattern."""
    start = 90 + 180 / n if start is None else start
    return [(math.cos(math.radians(start + 360 * k / n)), math.sin(math.radians(start + 360 * k / n))) for k in range(n)]


def abstract_panel(pattern, n, **kw):
    """A pattern drawn as a coloured complete graph on the corners of a regular n-gon."""
    idx = {e: c for e, c in zip(combinations(range(n), 2), pattern)}
    return panel(polygon_xy(n), segments="all", colour_of=lambda i, j: idx[(i, j)], **kw)


def lattice_xy(points):
    """Triangular-lattice coordinates (a, b) -> the plane, with a + b*omega, omega = e^(i pi/3)."""
    return [(a + b / 2, b * math.sqrt(3) / 2) for a, b in points]
