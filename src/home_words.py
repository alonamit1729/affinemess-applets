#!/usr/bin/env python3
"""Write the Home tab's six face words, as Caveat outlines, into src/master.html.

The dodecahedron on the Home tab draws its words as outlines (so they follow the
faces in true perspective and no font file loads). This script regenerates the
one line between /*HOME-WORDS BEGIN*/ and /*HOME-WORDS END*/ in the master.

usage:  python3 src/home_words.py path/to/Caveat[wght].ttf
        (the variable Caveat from Google Fonts' repo:
         https://raw.githubusercontent.com/google/fonts/main/ofl/caveat/Caveat%5Bwght%5D.ttf
         Caveat is by Impallari Type, SIL Open Font License 1.1)
then:   python3 src/build_site.py

Format: {word: {w: width in em, c: [contour, ...]}}, y up, centred on the word,
0 at the middle of the lowercase letters. A contour is [x0, y0, op, ...] where
op 0 = line to (x, y) and op 1 = quadratic curve to (cx, cy, x, y).
"""
import json, os, re, sys
from fontTools.ttLib import TTFont
from fontTools.varLib.instancer import instantiateVariableFont
from fontTools.pens.basePen import BasePen

WORDS = ['Applets', 'Papers', 'Desmos', 'Media', 'Quora', 'About']   # keep in step with ORDER in AmHome
WEIGHT = 600
MASTER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'master.html')


class Pen(BasePen):
    def __init__(self, gs, dx):
        super().__init__(gs); self.dx = dx; self.cs = []
    def P(self, p): return [p[0] + self.dx, p[1]]
    def _moveTo(self, p): self.cur = self.P(p); self.cs.append(self.cur)
    def _lineTo(self, p): self.cur += [0] + self.P(p)
    def _qCurveToOne(self, c, p): self.cur += [1] + self.P(c) + self.P(p)
    def _curveToOne(self, c1, c2, p):   # not in a TrueType font; approximated if ever met
        self.cur += [1] + self.P(((c1[0] + c2[0]) / 2, (c1[1] + c2[1]) / 2)) + self.P(p)
    def _closePath(self): pass


def outlines(ttf):
    f = instantiateVariableFont(TTFont(ttf), {'wght': WEIGHT})
    gs, cmap, upm = f.getGlyphSet(), f.getBestCmap(), f['head'].unitsPerEm
    out = {}
    for w in WORDS:
        x, cs = 0, []
        for ch in w:
            g = cmap[ord(ch)]; pen = Pen(gs, x); gs[g].draw(pen); cs += pen.cs; x += gs[g].width
        W = x / upm
        def tr(c):
            o, i = [round(c[0] / upm - W / 2, 3), round(c[1] / upm - 0.25, 3)], 2
            while i < len(c):
                n = 2 if c[i] == 0 else 4; o.append(c[i])
                for k in range(0, n, 2):
                    o += [round(c[i + 1 + k] / upm - W / 2, 3), round(c[i + 2 + k] / upm - 0.25, 3)]
                i += 1 + n
            return o
        out[w] = {'w': round(W, 3), 'c': [tr(c) for c in cs]}
    s = json.dumps(out, separators=(',', ':'))
    return re.sub(r'(?<=[,\[:-])0\.', '.', s)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data = outlines(sys.argv[1])
    src = open(MASTER, encoding='utf-8').read()
    new, n = re.subn(r'(/\*HOME-WORDS BEGIN\*/).*?(/\*HOME-WORDS END\*/)',
                     lambda m: m.group(1) + 'var AM_HOME_WORDS = ' + data + ';' + m.group(2), src, flags=re.S)
    if n != 1:
        sys.exit('master.html needs exactly one /*HOME-WORDS BEGIN*/ ... /*HOME-WORDS END*/ line')
    open(MASTER, 'w', encoding='utf-8').write(new)
    print('wrote %d characters of word outlines into src/master.html' % len(data))


if __name__ == '__main__':
    main()
