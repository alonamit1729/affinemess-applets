"""Numerical check of a family item: residual of the equal-distance equations at the sample, tangent-space
dimension (modulo similarity), and a trace along the family (predictor/corrector) recording the symmetry
type along the way.  Equations: |p_i - p_j|^2 - s_{cls(ij)} = 0 for every pair, in variables (coords, s).
Gauge: p0 = (0,0), p1 = (x1, 0), s_0 = 1 (fixes the similarity group) -> family dim = (#vars - 4) - rank."""
import itertools, math, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry import classify

class Family:
    def __init__(self, P, cls):
        P = np.asarray(P, float); self.n = n = len(P); self.E = list(itertools.combinations(range(n), 2))
        self.cls = list(cls); self.d = max(cls) + 1
        # normalise: p0 at origin, p1 on x-axis, |p0p1| = 1... keep general similarity gauge instead:
        # variables v = all 2n coords + d values s; gauge rows added to the Jacobian (translation x2, rotation, scaling)
        D2 = np.array([np.sum((P[i] - P[j]) ** 2) for i, j in self.E])
        s = np.array([D2[[k for k, (i, j) in enumerate(self.E) if cls[k] == c]].mean() for c in range(self.d)])
        self.v = np.concatenate([P.ravel(), s])
    def res(self, v):
        n = self.n; P = v[:2 * n].reshape(n, 2); s = v[2 * n:]
        return np.array([np.sum((P[i] - P[j]) ** 2) - s[self.cls[k]] for k, (i, j) in enumerate(self.E)])
    def jac(self, v):
        n = self.n; P = v[:2 * n].reshape(n, 2); J = np.zeros((len(self.E), len(v)))
        for k, (i, j) in enumerate(self.E):
            g = 2 * (P[i] - P[j]); J[k, 2 * i:2 * i + 2] = g; J[k, 2 * j:2 * j + 2] = -g; J[k, 2 * n + self.cls[k]] = -1
        return J
    def gauge(self, v):
        """rows spanning the tangent of the similarity orbit: translations, rotation, scaling"""
        n = self.n; P = v[:2 * n].reshape(n, 2); s = v[2 * n:]
        tx = np.zeros(len(v)); tx[0:2 * n:2] = 1
        ty = np.zeros(len(v)); ty[1:2 * n:2] = 1
        rot = np.zeros(len(v)); rot[0:2 * n:2] = -P[:, 1]; rot[1:2 * n:2] = P[:, 0]
        sc = np.zeros(len(v)); sc[:2 * n] = P.ravel(); sc[2 * n:] = 2 * s
        return np.array([tx, ty, rot, sc])
    def scale(self, v): return max(1.0, abs(v[2 * self.n:]).max())
    def tangent_dim(self, v, tol=1e-7):
        J = np.vstack([self.jac(v), self.gauge(v)])
        sv = np.linalg.svd(J, compute_uv=False); sv = sv / sv[0]
        rank = int((sv > tol).sum())
        return len(v) - rank, sv
    def project(self, v, it=40):
        for _ in range(it):
            r = self.res(v)
            if abs(r).max() < 1e-13 * self.scale(v): return v, True
            J = np.vstack([self.jac(v), self.gauge(v)]); rhs = np.concatenate([-r, np.zeros(4)])
            dv = np.linalg.lstsq(J, rhs, rcond=1e-9)[0]; v = v + dv
        return v, abs(self.res(v)).max() < 1e-10 * self.scale(v)
    def tangents(self, v, tol=1e-7):
        J = np.vstack([self.jac(v), self.gauge(v)]); u, sv, vt = np.linalg.svd(J); sv = sv / sv[0]
        return [vt[i] for i in range(len(v)) if (sv[i] if i < len(sv) else 0) < tol]
    def points(self, v): return v[:2 * self.n].reshape(self.n, 2)
    def margin(self, v):
        s = v[2 * self.n:]; sc = abs(s).max()
        return min(min(abs(a - b) for a, b in itertools.combinations(s, 2)), min(abs(s))) / sc
    def trace(self, v0, direction, steps=400, h=2e-3, rtol=1e-6):
        """walk along the family; returns list of (t, symmetry type, margin) and the final point"""
        v = v0.copy(); t = direction / np.linalg.norm(direction); out = []
        for k in range(steps):
            T = self.tangents(v)
            if not T: out.append((k, 'singular', None)); break
            tt = max(T, key=lambda u: abs(u @ t))
            if tt @ t < 0: tt = -tt
            t = tt
            w, ok = self.project(v + h * np.linalg.norm(v) * t)
            if not ok: out.append((k, 'lost', None)); break
            v = w
            if self.margin(v) < 1e-5: out.append((k, 'boundary', self.margin(v))); break
            if k % 20 == 0: out.append((k, classify(self.points(v), rtol)['type'], self.margin(v)))
        return out, v

def check_item(coords, cls, rtol=1e-6, trace=True):
    P = np.array([[float(x), float(y)] for x, y in coords]); F = Family(P, cls)
    v = F.v; r0 = abs(F.res(v)).max() / F.scale(v)
    v, ok = F.project(v); r1 = abs(F.res(v)).max() / F.scale(v)
    td, sv = F.tangent_dim(v)
    info = {'res0': r0, 'res_proj': r1, 'proj_ok': ok, 'tangent_dim': td, 'margin': F.margin(v),
            'sym0': classify(P, rtol)['type'], 'sym_loose': classify(P, 1e-3)['type'], 'n': F.n, 'd': F.d}
    if trace and td > 0:
        T = F.tangents(v); paths = []
        for u in (T[0], -T[0]):
            out, vend = F.trace(v, u)
            paths.append({'syms': [o[1] for o in out], 'margins': [o[2] for o in out], 'end': classify(F.points(vend), rtol)['type'], 'end_pts': F.points(vend).tolist()})
        info['trace'] = paths
    return info

if __name__ == '__main__':
    import json
    HERE = os.path.dirname(os.path.abspath(__file__))
    rows = json.load(open(os.path.join(HERE, 'step1_families.json')))
    which = sys.argv[1:] or ['e<=0']
    for r in rows:
        if which == ['e<=0'] and not (r['e'] <= 0 and r['order'] == 1): continue
        if which != ['e<=0'] and f"{r['cell']}:{r['idx']}" not in which: continue
        info = check_item(r['coords'], r['cls'])
        print(f"{r['cell']} #{r['idx']} dim={r['dim']} e={r['e']} tangent_dim={info['tangent_dim']} res={info['res0']:.1e}->{info['res_proj']:.1e} margin={info['margin']:.3g} sym={info['sym0']} loose={info['sym_loose']}")
        for p in info.get('trace', []):
            print('   trace:', p['syms'], 'end', p['end'])
