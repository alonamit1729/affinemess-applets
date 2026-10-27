"""m = 14, 15: X_D ∈ {R14, R15 - v, R15} (Erdos-Fishburn), all seven distances of X already occur in X_D, so every
extra point q has ALL its distances to X_D among the chord lengths of X_D.  Such a q has >= 2 distances in Delta(X_D),
hence is an intersection of two circles centred at vertices with chord-length radii: extend.candidates(X_D, 0)
enumerates them all (and discards vertices).  Also records the axis argument's prediction: each chord value at most
twice (a circle about q meets the circumcircle in <= 2 points), 14 vertices / 7 values -> q on a symmetry axis;
15 vertices -> impossible unless q = centre, and the circumradius is not a chord length."""
from extend import *
from cores11 import regular
import time
for name, P in [('R14', regular(14)), ('R15-v', regular(15)[1:]), ('R15', regular(15))]:
    t0 = time.time()
    vals, kA, surv = candidates(P, 0)
    R = mp.mpf(1)
    print(name, '|X_D| =', len(P), '|Delta| =', kA, ' points with all distances chord lengths:', len(surv),
          ' radius is a chord:', vals.find(R * R) is not None, ' (%.1fs)' % (time.time() - t0))
