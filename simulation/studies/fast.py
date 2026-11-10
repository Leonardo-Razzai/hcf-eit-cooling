"""Fast-channel checks (anomalous-cooling to-do): Fock-state drift at s = 1, case A, b_p = 1/2.
  sign    red detuning, Delta = -3 Gamma: the narrow peak sits at delta = -nu and a positive offset d2
          puts the red sideband on it; d2 chosen by the strongest cooling drift from n = 2, then the
          drift against Delta = +3 Gamma; gray-molasses picture: the resolved regime cools, the fast channel
          becomes heating (capture range R_sat(1 - c n / n_c))       -> results/fast_sign.dat
  delta   drift at Delta = 2 and 5 Gamma: the onset of the fast channel should move as
          2 (nu_r / Gamma_n)^2 / s                                    -> results/fast_delta.dat
"""
import sys
import time

import numpy as np

from _common import E, CASES, save, log

FC = CASES['A']
N0S = [1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 120, 160, 200]


def sign():
    """Red detuning. A steady state is not a good test: atoms that recoil past n_c run away to the
    truncation, so P0 is bimodal. Instead choose d2 by the strongest cooling drift from n = 2 (small N),
    then compare the Fock drift at Delta = -3 and +3 Gamma (N = 260)."""
    rows = []
    best = None
    for d2 in np.arange(0.6, 2.61, 0.2):
        dr = E.banded_fock_drift(E.params(FC, s=1.0, Delta_G=-3.0, d2=d2, N=60, bp=0.5), [2], B=4)[0, 0]
        log('sign d2 %+.2f: drift from n = 2: %.0f /s' % (d2, dr * FC['om']))
        rows.append([-3.0, d2, -2, np.nan, np.nan, np.nan, dr * FC['om']])
        if best is None or dr < best[1]:
            best = (d2, dr)
    d2b = best[0]
    for DG, d2 in ((-3.0, d2b), (3.0, 0.0)):
        t = time.time()
        out = E.banded_fock_drift(E.params(FC, s=1.0, Delta_G=DG, d2=d2, N=260, bp=0.5), N0S, B=4)
        for n0, (dr, dv) in zip(N0S, out):
            rows.append([DG, d2, n0, np.nan, np.nan, np.nan, dr * FC['om']])
        log('sign Delta %+g d2 %+.2f drift /s:' % (DG, d2), ' '.join('%.0f' % (x * FC['om']) for x in out[:, 0]), '%.0f s' % (time.time() - t))
    save('fast_sign.dat', ['DeltaG', 'd2', 'n0', 'P0', 'n', 'W0', 'drift'], rows)


def delta():
    rows = []
    for DG in (2.0, 5.0):
        t = time.time()
        p = E.params(FC, s=1.0, Delta_G=DG, N=260, bp=0.5)
        out = E.banded_fock_drift(p, N0S, B=4)
        nc = 2 / (p['Gn']**2) / 1.0
        for n0, (dr, dv) in zip(N0S, out):
            rows.append([DG, p['Gn'], nc, n0, dr * FC['om'], dv * FC['om']])
        log('delta %g Gamma (n_c = %.1f): drift /s' % (DG, nc), ' '.join('%.0f' % (x * FC['om']) for x in out[:, 0]), '%.0f s' % (time.time() - t))
    save('fast_delta.dat', ['DeltaG', 'Gn', 'nc', 'n0', 'drift', 'diffusion'], rows)


STUDIES = {'sign': sign, 'delta': delta}

if __name__ == '__main__':
    for name in sys.argv[1:]:
        STUDIES[name]()
