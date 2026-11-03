"""Validation and calibration of the generalised rate model (GRM).

  steady   GRM (c_sat = 1) vs master-equation steady state  -> results/grmcheck.dat
  drift    master-equation phonon drift from Fock states vs GRM for several
           c_sat; prints the rms log error used to fix eitsim.rate_models.CSAT
                                                             -> results/drift_calib.dat
Runtime ~3 min. Usage: python calibration.py [steady] [drift]
"""
import sys

import numpy as np

from _common import E, CASES, NME, save, log

CSAT_GRID = (1, 1.5, 2, 2.25, 2.5, 3, 4)


def steady():
    rows = []
    for key, fc in CASES.items():
        for s in (0.01, 0.1, 0.3, 1.0):
            p = E.params(fc, s=s, N=NME)
            r = E.solve_me(p)
            ss = E.grm_steady_state(E.grm_matrix(p, N=60, csat=1.0))
            rows.append([0 if key == 'A' else 1, s, r['n'], np.arange(60) @ ss])
            log('steady', rows[-1])
    save('grmcheck.dat', ['case', 's', 'nME', 'nGRM'], rows)


def drift():
    n0s = [1, 3, 6, 10, 15]
    rows = []
    for key, fc in CASES.items():
        for s in (0.2, 1.0):
            me = E.fock_drift(E.params(fc, s=s, N=34), n0s)
            p = E.params(fc, s=s, N=60)
            g = {c: (np.arange(60) @ E.grm_matrix(p, N=60, csat=c))[n0s] for c in CSAT_GRID}
            for i, n0 in enumerate(n0s):
                rows.append([0 if key == 'A' else 1, s, n0, me[i] * fc['om']]
                            + [g[c][i] * fc['om'] for c in CSAT_GRID])
            log('drift', key, s, 'ME', np.round(me * fc['om']))
    save('drift_calib.dat', ['case', 's', 'n0', 'ME'] + ['c%g' % c for c in CSAT_GRID], rows)
    r = np.array(rows)
    for j, c in enumerate(CSAT_GRID):
        err = np.sqrt(np.mean(np.log(r[:, 4 + j] / r[:, 3])**2))
        log('c_sat = %g: rms log error %.3f' % (c, err))


STUDIES = {'steady': steady, 'drift': drift}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
