"""Robustness of radial EIT cooling (paper Fig. 4).

Master equation: converged eliminated/banded solver E.steady (numerics.toml).
Runtime several hours for case A (weak-probe points use N = 400).

  carrier    LP01 admixture eps and trap offset x_s     -> results/carrier_A.dat   (~2 min)
  twophoton  two-photon detuning error                  -> results/twophoton_A.dat (<1 min)
  mismatch   actual vs design radial trap frequency     -> results/mismatch_{A,B}.dat (~1 min)

Usage: python robustness.py [carrier] [twophoton] [mismatch]   (default: all)
"""
import sys

import numpy as np

from _common import E, CASES, NME, save, log


def carrier():
    fc = CASES['A']
    svals = np.logspace(-2, 0.5, 8)
    rows = []
    for kind, vals in (('eps', [0, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2]),
                       ('xs', [0.5, 1, 2, 4, 7])):
        for v in vals:
            best = (-1, 0, 0)
            for s in svals:
                kw = {'eps': v} if kind == 'eps' else {'xs': v}
                p = E.params(fc, s=s, N=NME, **kw)
                r = E.steady(fc, s, **kw)
                if r['W'] > best[0]:
                    best = (r['W'], r['n'], s)
            epseff = v if kind == 'eps' else v * p['etap']
            pu = E.params(fc, s=1, N=NME, eps=epseff)
            rows.append([0 if kind == 'eps' else 1, v, epseff, best[0] * fc['om'],
                         best[1], best[2], E.unified(pu) * fc['om']])
            log('carrier', rows[-1])
    save('carrier_A.dat', ['kind', 'val', 'epseff', 'WmaxME', 'nAtWmax', 'sopt',
                           'Wunified'], rows)


def purity():
    """Steady state against LP01 admixture eps at fixed s (converged solver), cases A and B:
    P0, <n>, W and W0 for eps / eta_perp = 0 ... 5 at s = 0.3, 1, 3 -> purity.dat."""
    rows = []
    for key in ('A', 'B'):
        fc = CASES[key]
        etap = E.params(fc, s=1.0)['etap']
        for s in (0.3, 1.0, 3.0):
            for x in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
                r = E.steady(fc, s, eps=x * etap)
                rows.append([0 if key == 'A' else 1, s, x, x * etap, r['P0'], r['n'], r['W'] * fc['om'], r['W0'] * fc['om'], r['N'], r['band']])
                log('purity', key, 's %g eps/eta %g: P0 %.4f n %.4f W0 %.0f' % (s, x, r['P0'], r['n'], r['W0'] * fc['om']))
    save('purity.dat', ['case', 's', 'eps_eta', 'eps', 'P0', 'n', 'W', 'W0', 'N', 'band'], rows)


def twophoton():
    fc = CASES['A']
    rows = []
    for d2 in np.linspace(-0.6, 0.6, 13):
        row = [d2]
        for eps in (0.0, 0.05):
            r = E.steady(fc, 0.2, eps=eps, d2=d2)
            row += [r['n'], r['W'] * fc['om']]
        sm = E.simple(E.params(fc, s=0.2, N=NME, d2=d2), 2)
        row += [sm['n'], sm['W'] * fc['om']]
        rows.append(row)
        log('twophoton', row)
    save('twophoton_A.dat', ['d2', 'n0', 'W0', 'n5', 'W5', 'nSM', 'WSM'], rows)


def mismatch():
    for key, fc in (('A', CASES['A']),):          # case B converged at N = 20: kept as is
        rows = []
        for nr in np.linspace(0.6, 1.6, 21):
            p = E.params(fc, s=0.2, N=NME, nu_rel=nr)
            r = E.steady(fc, 0.2, nu_rel=nr)
            sm = E.simple(p, 2)
            # rates converted with the actual trap frequency
            rows.append([nr, r['n'], r['W'] * fc['om'] * nr, sm['n'],
                         sm['W'] * fc['om'] * nr])
            log('mismatch', key, rows[-1])
        save(f'mismatch_{key}.dat', ['nurel', 'nME', 'WME', 'nSM', 'WSM'], rows)


STUDIES = {'carrier': carrier, 'twophoton': twophoton, 'mismatch': mismatch, 'purity': purity}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
