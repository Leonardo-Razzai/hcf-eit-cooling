"""Cooling rate and steady state versus probe strength (paper Fig. 3).

Master equation: converged eliminated/banded solver E.steady (N and band from
numerics.toml); simple model with kappa = 1 and 2 photons per cycle.
Writes results/sat_A.dat (case B, converged at N = 20 with the full ME, is not rerun;
results/sat_B.dat comes from the previous version). Runtime ~20 min (case A,
s < 0.3 uses N = 400).
"""
import numpy as np

from _common import E, CASES, NME, save, log


def run():
    for key, fc in (('A', CASES['A']),):          # case B converged at N = 20: sat_B.dat kept
        rows = []
        for s in np.logspace(-2, 1, 13):
            p = E.params(fc, s=s, N=NME)
            r = E.steady(fc, s)
            s1, s2 = E.simple(p, 1), E.simple(p, 2)
            rows.append([s, r['n'], r['W'] * fc['om'], s1['n'], s1['W'] * fc['om'],
                         s2['n'], s2['W'] * fc['om'], r['pgc'], r['N'], r['band']])
            log('saturation', key, '%.3g' % s, 'n %.5f W %.1f (N %d, band %d, %s)'
                % (r['n'], r['W'] * fc['om'], r['N'], r['band'], r['status']))
        save(f'sat_{key}.dat', ['s', 'nME', 'WME', 'nSM1', 'WSM1', 'nSM2', 'WSM2',
                                'pgc', 'N', 'band'], rows)


if __name__ == '__main__':
    run()
