"""Cooling rate and steady state versus probe strength (paper Fig. 3).

Writes results/sat_A.dat and results/sat_B.dat. Runtime ~1 min.
"""
import numpy as np

from _common import E, CASES, NME, save, log


def run():
    for key, fc in CASES.items():
        rows = []
        for s in np.logspace(-2, 1, 13):
            p = E.params(fc, s=s, N=NME)
            r = E.solve_me(p)
            s1, s2 = E.simple(p, 1), E.simple(p, 2)
            rows.append([s, r['n'], r['W'] * fc['om'], s1['n'], s1['W'] * fc['om'],
                         s2['n'], s2['W'] * fc['om'], r['pgc'], r['pe']])
            log('saturation', key, rows[-1][:3])
        save(f'sat_{key}.dat', ['s', 'nME', 'WME', 'nSM1', 'WSM1', 'nSM2', 'WSM2',
                                'pgc', 'pe'], rows)


if __name__ == '__main__':
    run()
