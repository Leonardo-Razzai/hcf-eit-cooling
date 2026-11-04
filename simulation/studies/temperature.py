"""Cooling versus initial radial temperature (paper Fig. 6).

  scan   calibrated GRM (Delta = 3 Gamma, s = 1, LP11) from thermal states with
         k_B T_r / U0 = 1e-3 ... 1e-1 (20 log steps), both cases; U0 is the
         lattice depth for nu_z = 300 kHz. The range is limited by the
         harmonic radial trap (anharmonic shift ~ (3/8) k_B T / U0 < 4 %).
         The Fock space extends to the trap depth (U0 / hbar nu_r = 1552 for A,
         362 for B): recoil makes the energy diffuse by ~ eta_sp sqrt(n) per
         photon, so a smaller truncation removes the hot tail artificially.
                                    -> results/temperature_{A,B}.dat (metrics)
                                       results/temptrace_{A,B}.dat   (<n>(t), P0(t) at
                                       k_B T_r / U0 = TRACE_TU; column n_0p004 = 0.004)
  drift  master-equation phonon drift from Fock states n0 = 20, 30 at s = 1
         against the GRM, i.e. beyond the calibration range n0 <= 15
         (~25 min)                  -> results/drift_hi.dat

Usage: python temperature.py [scan] [drift]   (default: both)
"""
import sys

import numpy as np
from scipy.linalg import expm

from _common import E, CASES, save, log

TU = np.logspace(-3, -1, 20)                    # k_B T_r / U0 (metrics)
TRACE_TU = (1e-3, 4e-3, 1e-2, 4e-2, 1e-1)       # k_B T_r / U0 of the stored <n>(t), P0(t)
P0_TARGET = 0.9                                 # ground-state criterion
NFOCK = {'A': 1600, 'B': 400}
TGRID = {'A': (400e-3, 2000), 'B': (25e-3, 2000)}   # (t_max, steps)


def bose(fc, tu):
    """Thermal occupation at k_B T = tu * U0."""
    return 1 / np.expm1(E.HBAR * fc['om'] / (tu * fc['U0']))


def first_crossing(t, y, level):
    """First time y(t) <= level, linearly interpolated (0 if already below)."""
    i = int(np.argmax(y <= level))
    if y[i] > level:
        return np.nan
    if i == 0:
        return 0.0
    return t[i - 1] + (t[i] - t[i - 1]) * (y[i - 1] - level) / (y[i - 1] - y[i])


def scan():
    for key, fc in CASES.items():
        N = NFOCK[key]
        tmax, steps = TGRID[key]
        p = E.params(fc, s=1.0, N=N)
        Q = E.grm_matrix(p, N=N, csat=E.CSAT)
        ss = E.grm_steady_state(Q)
        nss, p0ss = np.arange(N) @ ss, ss[0]
        nb = bose(fc, TU)
        allTU = np.concatenate([TU, TRACE_TU])      # grid for metrics + trace temperatures
        V = np.column_stack([E.thermal(x, N) for x in bose(fc, allTU)])
        log('scan', key, 'n0 %.2f..%.1f' % (nb[0], nb[-1]),
            'truncation tail %.1e' % V[-20:].sum(axis=0).max())
        R0 = -(np.arange(N) @ Q @ V[:, :TU.size]) * fc['om']   # phonons / s
        t = np.linspace(0, tmax, steps + 1)
        U = expm(Q * (t[1] * fc['om']))
        nt = np.empty((t.size, allTU.size))
        pt = np.empty((t.size, allTU.size))
        for k in range(t.size):
            nt[k] = np.arange(N) @ V
            pt[k] = V[0]
            V = U @ V
        te = [first_crossing(t, nt[:, j] - nss, (nb[j] - nss) / np.e) for j in range(TU.size)]
        t90 = [first_crossing(t, 1 - pt[:, j], 1 - P0_TARGET) for j in range(TU.size)]
        rows = np.column_stack([TU, TU * fc['U0'] / E.KB * 1e6, nb, np.full(TU.size, nss),
                                np.full(TU.size, p0ss), np.array(te) * 1e3,
                                np.array(t90) * 1e3, R0, nt[-1, :TU.size], pt[-1, :TU.size]])
        save(f'temperature_{key}.dat',
             ['TU0', 'TuK', 'n0', 'nss', 'p0ss', 'te_ms', 't90_ms', 'R0', 'nend', 'p0end'], rows)
        tags = [('%g' % x).replace('.', 'p') for x in TRACE_TU]        # 0.004 -> 0p004
        save(f'temptrace_{key}.dat', ['tms'] + ['n_' + g for g in tags] + ['p_' + g for g in tags],
             np.column_stack([t * 1e3, nt[:, TU.size:], pt[:, TU.size:]]))
        for r in rows:
            log('  T/U0 %.2e  n0 %7.2f  t_e %8.3f ms  t_90 %8.3f ms  R0 %.3g /s  n_end %.3f  P0_end %.3f'
                % (r[0], r[2], r[5], r[6], r[7], r[8], r[9]))


def drift():
    rows = []
    for key, fc in CASES.items():
        for n0, N in ((20, 45), (30, 60)):
            me = E.fock_drift(E.params(fc, s=1.0, N=N), [n0])[0]
            NG = 2 * N
            g = (np.arange(NG) @ E.grm_matrix(E.params(fc, s=1.0, N=NG), N=NG, csat=E.CSAT))[n0]
            rows.append([0 if key == 'A' else 1, n0, me * fc['om'], g * fc['om']])
            log('drift', key, n0, 'ME %.0f' % rows[-1][2], 'GRM %.0f' % rows[-1][3])
    save('drift_hi.dat', ['case', 'n0', 'ME', 'GRM'], rows)


STUDIES = {'scan': scan, 'drift': drift}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
