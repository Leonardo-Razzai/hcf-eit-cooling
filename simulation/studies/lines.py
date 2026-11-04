"""D1 versus D2: effective-Lambda comparison (stage 1).

  effective  every closed Lambda of 87Rb (line D1/D2, excited level F' = 1, 2,
             probe on ground level F = gp = 1, 2) as a three-level system:
             the line sets Gamma and eta_sp, the hyperfine structure sets the
             branching ratio bp (photons per cycle 1/bp) and the line
             strengths of the probe and coupling legs (laser powers).
             Master equation over s = 0.01 ... 10. "Wmax" is the largest rate
             among cold steady states (n_ss < N_COLD): with bp = 1/6 the
             recoil per removed phonon h = alpha eta_sp^2 / bp exceeds 1 and
             W only measures relaxation towards a hot state. The other excited
             levels enter only through the light-shift estimate
             atoms.neighbour_shift (stage 2 treats them exactly).
                                        -> results/lines_{A,B}.dat (summary)
                                           results/lines_scan_{A,B}.dat (W, n vs s)
                                           (~7 min)

Usage: python lines.py [effective]   (default: all)
"""
import sys

import numpy as np

from _common import E, NME, save, log

S_GRID = np.concatenate([[0.01], np.logspace(-1, 1, 9)])   # contains s = 1
N_COLD = 0.5
DELTA_G = 3.0


def configs():
    for line in ('D1', 'D2'):
        for Fe in E.excited_levels(line):
            for gp in E.GROUND_F:
                try:
                    yield E.lambda_config(line, Fe, gp)
                except ValueError:          # F' couples to only one ground level
                    pass


def effective():
    for key in 'AB':
        rows, scans = [], []
        for ic, cfg in enumerate(configs()):
            fc = E.paper_cases(cfg['line'])[key]
            om = fc['om']

            res = [E.solve_me(E.params(fc, s=s, N=NME, Delta_G=DELTA_G, bp=cfg['bp']))
                   for s in S_GRID]
            scans += [[ic, s, r['W'] * om, r['n']] for s, r in zip(S_GRID, res)]
            weak, one = res[0], res[int(np.argmin(abs(S_GRID - 1)))]
            cold = [(r['W'], s, r['n']) for s, r in zip(S_GRID, res) if r['n'] < N_COLD]
            Wmax, sopt, nopt = max(cold) if cold else (np.nan, np.nan, np.nan)
            p1 = E.params(fc, s=1.0, N=NME, Delta_G=DELTA_G, bp=cfg['bp'])
            sm = E.simple(p1, 1 / cfg['bp'])
            Pp = E.lp01_power(fc, p1['Wp'] * om, cfg['Sp'])
            Pc = E.lp01_power(fc, p1['Wc'] * om, cfg['Sc'])
            Delta = DELTA_G * fc['Gamma']
            shift = E.neighbour_shift(cfg['line'], cfg['Fe'], cfg['gc'], Delta)
            dmin = min(abs(Delta - 2 * np.pi * d) for f, d in cfg['others'].items()
                       if E.strength(cfg['line'], cfg['gc'], f) > 0) / fc['Gamma']
            h = (1 / cfg['bp']) * E.ALPHA * fc['eta_sp']**2
            rows.append([int(cfg['line'][1]), cfg['Fe'], cfg['gp'], cfg['bp'], cfg['Sp'],
                         cfg['Sc'], fc['Gamma'] / 2e6 / np.pi, fc['eta_sp'], h,
                         weak['n'], one['W'] * om, one['n'], Wmax * om, sopt, nopt,
                         sm['W'] * om, sm['n'], Pp * 1e9, Pc * 1e9, shift, dmin])
            log('effective', key, cfg['line'], 'F\'=%d gp=%d bp=%.3f' % (cfg['Fe'], cfg['gp'], cfg['bp']),
                'h %.2f  W(s=1) %.0f/s n %.3f  Wmax %.0f/s (s=%.2f)  Pp %.2f nW  shift %.3f'
                % (h, rows[-1][10], rows[-1][11], rows[-1][12], sopt, rows[-1][17], shift))
        save(f'lines_scan_{key}.dat', ['config', 's', 'W', 'n'], scans)
        save(f'lines_{key}.dat',
             ['line', 'Fe', 'gp', 'bp', 'Sp', 'Sc', 'GammaMHz', 'eta_sp', 'h',
              'nweak', 'W1', 'n1', 'Wmax', 'sopt', 'nopt', 'WSM1', 'nSM1',
              'Pp_nW', 'Pc_nW', 'nbshift', 'dmin_G'], rows)


STUDIES = {'effective': effective}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
