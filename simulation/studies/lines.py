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

  multilevel every Lambda with all excited hyperfine levels in the master
             equation (stage 2), at Delta = 3 Gamma, s = 1, for both signs of
             the Raman interference through the other levels: lasers as in
             the three-level design, and with the two-photon detuning d2
             re-optimised (simple-model optimum refined with the ME).
                                        -> results/lines_ml_{A,B}.dat (~30 min, 4 workers)
  delta      the four cooling configurations (bp >= 1/2) versus Delta/Gamma,
             multi-level ME with re-optimised d2, sign +1; W is reliable only
             where nopt < N_RATE (column 'ok')
                                        -> results/lines_delta_{A,B}.dat (~40 min, 4 workers)

  imperfections  case A, D1 F'=2, D1 F'=1, D2 F'=2, D2 F'=1 (both signs), all
             excited levels, two-photon detuning at its stage-2 optimum
             (lines_ml_A): LP01 admixture eps (units of eta_perp, best of
             s = 0.3, 1, 3), two-photon error around the optimum and
             trap-frequency mismatch with fixed lasers (s = 1)
                                        -> results/lines_imperf_A.dat (~20 min, 4 workers)

Usage: python lines.py [effective] [multilevel] [delta] [imperfections]   (default: all)
"""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import minimize_scalar

from _common import E, NME, save, log, load

S_GRID = np.concatenate([[0.01], np.logspace(-1, 1, 9)])   # contains s = 1
N_COLD = 0.5
# W (relax_rate, from |g_p, 1>) is well defined only if n_ss << 1; the stage-2
# optimisation of the two-photon detuning selects among n_ss < N_RATE.
N_RATE = 0.2
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


def optimise_d2(fc, cfg, Delta_G=DELTA_G, s=1.0, sign=1.0, step=0.15):
    """Two-photon detuning (units of nu_r) that maximises the multi-level ME
    rate among states with n_ss < N_RATE: the multi-level simple model locates
    the red sideband on the narrow resonance (strongest local maximum of the
    absorption within +-10 nu_r), a
    parabola through three ME points refines it. Returns (d2, ME result at d2,
    simple-model d2)."""
    lv = E.excited_manifold(cfg['line'], cfg['Fe'], cfg['gp'])

    def par(d2):
        return E.params(fc, s=s, N=NME, Delta_G=Delta_G, d2=d2, levels=lv, sign=sign)
    grid = np.linspace(-10, 10, 801)
    A = E.absorption(par(0.0), 1.0 + grid)
    # the narrow resonance is the strongest interior local maximum; when the
    # neighbours wash it out, a broad background rises towards the window edge
    peaks = [i for i in range(1, grid.size - 1) if A[i - 1] < A[i] > A[i + 1]]
    g0 = grid[max(peaks, key=lambda i: A[i])] if peaks else 0.0
    d2sm = minimize_scalar(lambda g: -E.absorption(par(0.0), 1.0 + g),
                           bounds=(g0 - 0.05, g0 + 0.05), method='bounded').x
    pts = {g: E.solve_me(par(g)) for g in (d2sm - step, d2sm, d2sm + step)}
    xs = np.array(sorted(pts))
    ys = np.array([pts[g]['W'] for g in xs])
    c2, c1, _ = np.polyfit(xs, ys, 2)
    if c2 < 0 and xs[0] <= -c1 / (2 * c2) <= xs[-1]:
        g = -c1 / (2 * c2)
        pts[g] = E.solve_me(par(g))
    cold = [(r['W'], g) for g, r in pts.items() if r['n'] < N_RATE]
    g = max(cold)[1] if cold else d2sm
    return g, pts[g], d2sm


WORKERS = 4        # parallel ME solves (SuperLU is single-threaded; ~1 GB each)


def _multilevel_task(args):
    key, cfg, sign = args
    fc = E.paper_cases(cfg['line'])[key]
    om = fc['om']
    lv = E.excited_manifold(cfg['line'], cfg['Fe'], cfg['gp'])
    three = E.solve_me(E.params(fc, s=1.0, N=NME, Delta_G=DELTA_G, bp=cfg['bp']))
    raw = E.solve_me(E.params(fc, s=1.0, N=NME, Delta_G=DELTA_G, levels=lv, sign=sign))
    d2, opt, d2sm = optimise_d2(fc, cfg, sign=sign)
    popt = E.params(fc, s=1.0, N=NME, Delta_G=DELTA_G, d2=d2, levels=lv, sign=sign)
    beff = float(np.dot(E.excited_weights(popt, 1.0 + d2), [l[3] for l in popt['levels']]))
    sm = E.simple(popt, 1 / beff)
    row = [int(cfg['line'][1]), cfg['Fe'], cfg['gp'], sign, len(lv),
           three['W'] * om, three['n'], raw['W'] * om, raw['n'],
           d2sm, d2, opt['W'] * om, opt['n'], sm['W'] * om, sm['n'], beff]
    log('multilevel', key, cfg['line'], "F'=%d gp=%d sign %+d" % (cfg['Fe'], cfg['gp'], sign),
        '3-level W %.0f n %.3f | raw W %.0f n %.3f | d2 %.2f: W %.0f n %.3f'
        % tuple(row[5:9] + row[10:13]))
    return key, row


def multilevel():
    tasks = [(key, cfg, sign) for key in 'AB' for cfg in configs() for sign in (1.0, -1.0)]
    with ProcessPoolExecutor(WORKERS) as ex:
        res = list(ex.map(_multilevel_task, tasks))
    for key in 'AB':
        save(f'lines_ml_{key}.dat',
             ['line', 'Fe', 'gp', 'sign', 'K', 'W3', 'n3', 'Wraw', 'nraw', 'd2sm', 'd2opt',
              'Wopt', 'nopt', 'WSMopt', 'nSMopt', 'beff'], [r for k, r in res if k == key])


DELTAS = (1.0, 2.0, 3.0, 5.0, 8.0, 12.0)


def _delta_task(args):
    key, ic, cfg, DG = args
    fc = E.paper_cases(cfg['line'])[key]
    om = fc['om']
    three = E.solve_me(E.params(fc, s=1.0, N=NME, Delta_G=DG, bp=cfg['bp']))
    d2, opt, _ = optimise_d2(fc, cfg, Delta_G=DG)
    row = [ic, int(cfg['line'][1]), cfg['Fe'], cfg['gp'], DG,
           three['W'] * om, three['n'], d2, opt['W'] * om, opt['n'],
           int(three['n'] < N_RATE), int(opt['n'] < N_RATE)]
    log('delta', key, cfg['line'], "F'=%d gp=%d Delta %.0f G" % (cfg['Fe'], cfg['gp'], DG),
        '3-level W %.0f n %.3f | ML d2 %.2f W %.0f n %.3f' % tuple(row[5:10]))
    return key, row


def delta():
    cool = [c for c in configs() if c['bp'] >= 0.5]
    tasks = [(key, ic, cfg, DG) for key in 'AB' for ic, cfg in enumerate(cool) for DG in DELTAS]
    with ProcessPoolExecutor(WORKERS) as ex:
        res = list(ex.map(_delta_task, tasks))
    for key in 'AB':
        save(f'lines_delta_{key}.dat',
             ['config', 'line', 'Fe', 'gp', 'DeltaG', 'W3', 'n3', 'd2opt', 'Wopt', 'nopt',
              'ok3', 'okopt'],
             [r for k, r in res if k == key])


IMP_CFGS = (('D1', 2, 1), ('D1', 1, 2), ('D2', 2, 1), ('D2', 1, 1))
IMP_EPS = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)                 # eps / eta_perp
IMP_S = (0.3, 1.0, 3.0)
IMP_D2 = np.linspace(-0.6, 0.6, 9)                       # error around the optimum
IMP_NUREL = np.linspace(0.7, 1.3, 7)


def _imperf_task(args):
    line, Fe, gp, sign, d2opt, kind, x = args
    fc = E.paper_cases(line)['A']
    lv = E.excited_manifold(line, Fe, gp)

    def me(s, **kw):
        kw.setdefault('d2', d2opt)
        return E.solve_me(E.params(fc, s=s, N=NME, Delta_G=DELTA_G, levels=lv, sign=sign, **kw))
    if kind == 0:
        etap = E.gfac(1) * fc['x0a']
        W, s, n = max((r['W'], s, r['n']) for s, r in ((s, me(s, eps=x * etap)) for s in IMP_S))
    elif kind == 1:
        s = 1.0
        r = me(s, d2=d2opt + x)
        W, n = r['W'], r['n']
    else:                                   # lasers fixed: d2 in actual units is d2opt / nu_rel
        s = 1.0
        r = me(s, d2=d2opt / x, nu_rel=x)
        W, n = r['W'] * x, r['n']
    row = [int(line[1]), Fe, gp, sign, kind, x, W * fc['om'], n, s]
    log('imperfections', line, "F'=%d gp=%d sign %+d kind %d x %.2f" % (Fe, gp, sign, kind, x),
        'W %.0f n %.3f' % (row[6], n))
    return row


def imperfections():
    ml = load('lines_ml_A.dat')
    tasks = []
    for line, Fe, gp in IMP_CFGS:
        for sign in (1.0, -1.0):
            r = ml[(ml['line'] == int(line[1])) & (ml['Fe'] == Fe) & (ml['gp'] == gp)
                   & (ml['sign'] == sign)][0]
            for kind, xs in ((0, IMP_EPS), (1, IMP_D2), (2, IMP_NUREL)):
                tasks += [(line, Fe, gp, sign, float(r['d2opt']), kind, float(x)) for x in xs]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(_imperf_task, tasks))
    save('lines_imperf_A.dat', ['line', 'Fe', 'gp', 'sign', 'kind', 'x', 'W', 'n', 's'], rows)


STUDIES = {'effective': effective, 'multilevel': multilevel, 'delta': delta,
           'imperfections': imperfections}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
