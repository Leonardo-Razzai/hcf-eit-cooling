"""Is EIT cooling on the D2 line worth a 2D study? One-axis Zeeman screening (all ground sublevels, all four
5P3/2 hyperfine levels, tensor light shift of the 1064 nm lattice), case A, against the D1 best.

Configurations with a dark state (F' <= F): D2 F'=1 with the probe on F=1 (b_p = 5/6, mirror of the D1 best; weak
coupling leg F=2 -> F'=1, S = 1/20) and D2 F'=2 with the probe on F=2 (b_p = 1/2, reference). F'=3 has no dark state;
F'=1 with the probe on F=2 has b_p = 1/6 (h > 1).

  screen   s = 1; polarisations (probe, coupling) = (pi, lin), (lin, pi) about the quantisation axis (= lattice
           polarisation, so the tensor shift is diagonal); Delta = 2, 3, 5 Gamma above the target level; two-photon
           detuning -1.5..1 nu_r; repumpers on both legs (complementary polarisation) gamma = 0.03, 0.1, 0.3; coupling
           scale 1, 2. N = 12, B = 2 (as the D1 Zeeman refine).            -> results/d2_screen.dat (4 workers)
  refine   best points of the screen: s = 0.3, 1, 3; convergence N = 28, B = 3; without the tensor shift; with the
           target level only                                                  -> results/d2_refine.dat
  gate     compare with the D1 one-axis Zeeman best (zeeman_refine.dat): continue to 2D only if a D2 point reaches
           P0 >= 0.95 with W0 >= 0.7 x the D1 best at the same P0 threshold.

Result (2026-11-13): GO. D2 F'=1 (probe lin on F=1, coupling pi on F=2, Delta = 2 Gamma, d = 0.15, gamma = 0.1,
s = 0.3): P0 = 0.951, W0 = 0.048 nu_r (2.2e3 s^-1), 1.76 x the D1 best at P0 >= 0.95; but at P0 >= 0.97 only 0.54 x
(max P0 0.972 at s = 0.1, gamma = 0.01). F'=2 (probe F=2) reaches P0 <= 0.75. The tensor light shift helps (P0 0.936
vs 0.919 without); the other 5P3/2 levels matter (target level only: 0.825). N = 12 -> 28: P0 -0.002, W0 -1 %.
2D: scheme lp11_d2 in eitsim/zeeman2d.py, cluster grid `zeeman2d.py list --grid d2`.

Usage: python d2_zeeman.py screen | refine | gate
"""
import itertools
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.constants as sc

from _common import E, RESULTS, save, log
from eitsim import zeeman as Z
from eitsim.atoms import tensor_shift

LINE = 'D2'
FES = [0, 1, 2, 3]
S2 = 2**-0.5
POL = {'pi': {0: 1.0}, 'lin': {1: S2, -1: S2}, 'sp': {1: 1.0}, 'sm': {-1: 1.0}}
COMP = {'pi': 'lin', 'lin': 'pi'}
CONFIGS = {0: (1, 1), 1: (2, 2)}          # code: (F', probe level)
PNAME = {'pi': 0, 'lin': 1}
U0_HZ = 530e-6 * sc.k / sc.h               # case A lattice depth (ground-state light shift), Hz
WORKERS = 4


def structure(Fe, gp, pp, pc, tensor=True, fes=None):
    st = Z.structure_lambda(LINE, Fe, gp, POL[pp], POL[pc], fes=fes or FES)
    if tensor:
        st['eps'] = st['eps'] + np.array([tensor_shift(LINE, f, m, U0_HZ) for f, m in st['excited']])
    return st


def repumpers(st, pp, pc, p, gamma):
    """Weak repumpers on both legs with the complementary polarisation, at the probe detuning (each excited level at
    its own detuning, see zeeman.ops), scattering rate gamma per sublevel on the target level."""
    out = []
    Wr = np.sqrt(4 * gamma * (p['D']**2 + p['G']**2 / 4) / p['G'])
    na, nE = len(st['ground']), len(st['excited'])
    for F, comp, i0 in ((st['gp'], COMP[pp], 0), (3 - st['gp'], COMP[pc], st['nP'])):
        g = Z.field_matrix(LINE, F, st['excited'], POL[comp]) / np.sqrt(Z._norm(LINE, F, st['Fe']))
        gm = np.zeros((na, nE), complex)
        gm[i0:i0 + g.shape[0]] = g
        out.append(dict(gm=gm, W=Wr, detuning=p['D']))
    return out


def point(a):
    cfg, pp, pc, DG, d2, gamma, wc, s, N, B, tensor, single = a
    Fe, gp = CONFIGS[cfg]
    fc = E.paper_cases(LINE)['A']
    p = Z.setup(fc, s, DG, d2, N, wc_scale=wc)
    st = structure(Fe, gp, pp, pc, tensor=tensor, fes=[Fe] if single else None)
    try:
        r = Z.solve(st, p, B=B, repump=repumpers(st, pp, pc, p, gamma))
    except Exception as exc:            # noqa: BLE001 (keep the scan going; record NaN)
        log('point failed', a, exc)
        r = dict(n=np.nan, P0=np.nan, W=np.nan, W0=np.nan)
    om = fc['om']
    return [cfg, PNAME[pp], PNAME[pc], DG, d2, gamma, wc, s, N, B, int(tensor), int(single), r['P0'], r['n'], r['W0'], r['W0'] * om]


COLS = ['cfg', 'pp', 'pc', 'DeltaG', 'd2', 'gamma', 'wc', 's', 'N', 'B', 'tensor', 'single', 'P0', 'n', 'W0', 'W0_per_s']


def screen():
    tasks = []
    for cfg, (pp, pc), DG, d2, gamma, wc in itertools.product(CONFIGS, (('pi', 'lin'), ('lin', 'pi')), (2.0, 3.0, 5.0),
                                                               np.round(np.linspace(-1.5, 1.0, 11), 3), (0.03, 0.1, 0.3), (1.0, 2.0)):
        tasks.append((cfg, pp, pc, DG, float(d2), gamma, wc, 1.0, 12, 2, True, False))
    log('screen: %d points' % len(tasks))
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(point, tasks, chunksize=4))
    save('d2_screen.dat', COLS, rows)
    best(np.array(rows))


def best(rows, label='screen'):
    rows = rows[np.isfinite(rows[:, 12])]
    for cfg in CONFIGS:
        if not np.any(rows[:, 0] == cfg):
            continue
        for thr in (0.95, 0.97):
            q = rows[(rows[:, 0] == cfg) & (rows[:, 12] >= thr)]
            if len(q):
                r = q[np.argmax(q[:, 14])]
                log("%s D2 F'=%d probe F=%d, P0 >= %.2f: best W0 %.4f nu_r (%.0f s^-1) at pol %d/%d, Delta %gG, d2 %+.2f, gamma %g, wc %g, s %g: P0 %.4f" % (
                    label, *CONFIGS[cfg], thr, r[14], r[15], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[12]))
            else:
                log("%s D2 F'=%d probe F=%d: no point with P0 >= %.2f (max P0 %.4f)" % (label, *CONFIGS[cfg], thr, rows[rows[:, 0] == cfg][:, 12].max()))


def refine():
    """D2 F'=1 (probe F=1) only (F'=2 reaches P0 <= 0.75 in the screen): the 3 fastest points with P0 >= 0.9 and the
    3 coldest points, at s = 0.1, 0.3, 1 and two-photon detuning +-0.15; weaker repumpers (0.01, 0.03) at s = 0.3
    for the coldest; convergence (N = 28, B = 3), no tensor shift and target level only at s = 1 for the coldest."""
    d = np.genfromtxt(RESULTS / 'd2_screen.dat', names=True)
    d = d[np.isfinite(d['P0']) & (d['cfg'] == 0)]
    inv = {v: k for k, v in PNAME.items()}
    key = lambda r: (0, inv[int(r['pp'])], inv[int(r['pc'])], float(r['DeltaG']), float(r['d2']), float(r['gamma']), float(r['wc']))
    fast = [key(r) for r in d[d['P0'] >= 0.9][np.argsort(-d[d['P0'] >= 0.9]['W0'])][:3]]
    cold = [key(r) for r in d[np.argsort(-d['P0'])][:3]]
    tasks = []
    for pk in dict.fromkeys(fast + cold):
        cfg, pp, pc, DG, d2, gamma, wc = pk
        for s in (0.1, 0.3, 1.0):
            for dd in (d2 - 0.15, d2, d2 + 0.15):
                tasks.append((cfg, pp, pc, DG, round(dd, 3), gamma, wc, s, 12, 2, True, False))
    for cfg, pp, pc, DG, d2, gamma, wc in cold:
        for g in (0.01, 0.03):
            for s in (0.1, 0.3):
                tasks.append((cfg, pp, pc, DG, d2, g, wc, s, 12, 2, True, False))
    cfg, pp, pc, DG, d2, gamma, wc = cold[0]
    tasks.append((cfg, pp, pc, DG, d2, gamma, wc, 1.0, 28, 3, True, False))       # convergence
    tasks.append((cfg, pp, pc, DG, d2, gamma, wc, 0.3, 28, 3, True, False))
    tasks.append((cfg, pp, pc, DG, d2, gamma, wc, 1.0, 12, 2, False, False))      # no tensor shift
    tasks.append((cfg, pp, pc, DG, d2, gamma, wc, 1.0, 12, 2, True, True))        # target level only
    tasks = list(dict.fromkeys(tasks))
    log('refine: %d points' % len(tasks))
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(point, tasks, chunksize=2))
    save('d2_refine.dat', COLS, rows)
    best(np.array([r for r in rows if r[8] == 12 and r[10] == 1 and r[11] == 0]), 'refine')


def gate():
    """D2 (screen + refine, N = 12, tensor shift, all levels) against the D1 one-axis Zeeman best."""
    z = np.genfromtxt(RESULTS / 'zeeman_refine.dat', names=True)
    z = z[np.isfinite(z['P0'])]
    d1 = {}
    for thr in (0.95, 0.97):
        q = z[(z['P0'] >= thr) & (np.abs(z['nu_wrec'] - 1.968) < 0.01)]
        d1[thr] = q['W0'].max() if len(q) else np.nan
    rows = []
    for f in ('d2_screen.dat', 'd2_refine.dat'):
        if (RESULTS / f).exists():
            x = np.genfromtxt(RESULTS / f, names=True)
            rows.append(x[(x['N'] == 12) & (x['tensor'] == 1) & (x['single'] == 0) & np.isfinite(x['P0'])])
    d2 = np.concatenate(rows)
    go = False
    for thr in (0.95, 0.97):
        for cfg in (0,):
            q = d2[(d2['cfg'] == cfg) & (d2['P0'] >= thr)]
            w = q['W0'].max() if len(q) else 0.0
            ratio = w / d1[thr]
            log("gate P0 >= %.2f: D1 best W0 %.4f nu_r; D2 F'=%d probe F=%d best W0 %.4f nu_r (ratio %.2f)" % (thr, d1[thr], *CONFIGS[cfg], w, ratio))
            if thr == 0.95 and ratio >= 0.7:
                go = True
    log('gate:', 'GO to 2D' if go else 'STOP: D2 not worth a 2D study')
    return go


if __name__ == '__main__':
    {'screen': screen, 'refine': refine, 'gate': gate}[sys.argv[1]]()
