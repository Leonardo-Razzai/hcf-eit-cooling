"""Atom- and fibre-independent results of radial EIT cooling at a guided-mode node.

The cold-atom problem depends on nu_r / omega_rec (= 1/eta_sp^2), the branching
ratio bp of the Lambda and the laser settings (Delta/Gamma, s); Gamma/nu_r and
x0/a drop out (fibre.reduced_case). Rates and times are in units of nu_r and
1/nu_r. Delta = 3 Gamma throughout.

  steady   ME steady state and rate vs nu_r/omega_rec, bp, s
                                        -> results/universal_steady.dat   (~5 min)
  drift    ME phonon drift vs the calibrated GRM away from the calibration
           point (eta_sp^2 = 0.12-0.51)    -> results/universal_drift.dat    (~10 min)
  capture  GRM cooling from thermal states with initial nbar_0 (loading via a
           common pre-trap conserves nbar_0)
                                        -> results/universal_capture.dat  (~15 min)
  robust   fibre imperfections of paper Fig. 4 in reduced form: LP01 admixture
           eps (units of eta_perp; rate maximised over s), two-photon error and
           trap-frequency mismatch (s = 0.2)
                                        -> results/universal_robust.dat   (~5 min)

Usage: python universal.py [steady] [drift] [capture] [robust]   (default: all)
"""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.linalg import expm

from _common import E, NME, save, log
from temperature import first_crossing

WORKERS = 4
NU_WREC = np.logspace(np.log10(0.3), 2, 16)          # nu_r / omega_rec
BP = (1 / 6, 1 / 2, 5 / 6)
S_GRID = (0.03, 0.1, 0.3, 1.0, 2.0, 5.0, 10.0)


def case(nw):
    return E.reduced_case(nw)


# ---------------------------------------------------------------- steady state
def _steady_task(args):
    nw, bp = args
    rows = []
    for s in S_GRID:
        r = E.solve_me(E.params(case(nw), s=s, N=NME, bp=bp))
        rows.append([nw, bp, s, r['W'], r['n'], r['pgc'], r['pe']])
    log('steady', 'nu/wrec %.3g bp %.3f' % (nw, bp),
        'W(s=1) %.4f n %.3f' % tuple(rows[S_GRID.index(1.0)][3:5]))
    return rows


def steady():
    tasks = [(nw, bp) for nw in NU_WREC for bp in BP]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = [r for rs in ex.map(_steady_task, tasks) for r in rs]
    save('universal_steady.dat', ['nu_wrec', 'bp', 's', 'W', 'n', 'pgc', 'pe'], rows)


# ---------------------------------------------------------------- GRM validation
DRIFT_NW = (1.0, 3.0, 10.0, 30.0)
DRIFT_N0 = (3, 10)


def _drift_task(args):
    nw, bp = args
    me = E.fock_drift(E.params(case(nw), s=1.0, N=34, bp=bp), list(DRIFT_N0))
    NG = 80
    g = (np.arange(NG) @ E.grm_matrix(E.params(case(nw), s=1.0, N=NG, bp=bp), N=NG,
                                      csat=E.CSAT))[list(DRIFT_N0)]
    log('drift', 'nu/wrec %.3g bp %.3f' % (nw, bp), 'ME', np.round(me, 4), 'GRM', np.round(g, 4))
    return [[nw, bp, n0, m, x] for n0, m, x in zip(DRIFT_N0, me, g)]


def drift():
    tasks = [(nw, bp) for nw in DRIFT_NW for bp in (1 / 2, 5 / 6)]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = [r for rs in ex.map(_drift_task, tasks) for r in rs]
    save('universal_drift.dat', ['nu_wrec', 'bp', 'n0', 'ME', 'GRM'], rows)


# ---------------------------------------------------------------- capture from nbar_0
CAP_NW = (1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0)
NBAR0 = np.logspace(-1, np.log10(150), 16)
NCAP = 1600                       # Fock states (Rb case A trap depth, 1552 quanta)
P0_TARGET = 0.9


def _capture_task(args):
    nw, bp = args
    N = NCAP
    Q = E.grm_matrix(E.params(case(nw), s=1.0, N=N, bp=bp), N=N, csat=E.CSAT)
    ss = E.grm_steady_state(Q)
    nss, p0ss = np.arange(N) @ ss, ss[0]
    V = np.column_stack([E.thermal(x, N) for x in NBAR0])
    R0 = -(np.arange(N) @ Q @ V)                      # phonons per 1/nu_r
    p0init = V[0].copy()
    # fine steps first (fast cold states), then coarse steps (hot states)
    dt1, n1, every = 0.5, 400, 20
    U1 = expm(Q * dt1)
    U2 = np.linalg.matrix_power(U1, every)
    t, ts, nt, pt = 0.0, [], [], []
    target = p0ss - 0.1 * (p0ss - p0init)              # 90 % of the way to P0_ss
    # run until both criteria are met (P0 >= 0.9 only where the steady state allows it)
    stop = np.maximum(target, P0_TARGET) if p0ss > P0_TARGET + 0.01 else target
    k = 0
    while True:
        ts.append(t)
        nt.append(np.arange(N) @ V)
        pt.append(V[0].copy())
        if (pt[-1] >= stop).all() and t > 0 or t > 2e5:
            break
        if k < n1:
            V, t = U1 @ V, t + dt1
        else:
            V, t = U2 @ V, t + dt1 * every
        k += 1
    ts, nt, pt = np.array(ts), np.array(nt), np.array(pt)
    rows = []
    for j, nb in enumerate(NBAR0):
        t90 = first_crossing(ts, 1 - pt[:, j], 1 - P0_TARGET)
        trel = first_crossing(ts, target[j] - pt[:, j], 0.0)
        te = first_crossing(ts, nt[:, j] - nss, (nt[0, j] - nss) / np.e)
        rows.append([nw, bp, nb, nss, p0ss, t90, trel, te, R0[j]])
    log('capture', 'nu/wrec %.3g bp %.3f' % (nw, bp), 'nss %.3f P0ss %.3f' % (nss, p0ss),
        't90(nbar0=%.0f) %.0f /nu_r' % (NBAR0[-1], rows[-1][5]))
    return rows


def capture():
    tasks = [(nw, bp) for nw in CAP_NW for bp in (1 / 2, 5 / 6)]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = [r for rs in ex.map(_capture_task, tasks) for r in rs]
    save('universal_capture.dat',
         ['nu_wrec', 'bp', 'nbar0', 'nss', 'p0ss', 't90', 'trel', 'te', 'R0'], rows)


# ---------------------------------------------------------------- imperfections
ROB_NW = (1.5, 2.0, 3.0, 5.0, 10.0, 30.0)
EPS_ETA = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)                 # eps / eta_perp
ROB_S = (0.1, 0.3, 1.0, 3.0, 10.0)
D2_GRID = np.linspace(-0.6, 0.6, 13)
NU_REL = np.linspace(0.6, 1.6, 11)


def _robust_task(args):
    nw, bp = args
    c = case(nw)
    etap = E.gfac(1) * c['x0a']
    rows = []
    for e in EPS_ETA:                                    # kind 0: admixture, best s
        res = [E.solve_me(E.params(c, s=s, N=NME, bp=bp, eps=e * etap)) for s in ROB_S]
        W, s, n = max((r['W'], s, r['n']) for s, r in zip(ROB_S, res))
        rows.append([0, nw, bp, e, W, n, s])
    for d2 in D2_GRID:                                   # kind 1: two-photon error
        r = E.solve_me(E.params(c, s=0.2, N=NME, bp=bp, d2=d2))
        rows.append([1, nw, bp, d2, r['W'], r['n'], 0.2])
    for nr in NU_REL:                                    # kind 2: trap mismatch
        r = E.solve_me(E.params(c, s=0.2, N=NME, bp=bp, nu_rel=nr))
        rows.append([2, nw, bp, nr, r['W'] * nr, r['n'], 0.2])   # W in units of design nu_r
    log('robust', 'nu/wrec %.3g bp %.3f' % (nw, bp), 'W(eps) / W(0):',
        np.round([r[4] / rows[0][4] for r in rows[:len(EPS_ETA)]], 2))
    return rows


def robust():
    tasks = [(nw, bp) for nw in ROB_NW for bp in (1 / 2, 5 / 6)]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = [r for rs in ex.map(_robust_task, tasks) for r in rs]
    save('universal_robust.dat', ['kind', 'nu_wrec', 'bp', 'x', 'W', 'n', 's'], rows)


STUDIES = {'steady': steady, 'drift': drift, 'capture': capture, 'robust': robust}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
