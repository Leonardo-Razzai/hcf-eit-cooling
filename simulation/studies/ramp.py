"""Optimal saturation ramp s(t) from thermal states (s-ramp step 3; theory in
notes/anomalous-cooling/s_ramp_theory.md). Case A, b_p = 1/2, Delta = 3 Gamma, one axis.

Subcommands:
  drift     master-equation drift from Fock states n0 = 1, 2, 3 at s = 3, 5, 10 (fills the gap
            below rates_hi.dat) -> results/ramp_drift_low.dat            (~20 min)
  chain     population chain Q(s): GRM Franck-Condon recoil kernel, each Fock column scaled so
            that its drift equals the master-equation drift R(n, s); validated against the
            stored constant-s master-equation runs from thermal states   -> results/ramp_chain_check.dat
  optimise[=smax]  minimum-time s(t) to P0 >= P0* by forward-backward (Pontryagin) sweeps on the chain,
            compared with the best constant s                           -> results/ramp_optimal.dat
            (optimise=3, optimise=5: s capped, nbar0 = 30 -> ramp_optimal_smax{3,5}.dat)
  twostage  s_max then s = 1 with the best switch time, on the chain    -> results/ramp_twostage.dat
  validate  step 4: the two-stage protocols in the banded master equation -> results/ramp_validate.dat
  tailcheck chain steady states at s = 0.1, 0.3, 1, 5 (minimal-model check) -> results/ramp_tail_chain.dat
"""
import sys
import time

import numpy as np
from scipy.linalg import expm

from _common import E, CASES, save, load, log

FC = CASES['A']
BP = 0.5
S_GRID = np.array([0.01, 0.1, 1.0, 3.0, 5.0, 10.0])
NCH = 400                        # chain Fock space (thermal tail of nbar0 = 30 negligible above)


# ---------------------------------------------------------------- 1. missing drift data
def drift():
    rows = []
    for s in (3.0, 5.0, 10.0):
        t = time.time()
        p = E.params(FC, s=s, N=260, bp=BP)
        out = E.banded_fock_drift(p, [1, 2, 3], B=4)
        for n0, (dr, dv) in zip((1, 2, 3), out):
            rows.append([s, 4, n0, dr * FC['om'], dv * FC['om']])
        log('drift s %g' % s, ' '.join('%.0f' % r[3] for r in rows[-3:]), '%.0f s' % (time.time() - t))
    save('ramp_drift_low.dat', ['s', 'band', 'n0', 'drift', 'diffusion'], rows)


# ---------------------------------------------------------------- 2. chain
def me_drift_table():
    """Master-equation drift -d<n>/dt (units of nu_r) on (S_GRID, n0)."""
    a = load('convergence/rates_vs_n.dat')
    b = load('convergence/rates_hi.dat')
    c = load('ramp_drift_low.dat')
    pts = {}
    for r in a:
        pts[(float(r['s']), int(r['n0']))] = -r['drift'] / FC['om']
    for d in (b, c):
        for r in d:
            if int(r['band']) == 4:
                pts[(float(r['s']), int(r['n0']))] = -r['drift'] / FC['om']
    return pts


def base_chain(s, N=NCH):
    """Unsaturated GRM chain at saturation s (Franck-Condon absorption, exact recoil kernel)."""
    p = E.params(FC, s=s, N=N, bp=BP)
    return E.grm_matrix(p, N=N, sat=False)


def chain_drift(Q):
    n = np.arange(Q.shape[0])
    return -(n @ Q)                      # -d<n>/dt from each Fock column


N_SWITCH = 5     # below: calibrated GRM saturation (validated for n <= 15, b_p = 1/2); from here: ME drift


def build_chain(s, pts, N=NCH):
    """Q(s) with column n scaled so that its drift equals the master-equation drift R(n, s) for
    n >= N_SWITCH; below, the calibrated GRM saturation 1/(1 + c_sat s_n). At low n the ME drift
    mixes removal with tail building (it even changes sign at n = 1 for s >= 5), so a drift match
    is not meaningful there.

    R is interpolated in log n at the two nearest grid values of s (log-log in s), and held
    constant beyond the largest measured n0. Column 0 (heating only) is left as in the GRM."""
    p = E.params(FC, s=s, N=N, bp=BP)
    Q0 = E.grm_matrix(p, N=N, sat=False)
    Qs = E.grm_matrix(p, N=N, sat=True)
    g_grm = np.diag(Qs) / np.diag(Q0)
    d0 = chain_drift(Q0)
    i = int(np.clip(np.searchsorted(S_GRID, s) - 1, 0, S_GRID.size - 2))
    s1, s2 = S_GRID[i], S_GRID[i + 1]
    w = np.clip(np.log(s / s1) / np.log(s2 / s1), 0, 1)
    n = np.arange(N)
    logR = []
    for sg in (s1, s2):
        ns = np.array(sorted(m for (ss, m) in pts if ss == sg and m >= N_SWITCH))
        R = np.array([pts[(sg, m)] for m in ns])
        logR.append(np.interp(np.log(np.maximum(n, 1)), np.log(ns), np.log(R)))
    R = np.exp((1 - w) * logR[0] + w * logR[1])
    g = g_grm.copy()
    g[N_SWITCH:] = R[N_SWITCH:] / d0[N_SWITCH:]
    return Q0 * g[None, :], g


def evolve(Q, P0vec, times):
    w, V = np.linalg.eig(Q)
    c = np.linalg.solve(V, P0vec)
    return np.real((V * c[None, :]) @ np.exp(np.outer(w, times)))


def chain():
    pts = me_drift_table()
    th = load('universal_thermal.dat')
    hs = load('universal_hotstage.dat')
    rows = []
    for s, src in ((0.3, 'thermal'), (1.0, 'thermal'), (3.0, 'hot'), (5.0, 'hot')):
        Q, g = build_chain(s, pts)
        if src == 'thermal':
            cfg = 2 if s == 0.3 else 0
            r = th[(th['config'] == cfg) & (th['nbar0'] == 30)]
        else:
            r = hs[(hs['bp'] == BP) & (hs['s'] == s) & (hs['N'] == 400)]
        t = r['t_s']
        P = evolve(Q, E.thermal(30.0, NCH), t * FC['om'])
        p0, nm = P[0], np.arange(NCH) @ P
        for k in range(t.size):
            rows.append([s, t[k], r['P0'][k], p0[k], r['n'][k], nm[k]])
        tme = [np.interp(f, r['P0'], t) * 1e3 for f in (0.5, 0.7, 0.8)]
        tch = [np.interp(f, p0, t) * 1e3 for f in (0.5, 0.7, 0.8)]
        log('chain s %g: t(P0=0.5,0.7,0.8) ME %s ms, chain %s ms; P0 end ME %.3f chain %.3f; column scale g(1,10,50) %s' % (
            s, ' '.join('%.2f' % x for x in tme), ' '.join('%.2f' % x for x in tch), r['P0'][-1], p0[-1],
            ' '.join('%.2f' % g[k] for k in (1, 10, 50))))
    save('ramp_chain_check.dat', ['s', 't_s', 'P0_me', 'P0_chain', 'n_me', 'n_chain'], rows)


# ---------------------------------------------------------------- 3. optimal control
S_CTRL = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])


class Propagator:
    """exp(Q dt) acting on vectors, from one eigendecomposition of Q (checked against expm)."""

    def __init__(self, Q):
        self.w, self.V = np.linalg.eig(Q)
        self.Vi = np.linalg.inv(self.V)
        dt = 0.3 * 1 / max(1e-9, np.max(np.abs(self.w.real)))
        x = np.zeros(Q.shape[0])
        x[30] = 1.0
        err = np.max(np.abs(self.fwd(x, 50 * dt) - expm(Q * 50 * dt) @ x))
        assert err < 1e-8, 'eigendecomposition inaccurate: %g' % err

    def fwd(self, x, dt):
        return np.real(self.V @ (np.exp(self.w * dt) * (self.Vi @ x)))

    def bwd(self, lam, dt):                    # exp(Q dt)^T lam
        return np.real(self.Vi.T @ (np.exp(self.w * dt) * (self.V.T @ lam)))


def optimise_fixed_T(props, Pi, T, K=30, sweeps=30, sctrl=S_CTRL):
    """Maximise P0(T) over piecewise-constant s on K log-spaced stages (forward-backward sweeps)."""
    edges = np.concatenate([[0.0], np.geomspace(T * 1e-3, T, K)])
    dts = np.diff(edges)
    ctrl = np.full(K, len(sctrl) - 1)           # start from s_max
    for sweep in range(sweeps):
        Ps = [Pi]
        for k in range(K):
            Ps.append(props[ctrl[k]].fwd(Ps[-1], dts[k]))
        lam = np.zeros(NCH)
        lam[0] = 1.0
        changed = 0
        for k in range(K - 1, -1, -1):          # backward: best s for stage k given later stages
            vals = [lam @ pr.fwd(Ps[k], dts[k]) for pr in props]
            jbest = int(np.argmax(vals))
            if jbest != ctrl[k]:
                changed += 1
                ctrl[k] = jbest
            lam = props[ctrl[k]].bwd(lam, dts[k])
        if changed == 0:
            break
    P = Pi
    for k in range(K):
        P = props[ctrl[k]].fwd(P, dts[k])
    return P[0], edges, sctrl[ctrl]


def optimise(smax=10.0):
    """Optimal ramps with s <= smax (smax < 10: only nbar0 = 30, written to ramp_optimal_smax*.dat)."""
    sctrl = S_CTRL[S_CTRL <= smax]
    pts = me_drift_table()
    Qs = [build_chain(s, pts)[0] for s in sctrl]
    props = [Propagator(Q) for Q in Qs]
    rows = []
    for nb0 in ((10.0, 30.0) if smax >= 10 else (30.0,)):
        Pi = E.thermal(nb0, NCH)
        for target in (0.5, 0.8, 0.9):
            # best constant s
            tc = []
            tt = np.geomspace(1e-5, 0.5, 600) * FC['om']
            for j, s in enumerate(sctrl):
                p0 = np.array([props[j].fwd(Pi, x)[0] for x in tt])
                ok = np.nonzero(p0 >= target)[0]
                tc.append(tt[ok[0]] / FC['om'] if ok.size else np.inf)
            jc = int(np.argmin(tc))
            # optimal ramp: bisect T
            lo, hi = 0.2 * tc[jc] * FC['om'], tc[jc] * FC['om']
            best = None
            for _ in range(11):
                mid = np.sqrt(lo * hi)
                p0T, edges, sk = optimise_fixed_T(props, Pi, mid, sctrl=sctrl)
                if p0T >= target:
                    hi, best = mid, (edges, sk)
                else:
                    lo = mid
            if best is None:
                best = optimise_fixed_T(props, Pi, hi, sctrl=sctrl)[1:]
            Topt = hi / FC['om']
            log('nbar0 %g P0*=%.2f: best constant s=%g, %.2f ms; ramp %.2f ms (gain %.2f); s(t) by stage: %s' % (
                nb0, target, sctrl[jc], tc[jc] * 1e3, Topt * 1e3, tc[jc] / Topt,
                ' '.join('%g' % x for x in best[1][::4])))
            for k in range(len(best[1])):
                rows.append([nb0, target, sctrl[jc], tc[jc], Topt, best[0][k] / FC['om'], best[0][k + 1] / FC['om'], best[1][k]])
    save('ramp_optimal.dat' if smax >= 10 else 'ramp_optimal_smax%g.dat' % smax, ['nbar0', 'P0target', 's_const', 't_const', 't_ramp', 't_start', 't_end', 's'], rows)


def twostage():
    """Two stages on the chain: s_max until t_s, then s = 1; best t_s, from nbar0 = 30."""
    pts = me_drift_table()
    props = {s: Propagator(build_chain(s, pts)[0]) for s in (1.0, 3.0, 5.0, 10.0)}
    Pi = E.thermal(30.0, NCH)
    tt = np.linspace(0, 15e-3, 1501)
    rows = []
    for target in (0.8, 0.9):
        for smax in (3.0, 5.0, 10.0):
            best = (np.inf, None)
            for ts in np.linspace(0.5e-3, 8e-3, 76):
                P = props[smax].fwd(Pi, ts * FC['om'])
                p0 = np.array([props[1.0].fwd(P, x * FC['om'])[0] for x in tt])
                ok = np.nonzero(p0 >= target)[0]
                if ok.size and ts + tt[ok[0]] < best[0]:
                    best = (ts + tt[ok[0]], ts)
            rows.append([target, smax, best[1], best[0]])
            log('two-stage P0*=%.1f: s=%g until %.2f ms, then s=1: T = %.2f ms' % (target, smax, best[1] * 1e3, best[0] * 1e3))
    save('ramp_twostage.dat', ['P0target', 'smax', 't_switch', 'T'], rows)


def validate():
    """s-ramp step 4: the two-stage protocols of twostage() in the banded master equation (N = 400,
    band 4, as universal.py hotstage), from nbar0 = 30, against the best constant s; target P0 = 0.9."""
    tw = load('ramp_twostage.dat')
    tw = tw[tw['P0target'] == 0.9]
    protocols = [('constant s = 2', [(2.0, 0.025)])]
    for smax in (10.0, 5.0, 3.0):
        ts = float(tw[tw['smax'] == smax]['t_switch'][0])
        protocols.append(('s = %g until %.2f ms, then s = 1' % (smax, ts * 1e3), [(smax, ts), (1.0, 0.025 - ts)]))
    rows = []
    for k, (name, st) in enumerate(protocols):
        t0 = time.time()
        stages = [(E.params(FC, s=s, N=400, bp=BP), dt * FC['om']) for s, dt in st]
        tt, n, P0 = E.banded_evolve_piecewise(stages, E.thermal(30.0, 400), B=4, h=0.5, n_out=400)
        t90 = np.interp(0.9, np.maximum.accumulate(P0), tt) / FC['om'] if P0.max() >= 0.9 else np.nan
        rows += [[k, x / FC['om'], y, z] for x, y, z in zip(tt, n, P0)]
        log('validate %-34s t(P0 = 0.9) = %.2f ms, P0 end %.4f, <n> end %.3f (%.0f s)' % (name, t90 * 1e3, P0[-1], n[-1], time.time() - t0))
    save('ramp_validate.dat', ['protocol', 't_s', 'n', 'P0'], rows)


def tailcheck():
    """Steady state of the chain (exact recoil kernel, master-equation removal rate) at s = 0.1, 0.3, 1, 5,
    for the comparison with the master-equation tails in notes/anomalous-cooling (minimal model)."""
    pts = me_drift_table()
    rows = []
    for s in (0.1, 0.3, 1.0, 5.0):
        Q, _ = build_chain(s, pts)
        w, V = np.linalg.eig(Q)
        ss = np.real(V[:, np.argmin(np.abs(w))])
        ss /= ss.sum()
        rows += [[s, n, ss[n]] for n in range(NCH)]
        log('tailcheck s %g: P0 %.4f' % (s, ss[0]))
    save('ramp_tail_chain.dat', ['s', 'n', 'P'], rows)


STUDIES = {'drift': drift, 'chain': chain, 'optimise': optimise, 'twostage': twostage, 'validate': validate, 'tailcheck': tailcheck}

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        name, _, val = arg.partition('=')
        STUDIES[name](*([float(val)] if val else []))
