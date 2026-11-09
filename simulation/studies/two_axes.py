"""Both radial axes: rate chain on P(n_x, n_y), and a small exact check.

The chain has no free parameter beyond those fixed in 1D (units of nu_r), for a cycle of
kappa = 1/bp photons:
  red sideband     R(n) = s n Gn/2 / (1 + c s n), c = 1.125/bp (= c_sat 2.25 at bp = 1/2)
  blue sideband    B(n) = xi (n + 1) s Gn/2, xi = (Gamma/4 Delta)^2 = (Gn/4)^2
  fast channel     H(n) = c_h Rnet max(n - 0.8 n_c, 0)/n_c, n_c = 2/(Gn^2 s), c_h = 1.5, removal
                   without a photon (measured: more energy per cycle, not more cycles)
  recoil           T = bp P (1 - (1 - bp) P)^-1 per cycle on EACH axis, P the one-photon kernel
                   (axes independent: the correlation of u_x and u_y within a photon is neglected)
Time sharing: each axis is addressed for a fraction 1/2 of the time (alternating lobes, or fast
rotation). Only the addressed axis loses phonons, recoil acts on both.

  validate  1D chain against the master-equation map     -> results/two_axes_validate.dat
  steady    steady states, one axis and both, on the map -> results/two_axes_steady.dat (4 workers)
  thermal   time to P0 from thermal states, both axes     -> results/two_axes_thermal.dat
  exact     eliminated ME on (n_x, n_y), N = 8 per axis, against the chain
                                                         -> results/two_axes_exact.dat (~1 h)

Usage: python two_axes.py [validate] [steady] [thermal] [exact]
"""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import E, RESULTS, save, log
from eitsim.operators import oscillator

C_H = 1.5
WORKERS = 4


def kernel(nu_wrec, bp, N):
    """Recoil of one cycle on one axis, normalised column-wise."""
    b, w, V = oscillator(N)
    eta = nu_wrec**-0.5
    P = np.zeros((N, N))
    for u, wk in zip(E.NODES, E.WEIGHTS):
        Dk = (V * np.exp(1j * eta * u * w)) @ V.conj().T
        P += wk * np.abs(Dk)**2
    P /= P.sum(axis=0, keepdims=True)
    return bp * P @ np.linalg.inv(np.eye(N) - (1 - bp) * P)


def rates(nu_wrec, bp, DG, s, N):
    """Absorption matrix K[m, n], loss d(n), hill matrix Hm[m, n] (with its diagonal loss)."""
    Gn = 1.0 / DG
    h = E.ALPHA / (bp * nu_wrec)
    c = 1.125 / bp
    n = np.arange(N, dtype=float)
    R = s * n * Gn / 2 / (1 + c * s * n)
    B = (Gn / 4)**2 * (n + 1) * s * Gn / 2
    B[-1] = 0.0
    K = np.diag(R[1:], 1) + np.diag(B[:-1], -1)
    d = R + B
    nc = 2 / (Gn**2 * s)
    Rnet = max(1 - h, 0.0) * Gn / (2 * c)
    Hr = C_H * Rnet * np.maximum(n - 0.8 * nc, 0) / nc
    Hm = np.diag(Hr[1:], 1) - np.diag(Hr)
    return K, d, Hm


def chain_1d(nu_wrec, bp, DG, s, N=300):
    T = kernel(nu_wrec, bp, N)
    K, d, Hm = rates(nu_wrec, bp, DG, s, N)
    Q = T @ K - np.diag(d) + Hm
    A = Q.copy()
    A[0] = 1.0
    rhs = np.zeros(N)
    rhs[0] = 1.0
    p = np.linalg.solve(A, rhs)
    return p


def moments(p):
    n = np.arange(p.size)
    return p[0], p @ n, p @ n**2 - (p @ n)**2


def validate():
    D = np.genfromtxt(RESULTS / 'universal_map.dat', names=True)
    rows = []
    for r in D:
        if r['N'] == 250:                # keep the band-6 runs at N <= 120 and N = 60
            continue
        p = chain_1d(r['nu_wrec'], r['bp'], r['DeltaG'], r['s'], N=300)
        P0, m, _ = moments(p)
        rows.append([r['nu_wrec'], r['bp'], r['h'], r['DeltaG'], r['s'], r['P0'], P0, r['n'], m])
    rows = np.array(rows)
    save('two_axes_validate.dat', ['nu_wrec', 'bp', 'h', 'DeltaG', 's', 'P0_me', 'P0_chain', 'n_me', 'n_chain'], rows)
    for hmax in (0.5, 1.0, 9):
        k = rows[:, 2] < hmax
        e0 = (1 - rows[k, 6]) / (1 - rows[k, 5])
        en = rows[k, 8] / rows[k, 7]
        log('validate h < %g (%d points):' % (hmax, k.sum()), '(1-P0) chain/ME median %.2f, range %.2f-%.2f; <n> chain/ME median %.2f, range %.2f-%.2f'
            % (np.median(e0), e0.min(), e0.max(), np.median(en), en.min(), en.max()))


def evolve_2d(nu_wrec, bp, DG, s, N, fx, fy, P, tmax, tol=None, out=None):
    """Euler steps of dP/dt for matrix P[n_x, n_y]; returns P, t, and the trace of (t, <n_x>, <n_y>, P00)."""
    T = kernel(nu_wrec, bp, N)
    K, d, Hm = rates(nu_wrec, bp, DG, s, N)
    M = T @ K
    Tt, Mt, Hmt = T.T, M.T, Hm.T
    lam = (fx + fy) * (d.max() + abs(Hm).max())
    dt = 0.9 / lam
    n = np.arange(N, dtype=float)
    t, last = 0.0, None
    trace = []
    every = int(max(1, (out or tmax / 200) / dt))
    k = 0
    while t < tmax:
        dP = np.zeros_like(P)
        if fx:
            dP += fx * (M @ P @ Tt - d[:, None] * P + Hm @ P)
        if fy:
            dP += fy * (T @ P @ Mt - P * d[None, :] + P @ Hmt)
        P = P + dt * dP
        t += dt
        k += 1
        if k % every == 0:
            nx, ny = n @ P.sum(axis=1), n @ P.sum(axis=0)
            trace.append((t, nx, ny, P[0, 0], P[0].sum()))
            if tol and last is not None and abs(nx + ny - last[0]) < tol * (nx + ny) and abs(P[0, 0] - last[1]) < tol:
                break
            last = (nx + ny, P[0, 0])
    return P, np.array(trace)


def _steady_task(a):
    nw, bp, DG, s, N = a
    p = chain_1d(nw, bp, DG, s, N)
    P0_1, n_1, V_1 = moments(p)
    P = np.outer(p, p)
    P, tr = evolve_2d(nw, bp, DG, s, N, 0.5, 0.5, P, 12000.0, tol=2e-4, out=100.0)
    x, y = P.sum(axis=1), P.sum(axis=0)
    nn = np.arange(N)
    nx = nn @ x
    h = E.ALPHA / (bp * nw)
    # n at the truncation edge measures how much the tail reaches it
    edge = x[-20:].sum()
    log('steady nu/wrec %.2f bp %.2f h %.3f D %g s %g N %d' % (nw, bp, h, DG, s, N),
        '1D P0 %.4f n %.4f | 2D P00 %.4f P0x %.4f nx %.4f ratio %.1f edge %.1e t %.0f' % (P0_1, n_1, P[0, 0], x[0], nx, nx / (1 - x[0]), edge, tr[-1, 0]))
    return [nw, bp, h, DG, s, N, P0_1, n_1, P[0, 0], x[0], nx, edge, tr[-1, 0]]


def map_values():
    pc = E.paper_cases()
    return (1.4, pc['A']['eta_sp']**-2, 3.0, 5.0, pc['B']['eta_sp']**-2, 15.0)


def steady():
    tasks = [(nw, bp, DG, s, 120) for nw in map_values() for bp in (1 / 3, 1 / 2, 5 / 6)
             for DG in (3.0, 6.0) for s in (0.3, 1.0)]
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(_steady_task, tasks))
    save('two_axes_steady.dat', ['nu_wrec', 'bp', 'h', 'DeltaG', 's', 'N', 'P0_1D', 'n_1D', 'P00_2D', 'P0x_2D', 'nx_2D',
                                 'edge', 't_end'], rows)


STUDIES = {'validate': validate, 'steady': steady}

if __name__ == '__main__':
    for name in sys.argv[1:]:
        STUDIES[name]()
