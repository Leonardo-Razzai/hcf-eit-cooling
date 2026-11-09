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
  exact     eliminated ME on (n_x, n_y), N = 6 per axis, against the chain at the same N; switching
            time tau; slow alternation needed, tau ~ 10/Gn  -> results/two_axes_exact.dat (~15 min)

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


def thermal():
    """Time to P00 >= 0.5, 0.8, 0.9 of the steady value from thermal n_x = n_y = nbar_0, both axes."""
    pc = E.paper_cases()
    rows = []
    for nw, bp, DG, s, lab in ((pc['A']['eta_sp']**-2, 1 / 2, 3.0, 1.0, 'A 1/2'), (pc['A']['eta_sp']**-2, 5 / 6, 3.0, 1.0, 'A 5/6'),
                               (3.0, 5 / 6, 3.0, 1.0, 'nu=3 5/6'), (pc['B']['eta_sp']**-2, 1 / 2, 3.0, 1.0, 'B 1/2')):
        for nb in (3.0, 10.0, 30.0):
            N = 160
            th = E.thermal(nb, N)
            P, tr = evolve_2d(nw, bp, DG, s, N, 0.5, 0.5, np.outer(th, th), 6000.0, out=5.0)
            P00 = tr[:, 3]
            end = P00[-1]
            ts = [tr[np.argmax(P00 >= f * end), 0] if (P00 >= f * end).any() else np.nan for f in (0.5, 0.8, 0.9)]
            rows.append([nw, bp, nb, end, tr[-1, 1]] + ts)
            log('thermal', lab, 'nbar0 %g' % nb, 'P00(end) %.4f nx(end) %.3f; t(0.5, 0.8, 0.9 P00) in 1/nu_r: %.0f %.0f %.0f' % (end, tr[-1, 1], *ts))
    save('two_axes_thermal.dat', ['nu_wrec', 'bp', 'nbar0', 'P00_end', 'nx_end', 't50', 't80', 't90'], rows)


def exact_2d(nw, bp, DG, s, N=6, tmax=600.0, dt=0.1, tau=32.0, only_x=False, recoil_y=True, quiet=False):
    """Eliminated master equation on (n_x, n_y): the x and y probes alternate, L_x for tau and
    L_y for tau, with 1/Gn << tau << 1/W; RK4 for the dense 2 N^2 x 2 N^2 density matrix.
    Averaging the two Liouvillians instead would add the probes coherently in H: a single
    probe at 45 degrees at half amplitude, which is not alternation (it does not cool)."""
    import scipy.sparse as sp
    from eitsim.operators import profile_ops
    fc = E.reduced_case(nw)
    p = E.params(fc, s=s, N=N, bp=bp, Delta_G=DG)
    b, w, V = oscillator(N)
    fp, fcp = profile_ops(p, w, V)
    fp, fcp = fp.toarray(), fcp.toarray()
    I = np.eye(N)
    n1 = np.diag(np.arange(N, dtype=float))
    Nn = np.kron(n1, I) + np.kron(I, n1)
    d = 2 * N * N
    Zero = np.zeros((N * N, N * N), complex)
    Hg = np.block([[Nn, Zero], [Zero, Nn - p['d2'] * np.eye(N * N)]])
    ck = 1.0 / (-p['D'] - 0.5j * p['G'])
    fc2 = np.kron(fcp, fcp)
    As = []
    for fp2 in (np.kron(fp, I), np.kron(I, fp)):
        As.append(np.hstack([0.5 * p['Wp'] * fp2, 0.5 * p['Wc'] * fc2]))        # N^2 x 2 N^2
    # emission along (u_x, u_y): dipole along z, sin^2(theta) = 4/5, four azimuths (u and -u are equivalent)
    us = [np.sqrt(0.8) * np.array([np.cos(ph), np.sin(ph)]) for ph in (np.pi / 8, 3 * np.pi / 8, 5 * np.pi / 8, 7 * np.pi / 8)]
    eta = p['eta_sp']
    stages = []
    Zr = np.zeros_like(As[0])
    for Am in As:
        H = Hg - ck.real * (Am.conj().T @ Am)
        Cs = []
        for ux, uy in us:
            Dx = (V * np.exp(1j * eta * ux * w)) @ V.conj().T
            Dy = (V * np.exp(1j * eta * uy * w)) @ V.conj().T if recoil_y else np.eye(N)
            DA = np.kron(Dx, Dy) @ Am
            for j, bj in ((0, bp), (1, 1 - bp)):
                if bj:
                    blk = [DA, Zr] if j == 0 else [Zr, DA]
                    Cs.append(np.sqrt(p['G'] * bj / 4) * ck * np.vstack(blk))
        G = -1j * H - 0.5 * sum(C.conj().T @ C for C in Cs)
        stages.append((G, G.conj().T, Cs, [C.conj().T for C in Cs]))

    def rhs(r, st):
        G, Gd, Cs, Cd = st
        out = G @ r + r @ Gd
        for C, Ch in zip(Cs, Cd):
            out += C @ r @ Ch
        return out
    th = E.thermal(0.3, N)
    r = np.zeros((d, d), complex)
    idx = np.arange(N * N)
    r[idx, idx] = np.kron(th, th)
    trace = []
    nsteps = int(tmax / dt)
    per = int(tau / dt)
    for k in range(nsteps):
        st = stages[0 if only_x else (k // per) % 2]            # L_x for tau, then L_y for tau
        k1 = rhs(r, st)
        k2 = rhs(r + 0.5 * dt * k1, st)
        k3 = rhs(r + 0.5 * dt * k2, st)
        k4 = rhs(r + dt * k3, st)
        r = r + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k % (2 * per) == 2 * per - 1:
            pop = np.real(np.diag(r))
            pp = (pop[:N * N] + pop[N * N:]).reshape(N, N)
            x, y = pp.sum(axis=1), pp.sum(axis=0)
            trace.append((pp[0, 0], 0.5 * (x[0] + y[0]), 0.5 * (np.arange(N) @ x + np.arange(N) @ y)))
            if not quiet:
                log('exact', 'nu %.2f bp %.2f t %.0f' % (nw, bp, (k + 1) * dt), 'P00 %.5f P0x %.5f nx %.5f trace %.6f' % (*trace[-1], pop.sum()))
            x0 = x[0]
    return (trace[-1], trace[-2], x0 if only_x else None, np.arange(N) @ x if only_x else None)


def exact():
    """Exact 2D master equation at N = 6 per axis against the chain at the same N: switching time
    tau at one point, then slow alternation at larger h. The one-axis limit of this code gives the
    1D master equation (0.97207 against 0.97211 in P0)."""
    rows = []
    B = E.paper_cases()['B']['eta_sp']**-2
    A = E.paper_cases()['A']['eta_sp']**-2
    N = 6
    runs = [(B, 1 / 2, tau, tm) for tau, tm in ((4.0, 300.0), (8.0, 300.0), (16.0, 400.0), (32.0, 600.0), (64.0, 800.0))]
    runs += [(5.0, 1 / 2, 32.0, 800.0), (A, 5 / 6, 32.0, 800.0), (3.0, 1 / 2, 32.0, 1000.0)]
    for nw, bp, tau, tm in runs:
        DG, s = 3.0, 1.0
        (P00, P0x, nx), prev, _, _ = exact_2d(nw, bp, DG, s, N=N, tmax=tm, tau=tau, quiet=True)
        p1 = chain_1d(nw, bp, DG, s, N)
        P2, tr = evolve_2d(nw, bp, DG, s, N, 0.5, 0.5, np.outer(p1, p1), 12000.0, tol=1e-4, out=100.0)
        x = P2.sum(axis=1)
        r1 = E.solve_banded(E.params(E.reduced_case(nw), s=s, N=N, bp=bp, Delta_G=DG), B=N)
        h = E.ALPHA / (bp * nw)
        rows.append([nw, bp, h, tau, N, P00, P0x, nx, abs(P0x - prev[1]), P2[0, 0], x[0], np.arange(N) @ x, r1['P0'], r1['n'], p1[0]])
        log('exact nu %.2f bp %.2f h %.3f tau %g:' % (nw, bp, h, tau), 'ME P00 %.4f P0x %.4f nx %.4f | chain P00 %.4f P0x %.4f nx %.4f | 1D ME P0 %.4f, chain %.4f'
            % (P00, P0x, nx, P2[0, 0], x[0], rows[-1][11], r1['P0'], p1[0]))
    save('two_axes_exact.dat', ['nu_wrec', 'bp', 'h', 'tau', 'N', 'P00_me', 'P0x_me', 'nx_me', 'dP0x_last', 'P00_chain', 'P0x_chain', 'nx_chain',
                                'P0_1D_me', 'n_1D_me', 'P0_1D_chain'], rows)


STUDIES = {'validate': validate, 'steady': steady, 'thermal': thermal, 'exact': exact}

if __name__ == '__main__':
    for name in sys.argv[1:]:
        STUDIES[name]()
