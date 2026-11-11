"""Real fibre modes: vector structure of the LP11 group, longitudinal field at the node, beating
along z, and one vector mode for both radial axes.

  fields   E_z of TM01 / linear LP11 on the axis as a carrier eps_z (in the units of the purity
           study) against eta_perp; TE01-TM01-HE21 splitting and loss in the tube model of
           Zeisberger & Schmidt, Sci. Rep. 7, 11761 (2017), Eqs. (25), (26), against the
           antiresonance phase; beat lengths                       -> results/vector_fields.dat
  carrier  1D master equation with the on-axis E_z on the probe transition (worst case:
           quantisation axis not along the probe polarisation): eps = i eps_z (quadrature,
           the node is filled, not shifted), against eps = eps_z (real, shifted node) and 0
                                                                   -> results/vector_carrier.dat
  twoaxes  eliminated 2D master equation (N per axis) with several Lambda channels sharing g_p:
           'alt'   one channel, x and y probes alternating (reference, two_axes.switch)
           'coh'   one channel with profile (x + i y)/sqrt 2 (circular LP11 / OAM on one transition)
           'te01'  two channels, profile x on transition 1 and y on transition 2 (the two
                   transverse components of TE01 on different excited sublevels), shared g_c
           'te01s' the same with a separate g_c per channel
           P00(t) from |1,1>, steady P00 and W00 = P00ss / int (P00ss - P00) dt
                                                                   -> results/vector_twoaxes.dat
  shared   'te01' with the coupling of each channel scaled by 1/sqrt 2 (total light shift of g_c
           at the design value)                                    -> results/vector_shared.dat
  twoaxes_n  separate-g_c scheme at N = 8 per axis (convergence)   -> results/vector_twoaxes_n.dat

Usage: python vector_modes.py [fields] [carrier] [twoaxes]
"""
import sys

import numpy as np

from _common import E, RESULTS, save, log
from eitsim.operators import oscillator, profile_ops

EPS_SILICA = 1.45**2


# ---------------------------------------------------------------- fields
def eps_z(fc):
    """On-axis E_z of the linear LP11 (and TM01) in units of the LP01 axis field at equal power:
    E_z = (i / beta) div E_t, E_x = A J1(u1 r / a) cos(phi) -> eps_z = A u1 / (2 k a)."""
    k = 2 * np.pi / fc['lam']
    return E.amp_ratio(1) * E.U1[0] / (2 * k * fc['a'])


def tube_neff(kR, phi, mode, eps=EPS_SILICA):
    """Re(n_eff) - 1 and loss factor (alpha in units of 1/R) of the tube model, Eqs. (25), (26)."""
    from scipy.special import jn_zeros
    j = {'HE11': jn_zeros(0, 1)[0], 'TE01': jn_zeros(1, 1)[0], 'TM01': jn_zeros(1, 1)[0], 'HE21': jn_zeros(1, 1)[0]}[mode]
    fac = {'TE01': 1.0, 'TM01': eps, 'HE11': (eps + 1) / 2, 'HE21': (eps + 1) / 2}[mode]
    lfac = {'TE01': 1.0, 'TM01': eps**2, 'HE11': (eps**2 + 1) / 2, 'HE21': (eps**2 + 1) / 2}[mode]
    cot = 1 / np.tan(phi)
    dn = -j**2 / (2 * kR**2) - j**2 / kR**3 * cot / np.sqrt(eps - 1) * fac
    if mode == 'HE21':                    # sigma^4 term: HE21 - TE01 at antiresonance (c0 differs by -j^2 m s / 2 = +j^2)
        dn += j**2 / kR**4
    loss = (1 + cot**2) / (eps - 1) * j**3 / kR**3 * lfac        # alpha R
    return dn, loss


def fields():
    rows = []
    for tag, fc in E.paper_cases().items():
        ez = eps_z(fc)
        eta_perp = E.gfac(1) * fc['x0a']
        k = 2 * np.pi / fc['lam']
        kR = k * fc['a']
        log('fields case %s: eps_z %.4f, eta_perp %.4f, eps_z/eta_perp %.2f (eta_sp %.3f), kR %.1f' % (
            tag, ez, eta_perp, ez / eta_perp, fc['eta_sp'], kR))
        for phi_off in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4):   # phi = pi/2 + offset: antiresonance at 0; the expansion needs |cot phi| << kR
            phi = np.pi / 2 + phi_off * np.pi
            n = {m: tube_neff(kR, phi, m) for m in ('HE11', 'TE01', 'TM01', 'HE21')}
            d_tm = n['TM01'][0] - n['TE01'][0]
            d_he = n['HE21'][0] - n['TE01'][0]
            Ltm = fc['lam'] / abs(d_tm) if d_tm else np.inf
            Lhe = fc['lam'] / abs(d_he)
            loss_rel = [n[m][1] / n['HE11'][1] for m in ('TE01', 'TM01', 'HE21')]
            rows.append([0 if tag == 'A' else 1, fc['a'], kR, ez, eta_perp, phi_off, d_tm, d_he, Ltm, Lhe] + loss_rel)
            log('  phi - pi/2 = %.2f pi: dn(TM-TE) %.2e dn(HE21-TE) %.2e, L_beat TM-TE %.3g m, HE21-TE %.3g m; loss / HE11: TE %.2f TM %.2f HE21 %.2f' % (
                phi_off, d_tm, d_he, Ltm, Lhe, *loss_rel))
    save('vector_fields.dat', ['case', 'a', 'kR', 'eps_z', 'eta_perp', 'phi_off', 'dn_TM_TE', 'dn_HE21_TE',
                               'Lbeat_TM_TE', 'Lbeat_HE21_TE', 'loss_TE', 'loss_TM', 'loss_HE21'], rows)


# ---------------------------------------------------------------- carrier
def carrier():
    rows = []
    for tag, fc in E.paper_cases().items():
        ez = eps_z(fc)
        for s in (0.3, 1.0, 3.0):
            for kind, eps in (('none', 0.0), ('real', ez), ('quad', 1j * ez)):
                r = E.steady(fc, s, eps=eps)
                om = fc['om']
                rows.append([0 if tag == 'A' else 1, s, {'none': 0, 'real': 1, 'quad': 2}[kind], abs(eps),
                             r['P0'], r['n'], r['W'] * om, r['W0'] * om])
                log('carrier %s s %g %-4s eps %.4f: P0 %.4f n %.4f W0 %.0f s^-1' % (tag, s, kind, abs(eps), r['P0'], r['n'], r['W0'] * om))
    save('vector_carrier.dat', ['case', 's', 'kind', 'eps', 'P0', 'n', 'W', 'W0'], rows)


# ---------------------------------------------------------------- two axes
def evolve_2d(nw, bp, DG, s, channels, n_gc, N=6, tmax=900.0, dt=0.1, tau=None, init=(1, 1)):
    """Eliminated ME on (g_p, g_c1 .. g_c{n_gc}) x (n_x, n_y). channels: list of (profile, gc index),
    profile 'x', 'y' or 'xy' ((x + i y)/sqrt 2); each channel has its own excited level (incoherent
    with the others). tau: if given, the channels are switched on one at a time for tau each.
    Returns t, P00(t)."""
    fc = E.reduced_case(nw)
    p = E.params(fc, s=s, N=N, bp=bp, Delta_G=DG)
    b, w, V = oscillator(N)
    fp, fcp = profile_ops(p, w, V)
    fp, fcp = fp.toarray(), fcp.toarray()
    I = np.eye(N)
    n1 = np.diag(np.arange(N, dtype=float))
    M = N * N
    ng = 1 + n_gc
    d = ng * M
    Nn = np.kron(n1, I) + np.kron(I, n1)
    Hg = np.kron(np.eye(ng), Nn)
    for g in range(1, ng):
        Hg[g * M:(g + 1) * M, g * M:(g + 1) * M] -= p['d2'] * np.eye(M)
    ck = 1.0 / (-p['D'] - 0.5j * p['G'])
    prof = {'x': np.kron(fp, I), 'y': np.kron(I, fp)}
    prof['xy'] = (prof['x'] + 1j * prof['y']) / np.sqrt(2)
    fc2 = np.kron(fcp, fcp)
    us = [np.sqrt(0.8) * np.array([np.cos(ph), np.sin(ph)]) for ph in (np.pi / 8, 3 * np.pi / 8, 5 * np.pi / 8, 7 * np.pi / 8)]
    eta = p['eta_sp']
    Ds = []
    for ux, uy in us:
        Dx = (V * np.exp(1j * eta * ux * w)) @ V.conj().T
        Dy = (V * np.exp(1j * eta * uy * w)) @ V.conj().T
        Ds.append(np.kron(Dx, Dy))
    parts = []
    for ch in channels:
        pr, g = ch[:2]
        cs = ch[2] if len(ch) > 2 else 1.0                           # coupling Rabi scale of this channel
        A = np.zeros((M, d), complex)                                # e_k <- ground
        A[:, :M] = 0.5 * p['Wp'] * prof[pr]
        A[:, g * M:(g + 1) * M] = 0.5 * p['Wc'] * cs * fc2
        H = -ck.real * (A.conj().T @ A)
        Cs = []
        for Dk in Ds:
            DA = Dk @ A
            for tgt, bj in ((0, bp), (g, 1 - bp)):
                C = np.zeros((d, d), complex)
                C[tgt * M:(tgt + 1) * M, :] = DA
                Cs.append(np.sqrt(p['G'] * bj / 4) * ck * C)
        parts.append((H, Cs))

    def generator(idx):
        H = Hg + sum(parts[k][0] for k in idx)
        Cs = [C for k in idx for C in parts[k][1]]
        G = -1j * H - 0.5 * sum(C.conj().T @ C for C in Cs)
        return G, G.conj().T, Cs, [C.conj().T for C in Cs]

    stages = [generator(range(len(parts)))] if tau is None else [generator([k]) for k in range(len(parts))]

    def rhs(r, st):
        G, Gd, Cs, Cd = st
        out = G @ r + r @ Gd
        for C, Ch in zip(Cs, Cd):
            out += C @ r @ Ch
        return out
    r = np.zeros((d, d), complex)
    i0 = init[0] * N + init[1]
    r[i0, i0] = 1.0
    per = None if tau is None else int(round(tau / dt))
    nsteps = int(round(tmax / dt))
    P00 = []
    for k in range(nsteps):
        st = stages[0] if tau is None else stages[(k // per) % len(stages)]
        k1 = rhs(r, st)
        k2 = rhs(r + 0.5 * dt * k1, st)
        k3 = rhs(r + 0.5 * dt * k2, st)
        k4 = rhs(r + dt * k3, st)
        r = r + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        P00.append(sum(np.real(r[g * M, g * M]) for g in range(ng)))
    pop = np.real(np.diag(r)).reshape(ng, N, N).sum(axis=0)
    nx, ny = np.arange(N) @ pop.sum(axis=1), np.arange(N) @ pop.sum(axis=0)
    return np.array(P00), p, nx, ny


def twoaxes():
    rows = []
    B = E.paper_cases()['B']['eta_sp']**-2
    A = E.paper_cases()['A']['eta_sp']**-2
    N, DG, dt = 6, 3.0, 0.1
    schemes = {'alt': ([('x', 1), ('y', 1)], 1, 32.0), 'coh': ([('xy', 1)], 1, None),
               'te01': ([('x', 1), ('y', 1)], 1, None), 'te01s': ([('x', 1), ('y', 2)], 2, None)}
    for nw, bp in ((B, 0.5), (A, 5 / 6)):
        for s in (0.5, 1.0):
            p = E.params(E.reduced_case(nw), s=s, N=N, bp=bp, Delta_G=DG)
            r1 = E.solve_banded(p, B=N)
            for k, (name, (ch, ngc, tau)) in enumerate(schemes.items()):
                tmax = 1200.0
                tr, p2, nx, ny = evolve_2d(nw, bp, DG, s, ch, ngc, N=N, tmax=tmax, dt=dt, tau=tau)
                navg = int(round((64.0 if tau else 50.0) / dt))
                pss = tr[-navg:].mean()
                W00 = pss / (np.sum(pss - tr) * dt)
                drift = tr[-1] - tr[-navg - 1]
                rows.append([nw, bp, s, k, pss, W00, nx, ny, r1['W0'], r1['P0'], drift])
                log('twoaxes nu %.2f bp %.2f s %g %-5s: P00 %.4f W00 %.4f (W00/W0 %.3f) nx %.3f ny %.3f; 1D W0 %.4f P0 %.4f; drift %.1e' % (
                    nw, bp, s, name, pss, W00, W00 / r1['W0'], nx, ny, r1['W0'], r1['P0'], drift))
    save('vector_twoaxes.dat', ['nu_wrec', 'bp', 's', 'scheme', 'P00', 'W00', 'nx', 'ny', 'W0_1D', 'P0_1D', 'drift'], rows)


def shared():
    """Shared g_c with the coupling Rabi frequency of each channel scaled by 1/sqrt 2, so that the
    total light shift of g_c (sum over both excited levels) keeps the design value."""
    rows = []
    B = E.paper_cases()['B']['eta_sp']**-2
    A = E.paper_cases()['A']['eta_sp']**-2
    N, DG, dt = 6, 3.0, 0.1
    c = 1 / np.sqrt(2)
    for nw, bp in ((B, 0.5), (A, 5 / 6)):
        for s in (0.5, 1.0):
            p = E.params(E.reduced_case(nw), s=s, N=N, bp=bp, Delta_G=DG)
            r1 = E.solve_banded(p, B=N)
            tr, p2, nx, ny = evolve_2d(nw, bp, DG, s, [('x', 1, c), ('y', 1, c)], 1, N=N, tmax=1200.0, dt=dt)
            navg = int(round(50.0 / dt))
            pss = tr[-navg:].mean()
            W00 = pss / (np.sum(pss - tr) * dt)
            drift = tr[-1] - tr[-navg - 1]
            rows.append([nw, bp, s, pss, W00, nx, ny, r1['W0'], r1['P0'], drift])
            log('shared nu %.2f bp %.2f s %g: P00 %.4f W00 %.4f (W00/W0 %.3f) nx %.3f ny %.3f; drift %.1e' % (
                nw, bp, s, pss, W00, W00 / r1['W0'], nx, ny, drift))
    save('vector_shared.dat', ['nu_wrec', 'bp', 's', 'P00', 'W00', 'nx', 'ny', 'W0_1D', 'P0_1D', 'drift'], rows)


def twoaxes_n():
    """Convergence of the separate-g_c two-channel scheme in N: N = 8 per axis, s = 1, both cases."""
    rows = []
    B = E.paper_cases()['B']['eta_sp']**-2
    A = E.paper_cases()['A']['eta_sp']**-2
    DG, dt, s = 3.0, 0.1, 1.0
    for nw, bp in ((B, 0.5), (A, 5 / 6)):
        for N in (8,):
            p = E.params(E.reduced_case(nw), s=s, N=N, bp=bp, Delta_G=DG)
            r1 = E.solve_banded(p, B=N)
            tr, p2, nx, ny = evolve_2d(nw, bp, DG, s, [('x', 1), ('y', 2)], 2, N=N, tmax=1200.0, dt=dt)
            navg = int(round(50.0 / dt))
            pss = tr[-navg:].mean()
            W00 = pss / (np.sum(pss - tr) * dt)
            rows.append([nw, bp, s, N, pss, W00, nx, ny, r1['W0'], r1['P0']])
            log('twoaxes_n nu %.2f bp %.2f N %d: P00 %.4f W00 %.4f (W00/W0 %.3f) nx %.3f' % (nw, bp, N, pss, W00, W00 / r1['W0'], nx))
    save('vector_twoaxes_n.dat', ['nu_wrec', 'bp', 's', 'N', 'P00', 'W00', 'nx', 'ny', 'W0_1D', 'P0_1D'], rows)


STUDIES = {'fields': fields, 'carrier': carrier, 'twoaxes': twoaxes, 'shared': shared, 'twoaxes_n': twoaxes_n}

if __name__ == '__main__':
    for name in sys.argv[1:]:
        STUDIES[name]()
