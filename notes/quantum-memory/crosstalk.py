#!/usr/bin/env python3
"""Crosstalk between spin-wave modes in a lattice column: numbers for memory_theory.md, Section 7.

Light: about a minute. Writes crosstalk_output.txt and crosstalk.png next to this file.
Model (memory_theory.md): collective operators b_k = sum_j u_k*(z_j) sigma_gs^j, a storage map
b_k(T) = sum_k' M_kk' b_k'(0) + F_k and a read map; crosstalk is the off-diagonal weight.
"""
import os

import numpy as np
from scipy.special import jv, jn_zeros

import memory_budget as mb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, 'crosstalk_output.txt'), 'w')
RNG = np.random.default_rng(1)


def say(*args):
    line = ' '.join(str(a) for a in args)
    print(line)
    OUT.write(line + '\n')


SITE = mb.LAM_L / 2


# ---------------------------------------------------------------- A. commutator (Gram) floor
def gram_floor():
    say('== A. Overlap of plane-wave spin waves on a randomly filled lattice')
    say('O_kk\' = [b_k, b_k\'^dag] = (1/N) sum_j exp(i (k\' - k) z_j); orthogonal grid Delta k = 2 pi m / L')
    say('prediction for m != 0: E|O|^2 = (1 - f) / N')
    say('%8s %6s %8s %12s %12s' % ('sites', 'f', 'N', '<|O|^2>', '(1-f)/N'))
    for ns, f in ((2000, 1.0), (2000, 0.5), (2000, 0.1), (20000, 0.5)):
        z = np.arange(ns) * SITE
        L = ns * SITE
        vals = []
        for _ in range(400):
            c = RNG.random(ns) < f
            N = c.sum()
            for m in (1, 3, 17):
                vals.append(abs(np.sum(np.exp(1j * 2 * np.pi * m / L * z[c])))**2 / N**2)
        N = f * ns
        say('%8d %6.1f %8.0f %12.2e %12.2e' % (ns, f, N, np.mean(vals), (1 - f) / N))
    say('Aliasing: k and k + 2 pi / d (d = %.0f nm site spacing, 2 pi/d = %.3e 1/m) are the same spin wave.' % (
        SITE * 1e9, 2 * np.pi / SITE))
    k795 = 2 * np.pi / mb.LAM_S
    for name, k in (('forward write, co-propagating', 2 * np.pi * mb.HFS / mb.C),
                    ('backward write, counter-propagating', 2 * k795)):
        r = (k * SITE / (2 * np.pi)) % 1
        say('  %-36s k = %.3e 1/m, k d / 2 pi mod 1 = %.3f' % (name, k, r))


# ---------------------------------------------------------------- B. independent-atom phases
def independent_dephasing():
    say('\n== B. Independent random phases phi_j: transfer into another orthogonal mode')
    say('prediction E|M_kk\'|^2 = |chi|^2 |O_kk\'|^2 + (1 - |chi|^2) / N, chi = <exp(-i phi)>')
    say('motional example: phi_j = 2 pi eps n_j T, n_j thermal-like with P0 per atom')
    ns = 4000
    z = np.arange(ns) * SITE
    L = ns * SITE
    u0 = np.exp(1j * 0 * z) / np.sqrt(ns)
    u1 = np.exp(1j * 2 * np.pi / L * z) / np.sqrt(ns)
    for P0, x in ((0.95, 0.5), (0.95, 3.0), (0.8, 3.0)):
        q = 1 - P0
        vals, diag = [], []
        for _ in range(300):
            n = RNG.geometric(1 - q, ns) - 1          # P(n) = (1-q) q^n
            ph = np.exp(-1j * x * n)                   # x = 2 pi eps T
            vals.append(abs(np.sum(np.conj(u1) * u0 * ph))**2)
            diag.append(abs(np.sum(np.conj(u0) * u0 * ph))**2)
        chi = np.mean(np.exp(-1j * x * (RNG.geometric(1 - q, 200000) - 1)))
        say('P0=%.2f  2 pi eps T=%.1f  |M_kk|^2=%.4f (|chi|^2=%.4f)  N|M_kk\'|^2=%.3f (pred %.3f)' % (
            P0, x, np.mean(diag), abs(chi)**2, ns * np.mean(vals), 1 - abs(chi)**2))


# ---------------------------------------------------------------- C. correlated phases
def correlated_phases():
    say('\n== C. Phase smooth along z: coherent transfer between spin waves')
    say('Periodic phase phi(z) = phi0 cos(q z): amplitude into k + m q is J_m(phi0) (Jacobi-Anger).')
    say('Source: LP02 admixture p in the 1064 nm lattice, phi0 = 2 pi delta_nu T, delta_nu = 2 sqrt(p g) ddiff u')
    g = (jv(1, mb.U01) / jv(1, mb.U02))**2
    ddiff, _ = mb.lattice_shift_scatter(mb.trap(1)['U'])
    say('%8s %6s %12s %16s %16s' % ('p', 'u', 'dnu (Hz)', 'T for J1^2=1e-3', 'T for J1^2=1e-2'))
    rows = []
    for p in (1e-3, 1e-4, 1e-5):
        for u in (1.0, 0.01):
            dnu = 2 * np.sqrt(p * g) * ddiff * u
            t3 = 2 * np.arcsin(np.sqrt(1e-3)) / (2 * np.pi * dnu)   # J1(x) ~ x/2
            t2 = 2 * np.sqrt(1e-2) / (2 * np.pi * dnu)
            say('%8.0e %6.2f %12.3f %16.3e %16.3e' % (p, u, dnu, t3, t2))
            rows.append((p, u, dnu))
    say('Linear phase gradient beta z shifts every k by beta: crosstalk to the channel at k + Delta k')
    say('peaks when beta T = Delta k (a gradient echo memory uses exactly this, reversibly).')
    return rows


# ---------------------------------------------------------------- D. read-out phase matching
def readout_crosstalk():
    say('\n== D. Read-out of a spin wave with a phase mismatch phi = Delta k L (forward retrieval kernel)')
    say('X(phi) = sigma_phi^dag K sigma_phi, sigma = optimal retrieval spin wave of the channel')
    phis = np.concatenate([np.linspace(0, 60, 301), np.linspace(61, 1500, 1440)])
    curves = {}
    for d in (15, 50, 150):
        z, w, K = mb.kernel(d)
        vals, vecs = np.linalg.eigh(K)
        v = vecs[:, -1]
        X = []
        for phi in phis:
            s = v * np.exp(1j * phi * z)
            X.append((s.conj() @ K @ s).real)
        X = np.array(X)
        curves[d] = X
        first = lambda lev: ('%.0f' % phis[np.argmax(X / X[0] < lev)]) if np.any(X / X[0] < lev) else '>1500'
        say('d=%-4d X(0)=%.4f; X/X0 at phi = 2 pi, 10 pi, d, 2d, 4d: %s; phi for X/X0 < 1e-2: %s, < 1e-3: %s' % (
            d, X[0], ' '.join('%.1e' % (np.interp(p, phis, X) / X[0])
                              for p in (2 * np.pi, 10 * np.pi, d, 2 * d, 4 * d)), first(1e-2), first(1e-3)))
    say('so channels must differ by Delta k L >> d (EIT momentum window ~ d / L, Gorshkov 2007), not just 2 pi')
    return phis, curves


def mode_pair_channels(a=15e-6, L=0.01, M=4):
    say('\n== E. Mode-pair channels (signal LP0m, control LP0m\'): wavevector of the spin wave')
    say('kappa_mm\' = beta_m - beta_m\' = (u_m\'^2 - u_m^2) / (2 k a^2); a = %.1f um, L = %.0f mm' % (a * 1e6, L * 1e3))
    k = 2 * np.pi / mb.LAM_S
    u = jn_zeros(0, M)
    beta = k - u**2 / (2 * k * a**2)
    chans = [(m, mp) for m in range(M) for mp in range(M) if m != mp]
    chans.append((0, 0))      # all m = m' pairs give kappa = 0: one channel
    kap = {c: beta[c[0]] - beta[c[1]] for c in chans}
    say('%-10s %12s' % ('(m, m\')', 'kappa L (rad)'))
    for c in chans:
        say('LP0%d/LP0%d %12.1f' % (c[0] + 1, c[1] + 1, kap[c] * L))
    # Reading with control m' converts a spin wave kappa into a field at beta_m' + kappa.
    # (i) channels (x, m') with the same control are phase matched into their own mode x: read together;
    # (ii) other channels: mismatch to the nearest guided mode decides the crosstalk.
    say('Control LP0m\' reads every channel (x, m\') at once, each into its own signal mode LP0x:')
    say('  random access is per control mode (groups of M channels), unless the outputs are sorted by mode.')
    worst = []
    for mp in range(M):
        for c2 in chans:
            if c2[1] == mp or c2 == (0, 0):
                continue
            mism = min((abs(beta[mp] + kap[c2] - beta[m2]) * L, m2) for m2 in range(M))
            worst.append((mism[0], mp, c2, mism[1]))
    worst.sort()
    say('other groups, smallest mismatch: ' + '; '.join(
        '%.0f rad (control LP0%d, channel LP0%d/LP0%d -> LP0%d)' % (w[0], w[1] + 1, w[2][0] + 1, w[2][1] + 1, w[3] + 1)
        for w in worst[:4]))
    say('The kappa = 0 spin wave (any m = m\') is read by every control: keep it empty.')
    return worst


def make_figure(phis, curves, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    oi = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
    for i, (d, X) in enumerate(curves.items()):
        ax[0].loglog(np.maximum(phis, 1e-9) / d, np.maximum(X / X[0], 1e-7), color=oi[i], lw=1.5, label='OD = %d' % (2 * d))
    ax[0].set_xlabel('Phase mismatch ΔkL / d')
    ax[0].set_xlim(0.02, 10)
    ax[0].set_ylabel('Read-out crosstalk X(ΔkL) / X(0)')
    ax[0].set_ylim(1e-6, 1.5)
    ax[0].legend(frameon=False, loc='lower left')
    ax[0].set_title('a  Phase-matching tail', loc='left', fontweight='bold', fontsize=9)
    T = np.geomspace(1e-3, 10, 20000)
    for i, (p, u, dnu) in enumerate([r for r in rows if r[1] == 0.01] + [r for r in rows if r[1] == 1.0][:1]):
        ax[1].loglog(T, jv(1, 2 * np.pi * dnu * T)**2, color=oi[i], lw=1.5,
                     label='p = %.0e, depth ×%g' % (p, u))
    ax[1].set_ylim(1e-6, 1)
    ax[1].set_xlabel('Storage time T (s)')
    ax[1].set_ylabel('Transfer to k ± q, J₁²(φ₀)')
    ax[1].legend(frameon=False, loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=7)
    ax[1].set_title('b  Lattice LP₀₂ beat, no echo', loc='left', fontweight='bold', fontsize=9)
    fig.savefig(os.path.join(HERE, 'crosstalk.png'), dpi=200)


if __name__ == '__main__':
    gram_floor()
    independent_dephasing()
    rows = correlated_phases()
    phis, curves = readout_crosstalk()
    mode_pair_channels()
    make_figure(phis, curves, rows)
    OUT.close()
