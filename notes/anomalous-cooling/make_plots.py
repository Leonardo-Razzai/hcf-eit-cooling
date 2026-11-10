"""Figures for anomalous_cooling.tex (case A, 87Rb D1, Delta = 3 Gamma).

Uses eitsim only for parameters and the exact recoil kernel (no master-equation
solves) and the stored diagnostics in simulation/results/convergence/.
Writes figures/*.pdf. Usage (from this folder): python make_plots.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'simulation'))
sys.path.insert(0, str(ROOT / 'simulation' / 'plotting'))
import eitsim as E            # noqa: E402
import style                  # noqa: E402
from style import MM, WDOUBLE, WSINGLE, CA, CB, CG, CL, CO, letter   # noqa: E402

style.apply()
FIG = HERE / 'figures'
DATA = ROOT / 'simulation' / 'results' / 'convergence'
FC = E.paper_cases()['A']
P1 = E.params(FC, s=1.0)
H_A = (1 / P1['bp']) * E.ALPHA * FC['eta_sp']**2          # recoil quanta per removed phonon
SNC = 2 / P1['Gn']**2                                     # s n at loss of sideband resolution


def load(name):
    return np.genfromtxt(DATA / name, names=True)


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / name)
    plt.close(fig)


def fig_geometry():
    """LP11 node: field profiles, dark-state angle and the bright-state light shift
    (Sisyphus hill) compared with the trap and classical turning points."""
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    xa = np.linspace(-0.6, 0.6, 400)
    ax[0].plot(xa, jv(0, E.U01 * xa), 'k', label=r'LP$_{01}$ (coupling)')
    ax[0].plot(xa, E.amp_ratio(1) * jv(1, E.U1[0] * xa), color=CA, label=r'LP$_{11}$ (probe)')
    ax[0].plot(xa, E.gfac(1) * xa, '--', color=CA, lw=0.7, label=r'$G_{11}x/a$')
    for s, col in ((0.1, CL), (1.0, CO)):
        p = E.params(FC, s=s)
        th = np.arctan(p['Wp'] * E.amp_ratio(1) * jv(1, E.U1[0] * xa) / p['Wc'])
        ax[0].plot(xa, th, ':', color=col, lw=1.2, label=r'dark-state angle $\theta(x)$, $s=%g$' % s)
    ax[0].axvspan(-2 * np.sqrt(200) * FC['x0a'], 2 * np.sqrt(200) * FC['x0a'], color=CG, alpha=0.12, lw=0)
    ax[0].set_xlabel('$x/a$')
    ax[0].set_ylabel('field (LP$_{01}$ peak $=1$), angle (rad)')
    ax[0].set_ylim(-1.6, 1.6)
    ax[0].legend(loc='lower right', fontsize=5.5)
    ax[0].text(0.02, 0.97, 'shaded: radial excursion of $n=200$', transform=ax[0].transAxes, va='top')
    letter(ax[0], 'a')
    # bright-state light shift above its value on axis, in units of h nu_r, vs x / x0
    w = np.linspace(0, 35, 500)
    for s_, col in ((0.1, CL), (1.0, CO), (5.0, CB)):
        p = E.params(FC, s=s_)
        fp = E.amp_ratio(1) * jv(1, E.U1[0] * w * FC['x0a'])
        th = np.arctan(p['Wp'] * fp / p['Wc'])
        ax[1].plot(w, (p['Wp'] * fp)**2 / (4 * p['D']), color=col, label=r'$s=%g$' % s_)
        x45 = w[np.argmin(np.abs(th - np.pi / 4))]
        ax[1].plot(x45, (p['Wp'] * E.amp_ratio(1) * jv(1, E.U1[0] * x45 * FC['x0a']))**2 / (4 * p['D']),
                   'o', color=col, ms=3)
    for n in (1, 18, 180):
        A = np.sqrt(2 * (2 * n + 1))
        ax[1].axvline(A, color=CG, lw=0.5, ls=':')
        ax[1].text(A + 0.4, 0.97, '$n=%d$' % n, transform=ax[1].get_xaxis_transform(), va='top', color=CG)
    ax[1].set_xlim(0, 35)
    ax[1].set_ylim(0, 40)
    ax[1].set_xlabel(r'$x/x_0$')
    ax[1].set_ylabel(r'$\delta_B(x)-\delta_B(0)$ ($\hbar\nu_r$)')
    ax[1].legend(loc='center right', title='bright-state light shift')
    ax[1].text(0.02, 0.62, 'dots: $\\theta=\\pi/4$\ndotted: turning points', transform=ax[1].transAxes, va='top')
    letter(ax[1], 'b')
    save(fig, 'geometry.pdf')


def fig_recoil():
    """One photon from |n>: exact Franck-Condon distribution of the jump n -> n'
    (emission quadrature), compared with the single phonon removed per cycle."""
    N = 500
    from eitsim.operators import oscillator
    b, w, V = oscillator(N)
    P = np.zeros((N, N))
    for u, wk in zip(E.NODES, E.WEIGHTS):
        D = (V * np.exp(1j * FC['eta_sp'] * u * w)) @ V.conj().T
        P += wk * np.abs(D)**2
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 55 * MM), constrained_layout=True)
    for n, col in ((1, CA), (10, CL), (100, CO)):
        dn = np.arange(N) - n
        ok = np.abs(dn) <= 60
        ax[0].semilogy(dn[ok], np.clip(P[ok, n], 1e-12, None), 'o-', color=col, ms=1.5, lw=0.8,
                       label=r'$n=%d$' % n)
    ax[0].axvline(-1, color='k', lw=0.6, ls='--')
    ax[0].text(-58, 0.6, 'dashed: the single phonon\nremoved per cooling cycle', va='top')
    ax[0].set_xlim(-60, 60)
    ax[0].set_ylim(1e-8, 1.5)
    ax[0].set_xlabel(r"jump $n'-n$ after one spontaneous photon")
    ax[0].set_ylabel(r"$P(n'\,|\,n)$")
    ax[0].legend(loc='upper right')
    letter(ax[0], 'a')
    ns = np.arange(0, 300)
    k = np.arange(N)
    mean = np.array([(k - n) @ P[:, n] for n in ns])
    var = np.array([((k - n)**2) @ P[:, n] for n in ns]) - mean**2
    ax[1].plot(ns, np.sqrt(var), color=CA, label=r'spread $\sqrt{\mathrm{Var}(\Delta n)}$, one photon')
    ax[1].plot(ns, np.sqrt(2 * E.ALPHA * FC['eta_sp']**2 * (2 * ns + 1) / 2), '--', color=CG, lw=0.8,
               label=r'$\eta_{\rm sp}\sqrt{\alpha(2n+1)}$')
    ax[1].plot(ns, mean, color=CO, label=r'mean $\langle\Delta n\rangle=\alpha\eta_{\rm sp}^2$')
    ax[1].axhline(1, color='k', lw=0.6, ls=':')
    ax[1].set_xlabel('Fock state $n$')
    ax[1].set_ylabel(r'quanta (case A, $\eta_{\rm sp}=%.2f$)' % FC['eta_sp'])
    ax[1].legend(loc='upper left')
    letter(ax[1], 'b')
    save(fig, 'recoil.pdf')


def fig_regimes():
    d = load('rates_vs_n.dat')
    fig, ax = plt.subplots(figsize=(WSINGLE, 62 * MM), constrained_layout=True)
    for lo, hi, shade in ((1e-3, 1, '#EAF2FA'), (1, SNC, '#F4F4F4'), (SNC, 1e4, '#FDF0E6')):
        ax.axvspan(lo, hi, color=shade, lw=0)
    for x, lab in ((0.07, 'I: resolved\n$R\\propto n$'), (4, 'II: saturated\nplateau'),
                   (80, 'III: Sisyphus-\nlike')):
        ax.text(x, 1.4e5, lab, ha='center', va='top', fontsize=6)
    cols = {0.01: CG, 0.1: CL, 1.0: CA}
    for s, col in cols.items():
        r = d[d['s'] == s]
        ax.loglog(s * r['n0'], -r['drift'], 'o', color=col, ms=2.5, label='$s=%g$' % s)
    x = np.logspace(-2.5, 2.5, 200)
    Rsat, c = 2.0e3, E.CSAT                                 # c = calibrated GRM constant
    ax.loglog(x, Rsat * c * x / (1 + c * x), 'k-', lw=0.7, label=r'$R_{\rm sat}\,c\,sn/(1+c\,sn)$, $c=%.2f$' % E.CSAT)
    xs = x[x > SNC]
    ax.loglog(xs, Rsat * xs / SNC, 'k--', lw=0.7, label=r'$\propto sn$ (regime III)')
    ax.set_xlim(3e-3, 300)
    ax.set_ylim(20, 2e5)
    ax.set_xlabel(r'$s\,n_0$')
    ax.set_ylabel(r'removal rate $-d\langle n\rangle/dt$ (s$^{-1}$)')
    ax.legend(loc='lower right', fontsize=5.5)
    save(fig, 'regimes.pdf')


def fig_tail():
    t = load('tail_s0p1.dat')
    n, Pn = t['n'], t['P_N400']
    ok = np.isfinite(Pn) & (Pn > 0) & (n > 0)
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 58 * MM), constrained_layout=True)
    ax[0].loglog(n[ok], Pn[ok], color=CA, label=r'master equation, $s=0.1$, $N=400$')
    fit = (n >= 30) & (n <= 200) & ok
    mu, c0 = np.polyfit(np.log(n[fit]), np.log(Pn[fit]), 1)
    xx = np.logspace(np.log10(20), np.log10(300), 50)
    ax[0].loglog(xx, np.exp(c0) * xx**mu, '--', color=CO, label=r'fit $n^{%.2f}$ ($30\leq n\leq200$)' % mu)
    ref = np.exp(np.log(Pn[50]) + np.log(50) / H_A)
    ax[0].loglog(xx, ref * xx**(-1 / H_A), ':', color='k', label=r'$n^{-1/h}$, $h=%.2f$' % H_A)
    ax[0].axvline(1 / 0.1, color=CG, lw=0.5, ls=':')
    ax[0].axvline(SNC / 0.1, color=CG, lw=0.5, ls='--')
    ax[0].set_xlabel('Fock state $n$')
    ax[0].set_ylabel('$P(n)$, case A')
    ax[0].set_ylim(1e-9, 3e-2)
    ax[0].legend(loc='lower left', fontsize=5.5)
    ax[0].text(0.97, 0.97, 'dotted: $n=1/s$\ndashed: $n=2(\\nu_r/\\Gamma_n)^2/s$',
               transform=ax[0].transAxes, ha='right', va='top')
    letter(ax[0], 'a')
    h = np.linspace(0.02, 1.2, 400)
    ax[1].plot(h, 1 / h, 'k')
    for k, lab, col in ((1, 'normalisable', CB), (2, r'$\langle n\rangle$ finite', CO),
                        (3, r'Var$(n)$ finite', CA)):
        ax[1].axhline(k, color=col, lw=0.7, ls='--')
        ax[1].text(1.15, 11.5 - 0.9 * k, r'$\mu>%d$: ' % k + lab, color=col, ha='right', va='top')
    for hv, lab in ((H_A, 'A, $b_p=1/2$'), (H_A * 0.6, r'A, $b_p=5/6$'),
                    ((1 / 0.5) * E.ALPHA * E.paper_cases()['B']['eta_sp']**2, 'B, $b_p=1/2$')):
        ax[1].plot(hv, 1 / hv, 'o', color='k', ms=3)
        ax[1].annotate(lab, (hv, 1 / hv), xytext=(6, 4), textcoords='offset points')
    ax[1].set_xlim(0, 1.2)
    ax[1].set_ylim(0, 12)
    ax[1].set_xlabel(r'$h=\kappa\alpha\eta_{\rm sp}^2$ (recoil quanta per removed phonon)')
    ax[1].set_ylabel(r'tail exponent $\mu=1/h$, $P(n)\propto n^{-\mu}$')
    letter(ax[1], 'b')
    save(fig, 'tail.pdf')


def fig_phase():
    """h over (nu_r/omega_rec, b_p) with the thresholds, Rb cases, core-radius axis."""
    fig, ax = plt.subplots(figsize=(WSINGLE, 70 * MM), constrained_layout=True)
    x = np.logspace(np.log10(0.3), 2, 300)
    bp = np.linspace(0.1, 1.0, 300)
    X, Y = np.meshgrid(x, bp)
    Hm = E.ALPHA / (Y * X)
    cs = ax.contourf(X, Y, Hm, levels=[0, 1 / 3, 1 / 2, 1, 100],
                     colors=['#EAF2FA', '#D6E6F5', '#FDF0E6', '#F6D5C8'])
    ax.contour(X, Y, Hm, levels=[1 / 3, 1 / 2, 1], colors='k', linewidths=0.6)
    for xt, yt, lab in ((30, 0.55, '$h<1/3$\nmean and\nvariance finite'), (0.75, 0.5, '$1/2<h<1$\n$\\langle n\\rangle=\\infty$'), (0.42, 0.2, '$h>1$\nno cold\nstate')):
        ax.text(xt, yt, lab, ha='center', va='center', fontsize=6)
    cases = E.paper_cases()
    for xv, yv, lab in ((cases['A']['eta_sp']**-2, 0.5, 'A'), (cases['A']['eta_sp']**-2, 5 / 6, r'A, D1 $F^\prime=1$'),
                        (cases['B']['eta_sp']**-2, 0.5, 'B')):
        ax.plot(xv, yv, 'ko', ms=3)
        ax.annotate(lab, (xv, yv), xytext=(4, 3), textcoords='offset points')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu_r/\omega_{\rm rec}$')
    ax.set_ylabel(r'branching ratio $b_p$ ($\kappa=1/b_p$)')
    kl = 2 * np.pi / 1064e-9
    wrec = 2 * np.pi * 3632.5
    def nw_to_a(nw):
        return np.sqrt(2) * 2 * np.pi * 300e3 / (kl * 0.67 * np.asarray(nw) * wrec) * 1e6
    top = ax.secondary_xaxis('top', functions=(nw_to_a, lambda a: np.sqrt(2) * 2 * np.pi * 300e3
                                               / (kl * 0.67 * np.asarray(a) * 1e-6 * wrec)))
    top.set_xlabel(r'core radius $a$ ($\mu$m) for $^{87}$Rb D1, $\nu_z=300$ kHz')
    save(fig, 'phase.pdf')


def fig_scales():
    s = np.logspace(-2.5, 1, 200)
    fig, ax = plt.subplots(figsize=(WSINGLE, 62 * MM), constrained_layout=True)
    ax.fill_between(s, 1 / s, SNC / s, color='#F4F4F4', lw=0, label='power-law window (regime II)')
    ax.loglog(s, 1 / s, color=CA, label=r'$n^*=1/s$ (saturation)')
    ax.loglog(s, SNC / s, color=CO, label=r'$n_{\rm res}=2(\nu_r/\Gamma_n)^2/s$')
    nlin = 1 / (4 * (E.U1[0] * FC['x0a'])**2)
    ndep = FC['U0'] / (E.HBAR * FC['om'])
    nanh = (8 / 3) * P1['Gn'] * ndep
    for v, lab, ls in ((nlin, 'Bessel non-linearity', '--'), (nanh, r'anharmonic shift $=\Gamma_n$', '-.'),
                       (ndep, 'trap depth', ':')):
        ax.axhline(v, color='k', lw=0.6, ls=ls)
        ax.text(9, v * (1.15 if 'depth' in lab else 0.87), '%s ($n\\approx%.0f$)' % (lab, v), fontsize=5.5,
                ha='right', va='bottom' if 'depth' in lab else 'top')
    ax.axvline(0.1, color=CG, lw=0.5, ls=':')
    ax.set_xlim(3e-3, 10)
    ax.set_ylim(0.5, 5e4)
    ax.set_xlabel('saturation parameter $s$')
    ax.set_ylabel('Fock state $n$ (case A)')
    ax.legend(loc='lower left', fontsize=5.5)
    save(fig, 'scales.pdf')
    return nlin, nanh, ndep


def fp_from_measured(s):
    """Fokker-Planck stationary P(n) built only from the measured drift v and
    D = (1/2) dVar/dt: P = D^-1 exp(int v/D dn). Returns n0, r = -v n/D, n grid, P."""
    d = load('rates_vs_n.dat')
    d = d[d['s'] == s]
    n0, v, D = d['n0'], d['drift'], 0.5 * d['diffusion']
    r = -v * n0 / D
    n = np.arange(1, int(n0[-1]) + 1)
    rn = np.interp(np.log(n), np.log(n0), r)
    Dn = np.exp(np.interp(np.log(n), np.log(n0), np.log(D)))
    integ = np.concatenate([[0], np.cumsum(0.5 * (rn[1:] + rn[:-1]) * np.diff(np.log(n)))])
    return n0, r, n, np.exp(-integ) / Dn


def fig_fom():
    """Stored data against the theory of Sec. theory: P(n) from the measured v and D;
    local exponent r(n) against s n; coherence of a memory against storage time."""
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 58 * MM), constrained_layout=True)
    t1, t2 = load('tail.dat'), load('tail_s0p1.dat')
    cases = ((0.1, t2['P_N400'], CA, 12), (1.0, t1['P_A_s1'], CB, 5))
    for s, P, col, n_ref in cases:
        n = np.arange(P.size)
        ax[0].loglog(n[1:], P[1:], color=col, label=r'ME, $s=%g$' % s)
        _, _, nf, Pf = fp_from_measured(s)
        k = nf >= 3
        ax[0].loglog(nf[k], Pf[k] * P[n_ref] / Pf[n_ref - 1], '--', color='k', lw=0.7)
        fit = (nf >= 30) & (nf <= 200)
        print('s %g: FP slope 30-200 %.2f, ME slope %.2f' % (
            s, np.polyfit(np.log(nf[fit]), np.log(Pf[fit]), 1)[0],
            np.polyfit(np.log(n[30:161]), np.log(P[30:161]), 1)[0]))
    ax[0].plot([], [], '--', color='k', lw=0.7, label=r'from measured $v$, $D$')
    ax[0].set_xlim(1, 400)
    ax[0].set_ylim(1e-9, 1e-1)
    ax[0].set_xlabel('Fock state $n$')
    ax[0].set_ylabel('$P(n)$, case A')
    ax[0].legend(loc='lower left', fontsize=5.5)
    letter(ax[0], 'a')
    for s, col in ((0.01, CL), (0.1, CA), (1.0, CB)):
        n0, r, _, _ = fp_from_measured(s)
        ax[1].loglog(s * n0, r, 'o-', color=col, ms=2, label=r'$s=%g$' % s)
    ax[1].axhline((1 - H_A) / H_A, color='k', lw=0.7, ls='--')
    ax[1].axvline(SNC, color=CG, lw=0.5, ls=':')
    ax[1].text(0.02, 0.97, r'dashed: $(1-h)/h$; dotted: $sn=2(\nu_r/\Gamma_n)^2$',
               transform=ax[1].transAxes, va='top', fontsize=5.5)
    ax[1].set_ylim(0.5, 60)
    ax[1].set_xlabel(r'$s\,n_0$')
    ax[1].set_ylabel(r'$r(n)=-v\,n/D$')
    ax[1].legend(loc='center left', fontsize=5.5)
    letter(ax[1], 'b')
    eps = 800.0 / (2 * 1552.0)                 # Hz per quantum per axis, case A
    T = np.logspace(-3, 0.5, 300)              # up to T < 1 / eps: no revival shown
    for lab, P, col in ((r'$s=0.1$', t2['P_N400'], CA), (r'$s=1$', t1['P_A_s1'], CB),
                        (r'$s=5$', t1['P_A_s5'], CO)):
        P = np.nan_to_num(P)
        P = P / P.sum()
        n = np.arange(P.size)
        C = np.array([abs(np.sum(P * np.exp(2j * np.pi * eps * x * n)))**2 for x in T])
        ax[2].semilogx(T, C, color=col, label=lab)
        ax[2].axhline(P[0]**2, color=col, lw=0.5, ls='--')
    ax[2].set_ylim(0.6, 1.02)
    ax[2].set_xlabel('Storage time $T$ (s)')
    ax[2].set_ylabel(r'Coherence $|\phi(T)|^2$, two axes')
    ax[2].legend(loc='lower left', fontsize=5.5)
    letter(ax[2], 'c')
    save(fig, 'fom.pdf')


def closed_form(n, s, h=None, c=1.0, csat=2.25):
    """Minimal model, Fokker-Planck: v = -(1-h) R_b F, D = h (n + 1/2) R_b, R_b = csat s n / (1 + csat s n),
    F = 1 + c n / n_c (fast channel: extra friction, no extra diffusion). Unnormalised P(n), n >= 1."""
    h = H_A if h is None else h
    nc = SNC / s
    Rb = csat * s * n / (1 + csat * s * n)
    v = -(1 - h) * Rb * (1 + c * n / nc)
    D = h * (n + 0.5) * Rb
    g = np.concatenate([[0], np.cumsum(0.5 * (v[1:] / D[1:] + v[:-1] / D[:-1]) * np.diff(n))])
    return np.exp(g) / D


def fig_minimal():
    """Minimal model against the master equation: steady P(n), case A, s = 0.1 and 1. Closed form normalised
    at n = 5; chain (exact recoil kernel, ME removal rate) absolute."""
    tail = load('tail.dat')
    ch = np.genfromtxt(ROOT / 'simulation' / 'results' / 'ramp_tail_chain.dat', names=True)
    fig, axs = plt.subplots(1, 2, figsize=(WDOUBLE, 62 * MM))
    for ax, s, col, lab in ((axs[0], 0.1, 'P_A_s0p1', 'a'), (axs[1], 1.0, 'P_A_s1', 'b')):
        n = tail['n']
        P = tail[col]
        m = n >= 1
        ax.loglog(n[m], P[m], 'o', color=CA, ms=2.2, label='master equation')
        c = ch[ch['s'] == s]
        ax.loglog(c['n'][c['n'] >= 1], c['P'][c['n'] >= 1], '-', color=CO, lw=1.0, label='chain (exact kernel)')
        nn = n[m].astype(float)
        cf = closed_form(nn, s)
        i5 = np.argmin(np.abs(nn - 5))
        ax.loglog(nn, cf / cf[i5] * P[m][i5], '--', color=CG, lw=1.0, label='closed form (Fokker-Planck)')
        ax.loglog(nn, P[m][i5] * (nn / 5)**(-1 / H_A), ':', color=CL, lw=0.8, label=r'$n^{-1/h}$')
        ax.axvline(1 / (2.25 * s), color=CL, lw=0.5)
        ax.axvline(SNC / s, color=CL, lw=0.5, ls='--')
        ax.set_ylim(1e-10, 0.1)
        ax.set_xlim(0.9, 300)
        ax.set_xlabel('Fock state $n$')
        ax.set_ylabel('$P(n)$')
        ax.set_title('$s=%g$' % s, fontsize=7, loc='right')
        letter(ax, lab)
    axs[1].legend(loc='lower left', fontsize=6)
    fig.tight_layout()
    save(fig, 'minimal.pdf')


if __name__ == '__main__':
    fig_geometry()
    fig_recoil()
    fig_regimes()
    fig_tail()
    fig_phase()
    fig_fom()
    fig_minimal()
    print('scales n_lin %.0f n_anh %.0f n_depth %.0f; h_A %.3f, s n_res %.1f' % (fig_scales() + (H_A, SNC)))
