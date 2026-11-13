"""Figures for the TE01 two-axis scheme (background report). Writes results/vector/*.pdf:

  modes.pdf     the LP11 group near the axis: TE01, TM01, HE21e, HE21o and the two linear LP11
                (intensity and transverse polarisation; on-axis E_z marked)
  schemes.pdf   minimal model (two Lambda channels with separate final states), geometry with the
                field along the fibre, and the frequency components
  rb.pdf        87Rb D1 realisations: B along z (sigma+ coupling on F = 1, F' = 2) and B along x
                (pi coupling); squared couplings from eitsim.atoms; spectator sublevels
  lp11_hf.pdf   best LP11-switched scheme of the 2D Zeeman runs: 87Rb D1 F'=1 level scheme with couplings, repumpers
                and the LP11 orientation switching
  results.pdf   2D master equation: P00 and W00 / W0 for the four schemes (vector_twoaxes.dat,
                vector_shared.dat)
Usage: python te01_plots.py
"""
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from scipy.special import jv

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import eitsim as E  # noqa: E402
from eitsim.atoms import zeeman_coupling as zc  # noqa: E402

R = ROOT / 'results'
OUT = R / 'vector'
OUT.mkdir(exist_ok=True)
style.apply()
GREEN = '#009E73'
U1 = E.U1[0]


def frac(x):
    return str(Fraction(x).limit_denominator(48))


# ------------------------------------------------------------------ modes
def modes():
    n = 241
    X, Y = np.meshgrid(np.linspace(-1.08, 1.08, n), np.linspace(-1.08, 1.08, n))
    Rr, Ph = np.hypot(X, Y), np.arctan2(Y, X)
    q = np.linspace(-0.72, 0.72, 7)
    QX, QY = np.meshgrid(q, q)
    QR, QP = np.hypot(QX, QY), np.arctan2(QY, QX)
    J = jv(1, U1 * QR)
    ok = (QR < 0.9) & (QR > 0.05)
    # (title, direction angle of E_t as a function of phi, intensity weight, E_z on the axis, div E_t)
    M = [('TE$_{01}$', lambda p: p + np.pi / 2, None, False, r'$E_t\propto(-y,x)$'),
         ('TM$_{01}$', lambda p: p, None, True, r'$E_t\propto(x,y)$'),
         ('HE$_{21}^{e}$', lambda p: -p, None, False, r'$E_t\propto(x,-y)$'),
         ('HE$_{21}^{o}$', lambda p: np.pi / 2 - p, None, False, r'$E_t\propto(y,x)$'),
         ('LP$_{11}$, $\\hat x$ along lobes', 'x', 'cos', True, r'$=\frac{1}{2}$(TM$_{01}$+HE$_{21}^e$)'),
         ('LP$_{11}$, $\\hat x$ across lobes', 'x', 'sin', False, r'$=\frac{1}{2}$(HE$_{21}^o-$TE$_{01}$)')]
    fig, ax = plt.subplots(1, 6, figsize=(WDOUBLE, 42 * MM))
    for a, (title, d, wt, ez, sub) in zip(ax, M):
        f = jv(1, U1 * Rr)
        if wt is None:
            I = f**2
            th = d(QP)
            U, V = J * np.cos(th), J * np.sin(th)
        else:
            w = np.cos(Ph) if wt == 'cos' else np.sin(Ph)
            I = (f * w)**2
            wq = np.cos(QP) if wt == 'cos' else np.sin(QP)
            U, V = J * wq, 0 * QX
        I = np.where(Rr <= 1, I, np.nan)
        a.imshow(I, extent=(-1.08, 1.08, -1.08, 1.08), origin='lower', cmap='Blues')
        a.quiver(QX[ok], QY[ok], U[ok], V[ok], color='k', scale=3.2, width=0.014, headwidth=4, pivot='mid')
        a.add_patch(Circle((0, 0), 1, fill=False, lw=0.8, ec=CG))
        mk = dict(marker='o', ms=4.5, mfc=CO if ez else 'white', mec='k', mew=0.6)
        a.plot(0, 0, **mk)
        a.set_title(title, fontsize=6.5, pad=2)
        a.text(0.5, -0.1, sub, transform=a.transAxes, ha='center', va='top', fontsize=5.5)
        a.set_axis_off()
        a.set_aspect('equal')
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker='o', ls='', mfc=CO, mec='k', ms=4.5), Line2D([], [], marker='o', ls='', mfc='white', mec='k', ms=4.5)]
    fig.legend(hs, [r'$E_z\neq0$ on the axis ($\nabla\cdot E_t\neq0$)', r'node in all components'], loc='lower center',
               ncol=2, fontsize=6, bbox_to_anchor=(0.5, -0.06))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.2, wspace=0.08)
    fig.savefig(OUT / 'modes.pdf', bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------ helpers for level diagrams
def level(a, x, y, w=0.36, lab=None, c='k', lw=1.2, dy=-0.13, fs=5.5):
    a.plot([x - w / 2, x + w / 2], [y, y], color=c, lw=lw, solid_capstyle='butt')
    if lab:
        a.text(x, y + dy, lab, ha='center', va='top' if dy < 0 else 'bottom', fontsize=fs, color=c)


def arrow(a, p0, p1, c, ls='-', lw=1.1, alpha=1.0, z=3):
    a.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=6, color=c, lw=lw, ls=ls,
                                alpha=alpha, zorder=z, shrinkA=0, shrinkB=0))


# ------------------------------------------------------------------ generic scheme
def schemes():
    fig = plt.figure(figsize=(WDOUBLE, 62 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 0.95], wspace=0.3)
    # a: minimal model
    a = fig.add_subplot(gs[0, 0])
    level(a, -1.0, 2.0, lab=r'$|e_1\rangle$', dy=0.08)
    level(a, 1.0, 2.0, lab=r'$|e_2\rangle$', dy=0.08)
    level(a, 0.0, 0.0, lab=r'$|g_p\rangle$')
    level(a, -2.0, 0.35, lab=r'$|g_{c1}\rangle$')
    level(a, 2.0, 0.35, lab=r'$|g_{c2}\rangle$')
    arrow(a, (0.0, 0.0), (-1.0, 2.0), CA)
    arrow(a, (0.0, 0.0), (1.0, 2.0), CA, ls='--')
    arrow(a, (-2.0, 0.35), (-1.0, 2.0), CB)
    arrow(a, (2.0, 0.35), (1.0, 2.0), CB)
    a.text(-0.75, 0.95, r'$\Omega_p f(x)$' + '\n' + r'($E_y$)', color=CA, fontsize=5.5, ha='right')
    a.text(0.75, 0.95, r'$\Omega_p f(y)$' + '\n' + r'($E_x$)', color=CA, fontsize=5.5, ha='left')
    a.text(-1.75, 1.35, r'$\Omega_c$', color=CB, fontsize=5.5, ha='right')
    a.text(1.75, 1.35, r'$\Omega_c$', color=CB, fontsize=5.5, ha='left')
    a.text(0, -0.75, 'cools $x$ via $e_1\\to g_{c1}$,\n$y$ via $e_2\\to g_{c2}$', ha='center', fontsize=5.5)
    a.set_xlim(-2.6, 2.6)
    a.set_ylim(-1.1, 2.5)
    a.set_axis_off()
    letter(a, 'a', x=-0.02, y=0.98)
    # b: geometry
    b = fig.add_subplot(gs[0, 1])
    b.add_patch(Rectangle((0, -0.5), 6.0, 1.0, fc='0.95', ec=CG, lw=0.8))
    b.plot([0, 6.0], [0, 0], color=CG, lw=0.4, ls=':')
    for k in range(9):
        b.plot(0.7 + 0.6 * k, 0, 'o', ms=2.2, color='k')
    arrow(b, (-0.9, 0.32), (-0.05, 0.32), CA, lw=1.4)
    arrow(b, (-0.9, -0.05), (-0.05, -0.05), CB, lw=1.4)
    arrow(b, (-0.9, -0.38), (-0.05, -0.38), '#8C2D19', lw=1.4)
    b.text(-1.0, 0.32, r'probe TE$_{01}$ ($\omega_1,\omega_2$)', ha='right', va='center', fontsize=5.5, color=CA)
    b.text(-1.0, -0.05, r'coupling LP$_{01}$, $\sigma^+$', ha='right', va='center', fontsize=5.5, color=CB)
    b.text(-1.0, -0.38, r'lattice LP$_{01}$, 1064 nm', ha='right', va='center', fontsize=5.5, color='#8C2D19')
    arrow(b, (2.0, 0.85), (4.0, 0.85), 'k', lw=1.4)
    b.text(3.0, 0.95, r'$\mathbf{B}\parallel z$ (solenoid), $\mu_BB/h\approx0.1$–0.3 MHz', ha='center', va='bottom', fontsize=5.5)
    b.text(3.0, -0.62, 'atoms on the axis, at the TE$_{01}$ node;\n' + r'TE$_{01}$ = $\sigma^-$ part ($\propto x+iy$) + $\sigma^+$ part ($\propto x-iy$)', ha='center', va='top', fontsize=5.5)
    b.set_xlim(-3.9, 6.3)
    b.set_ylim(-1.6, 1.8)
    b.set_axis_off()
    letter(b, 'b', x=0.0, y=0.98)
    # c: frequencies
    c = fig.add_subplot(gs[0, 2])
    c.axhline(0, color='k', lw=0.6)
    for xv, h, col, lab in ((-1.0, 1.0, CA, r'$\omega_1$ ($\sigma^+$)'), (0.0, 1.25, CA, r'$\omega_2$ ($\sigma^-$)'),
                            (2.4, 0.8, CB, r'$\omega_c$' + '\n' + 'coupling')):
        c.plot([xv, xv], [0, h], color=col, lw=1.6)
        c.text(xv, h + 0.05, lab, ha='center', va='bottom', fontsize=5.5, color=col)
    c.annotate('', xy=(-1.0, 0.45), xytext=(0.0, 0.45), arrowprops=dict(arrowstyle='<->', lw=0.6))
    c.text(-1.15, 0.45, r'$\mu_BB/h$' + '\n' + r'$+\delta_1-\delta_2$', ha='right', va='center', fontsize=5)
    c.text(1.2, -0.12, r'$\omega_c-\omega_{1,2}=\omega_{\rm hfs}$ + Zeeman + $\delta_{1,2}$', ha='center', va='top', fontsize=5)
    c.text(1.2, -0.42, r'$\delta_k\approx$ light shift of channel $k$ ($\propto$ its $\Omega_c^2$)', ha='center', va='top', fontsize=5)
    c.text(1.2, -0.72, r'both $\approx\Delta=3\Gamma$ above $5P_{1/2},F^\prime=2$', ha='center', va='top', fontsize=5)
    c.set_xlim(-2.4, 3.2)
    c.set_ylim(-1.0, 1.6)
    c.set_axis_off()
    c.text(3.2, 1.5, 'frequency axis,\nnot to scale', ha='right', va='top', fontsize=5, color=CG)
    letter(c, 'c', x=0.0, y=0.98)
    fig.savefig(OUT / 'schemes.pdf', bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------ Rb87
def d2(F, m, Fe, me):
    if abs(m) > F or abs(me) > Fe:
        return 0.0
    return zc('D1', F, m, Fe, me, m - me)**2


def rb():
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 70 * MM))
    xs1, xs2, ye, y1, y2 = -4.2, 1.6, 2.6, 0.0, 0.55
    for k, (a, geom) in enumerate(zip(ax, ('z', 'x'))):
        for me in range(-2, 3):
            level(a, me, ye, lab=r'$%+d$' % me if me else '0', dy=0.08, fs=5)
        for m in range(-1, 2):
            level(a, xs1 + m, y1, lab=r'$%+d$' % m if m else '0', fs=5)
        for m in range(-2, 3):
            level(a, xs2 + 1.6 + m, y2, lab=r'$%+d$' % m if m else '0', fs=5)
        a.text(-2.9, ye, r"$5P_{1/2}\,F'=2$", ha='right', va='center', fontsize=6)
        a.text(xs1, y1 - 0.42, r'$5S_{1/2}\,F=1$ ($g_c$)', ha='center', fontsize=6)
        a.text(xs2 + 1.6, y2 - 0.42, r'$5S_{1/2}\,F=2$ ($g_p$)', ha='center', fontsize=6)
        gp = (xs2 + 1.6 + 1, y2)
        a.add_patch(Circle(gp, 0.13, fc=CO, ec='k', lw=0.5, zorder=5))

        def P(m):
            return (xs2 + 1.6 + m, y2)

        def C(m):
            return (xs1 + m, y1)

        def X(me):
            return (me, ye)
        if geom == 'z':
            # (m, q, me, mc, ls, probe label, probe label xy, coupling label xy)
            chans = [(1, +1, 2, 1, '-', r'$\sigma^+$, $\omega_1$', (3.5, 1.85), (-1.1, 0.5)),
                     (1, -1, 0, -1, '--', r'$\sigma^-$, $\omega_2$', (0.2, 1.3), (-2.95, 1.3))]
            spect = [(2, -1, 1, 0), (0, +1, 1, 0)]          # (m, q, me, mc): |2,2> sigma-, |2,0> sigma+
            ttl = r'a   $\mathbf{B}\parallel z$: probe $\sigma^\pm$ (TE$_{01}$), coupling $\sigma^+$ (LP$_{01}$)'
        else:
            chans = [(1, 0, 1, 1, '-', r'$\pi$ ($E_x\propto y$)', (3.5, 1.85), (-1.1, 0.5)),
                     (1, -1, 0, 0, '--', r'$\sigma^-$ ($E_y\propto x$)', (0.0, 1.25), (-2.45, 1.3))]
            spect = [(1, +1, 2, 2)]
            ttl = r'b   $\mathbf{B}\parallel x$: probe $\pi$ + $\sigma^\pm$ (TE$_{01}$), coupling $\pi$ (LP$_{01}$, $\hat x$)'
        for m, q, me, mc, ls, lab, pxy, cxy in chans:
            pd, cd = d2(2, m, 2, me), d2(1, mc, 2, me)
            arrow(a, P(m), X(me), CA, ls=ls, lw=1.2)
            arrow(a, C(mc), X(me), CB, ls=ls, lw=1.2)
            a.text(*pxy, '%s\n$d^2$=%s' % (lab, frac(pd)), fontsize=5, color=CA, ha='left', va='center')
            a.text(*cxy, '$d^2$=%s' % frac(cd), fontsize=5, color=CB, ha='right' if ls == '--' else 'left', va='center')
        for m, q, me, mc in spect:
            arrow(a, P(m), X(me), CA, ls=':', lw=0.7, alpha=0.6)
            if abs(mc) <= 1:
                arrow(a, C(mc), X(me), CB, ls=':', lw=0.7, alpha=0.6)
        a.set_title(ttl, fontsize=6.5, loc='left')
        a.set_xlim(-6.0, 5.4)
        a.set_ylim(-0.75, 3.1)
        a.set_axis_off()
    ax[0].text(-6.0, -0.72, r'Raman shifts ($g_1=-\frac{1}{2}$, $g_2=+\frac{1}{2}$): $|2,1\rangle\to|1,1\rangle$: $-\mu_BB$; $|2,1\rangle\to|1,-1\rangle$: 0.'
               '\nDotted: spectators $|2,2\\rangle$ ($\\sigma^-$ only) and $|2,0\\rangle$ ($\\sigma^+$ only), both via $|1,0\\rangle$.', fontsize=5, va='top')
    ax[1].text(-6.0, -0.72, r'Raman shifts: $|2,1\rangle\to|1,1\rangle$: $-\mu_BB$; $|2,1\rangle\to|1,0\rangle$: $-\frac{1}{2}\mu_BB$.'
               '\nDotted: $\\sigma^+$ to $|2^\\prime,2\\rangle$ has no partner in $F=1$ ($\\pi$ coupling).', fontsize=5, va='top')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.12, wspace=0.05)
    fig.savefig(OUT / 'rb.pdf', bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------ results
def results():
    d = np.genfromtxt(R / 'vector_twoaxes.dat', names=True)
    sh = np.genfromtxt(R / 'vector_shared.dat', names=True)
    labels = ['LP$_{11}$\nswitched', 'TE$_{01}$\none $\\Lambda$', 'TE$_{01}$, 2 $\\Lambda$\nshared $g_c$', 'TE$_{01}$, 2 $\\Lambda$\nseparate $g_c$']
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 55 * MM), constrained_layout=True)
    xs = np.arange(4)
    k = 0
    for nw, case, cols in ((8.43, 'B', (CB, '#F0A57A')), (1.97, 'A', (CA, CL))):
        for s, col in zip((1.0, 0.5), cols):
            r = d[(np.abs(d['nu_wrec'] - nw) < 0.05) & (d['s'] == s)]
            r = r[np.argsort(r['scheme'])]
            q = sh[(np.abs(sh['nu_wrec'] - nw) < 0.05) & (sh['s'] == s)][0]
            P = r['P00'].copy()
            P[2] = q['P00']
            W = r['W00'] / r['W0_1D']
            W[2] = q['W00'] / q['W0_1D']
            W[P < 0.85] = np.nan
            off = (k - 1.5) * 0.2
            ax[0].bar(xs + off, P, 0.19, color=col, label='case %s, $s=%g$' % (case, s))
            ax[1].bar(xs + off, np.nan_to_num(W), 0.19, color=col)
            k += 1
    for a in ax:
        a.axvspan(2.55, 3.45, color=GREEN, alpha=0.12, lw=0, zorder=0)
        a.set_xticks(xs, labels, fontsize=6)
        a.set_xlim(-0.55, 3.55)
    ax[0].set_ylim(0, 1.15)
    ax[0].set_ylabel(r'$P_{00}$ (both axes in $n=0$)')
    ax[1].set_ylim(0, 0.55)
    ax[1].set_ylabel(r'$W_{00}/W_0$ (only where $P_{00}\geq0.85$)')
    ax[0].legend(loc='upper left', ncol=2, fontsize=5.5)
    for x in (1, 2):
        ax[1].text(x, 0.01, 'not cold', ha='center', va='bottom', fontsize=5.5, color=CG)
    letter(ax[0], 'a', x=-0.12)
    letter(ax[1], 'b', x=-0.12)
    fig.savefig(OUT / 'results.pdf')
    plt.close(fig)


# ------------------------------------------------------------------ best LP11-switched scheme
def lp11hf():
    """Best LP11-switched scheme of the 2D Zeeman runs (lp11_hf): 87Rb D1 F'=1, probe pi on F=2 (LP11, polarised
    along the quantisation axis x, orientation switched), coupling linear perpendicular (sigma+ and sigma-) on
    F=1 (LP01, y), weak repumpers on both legs with the complementary polarisation; no bias field."""
    fig = plt.figure(figsize=(WDOUBLE, 72 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.3, 1.0], wspace=0.05)
    a = fig.add_subplot(gs[0, 0])
    ye, y1, y2, xs1, xs2 = 2.6, 0.0, 0.45, -3.6, 2.6
    for me in (-1, 0, 1):
        level(a, me, ye, lab=('%+d' % me) if me else '0', dy=0.08, fs=5)
    for m in (-1, 0, 1):
        level(a, xs1 + m, y1, lab=('%+d' % m) if m else '0', fs=5)
    for m in range(-2, 3):
        level(a, xs2 + m, y2, lab=('%+d' % m) if m else '0', fs=5, c=CG if abs(m) == 2 else 'k')
    a.text(-1.5, ye, r"$5P_{1/2}\,F'=1$", ha='right', va='center', fontsize=6)
    a.text(xs1, y1 - 0.42, r'$5S_{1/2}\,F=1$ ($g_c$)', ha='center', fontsize=6)
    a.text(xs2, y2 - 0.42, r'$5S_{1/2}\,F=2$ ($g_p$)', ha='center', fontsize=6)
    # probe pi: |2,m> -> |1',m>
    for m in (-1, 0, 1):
        arrow(a, (xs2 + m, y2), (m + 0.06 * m, ye), CA, lw=1.2)
    a.text(4.0, 1.75, r'probe $\pi$' + '\n' + r'(LP$_{11}$, $\hat x$)' + '\n' + r'$d^2$ = %s, %s, %s' % (frac(d2(2, -1, 1, -1)), frac(d2(2, 0, 1, 0)), frac(d2(2, 1, 1, 1))) + '\n' + r'for $m=-1,0,+1$',
           color=CA, fontsize=5.5, ha='left', va='center')
    # coupling sigma+ (solid) and sigma- (dashed): |1,m> -> |1',m+-1>
    for m, me in ((-1, 0), (0, 1)):
        arrow(a, (xs1 + m, y1), (me - 0.05, ye), CB, lw=1.0)
    for m, me in ((1, 0), (0, -1)):
        arrow(a, (xs1 + m, y1), (me + 0.05, ye), CB, ls='--', lw=1.0)
    a.text(-6.3, 1.75, r'coupling' + '\n' + r'$\sigma^+$ (solid), $\sigma^-$ (dashed)' + '\n' + r'(LP$_{01}$, $\hat y$)' + '\n' + r'$d^2$ = %s each' % frac(d2(1, 0, 1, 1)),
           color=CB, fontsize=5.5, ha='left', va='center')
    # repumpers (dotted): F=2 sigma+- from |2,+-2>, F=1 pi from |1,+-1>
    for m in (-2, 2):
        arrow(a, (xs2 + m, y2), (np.sign(m) * 1.0 + 0.12 * np.sign(m), ye), CG, ls=':', lw=0.7, alpha=0.8)
    for m in (-1, 1):
        arrow(a, (xs1 + m, y1), (m - 0.12 * m, ye), CG, ls=':', lw=0.7, alpha=0.8)
    for m in (-2, 2):
        a.text(xs2 + m, y2 + 0.08, 'dark', fontsize=4.5, color=CG, ha='center', va='bottom')
    a.set_xlim(-6.4, 6.4)
    a.set_ylim(-0.75, 3.0)
    a.set_axis_off()
    a.text(-6.4, 3.0, 'a', fontsize=8, fontweight='bold', va='top')
    # b: switching of the LP11 orientation
    n = 121
    X, Y = np.meshgrid(np.linspace(-1.05, 1.05, n), np.linspace(-1.05, 1.05, n))
    Rr, Ph = np.hypot(X, Y), np.arctan2(Y, X)
    sub = gs[0, 1].subgridspec(2, 1, hspace=0.35)
    for k, (w, ttl) in enumerate(((np.cos, r'stage 1 ($\tau$): lobes along $x$, cools $x$'), (np.sin, r'stage 2 ($\tau$): lobes along $y$, cools $y$'))):
        b = fig.add_subplot(sub[k, 0])
        b.imshow(np.where(Rr <= 1, (jv(1, U1 * Rr) * w(Ph))**2, np.nan), extent=(-1.05, 1.05, -1.05, 1.05), origin='lower', cmap='Blues')
        for xx, yy in ((-0.55, 0.0), (0.55, 0.0)) if k == 0 else ((0.0, -0.55), (0.0, 0.55)):
            b.annotate('', xy=(xx + 0.22, yy), xytext=(xx - 0.22, yy), arrowprops=dict(arrowstyle='<|-|>', lw=0.6, color='k', mutation_scale=5))
        b.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, ec=CG))
        b.plot(0, 0, 'o', ms=3, mfc=CO, mec='k', mew=0.5)
        b.set_title(ttl, fontsize=5.5, pad=2)
        b.set_axis_off()
        b.set_aspect('equal')
        if k == 0:
            b.text(-0.32, 1.22, 'b', transform=b.transAxes, fontsize=8, fontweight='bold', va='top')
    fig.savefig(OUT / 'lp11_hf.pdf', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    for f in (modes, schemes, rb, results, lp11hf):
        f()
        print('ok', f.__name__)
