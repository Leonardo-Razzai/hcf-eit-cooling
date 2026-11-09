"""Plots of the (bp, nu_r/omega_rec) cooling map (studies/universal.py map).

Writes results/universal/map.{pdf,png} (P0, W0, tail and N-sensitivity on the plane),
results/universal/collapse.{pdf,png} (is the tail a function of h alone?) and
results/universal/target.{pdf,png} (best rate that reaches a target P0, over the
six (Delta, s) settings of the map: coarse, so steps between cells are partly the grid).
Usage: python universal_plots.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

R = Path(__file__).resolve().parents[1] / 'results'
OUT = R / 'universal'
OUT.mkdir(exist_ok=True)
style.apply()
ALPHA = 0.397
P0_COLD = 0.85        # W0 is meaningful only above this, and for n_ss < 1
D = np.genfromtxt(R / 'universal_map.dat', names=True)
NU = np.unique(D['nu_wrec'])
BP = np.unique(D['bp'])
NU_A, NU_B = NU[2], NU[5]
EDGE_NU = np.sqrt(NU[1:] * NU[:-1])
EDGE_NU = np.r_[NU[0]**2 / EDGE_NU[0], EDGE_NU, NU[-1]**2 / EDGE_NU[-1]]
EDGE_BP = (BP[1:] + BP[:-1]) / 2
EDGE_BP = np.r_[2 * BP[0] - EDGE_BP[0], EDGE_BP, 2 * BP[-1] - EDGE_BP[-1]]


def grid(q, DG, s, largest_N=True):
    Z = np.full((BP.size, NU.size), np.nan)
    for i, b in enumerate(BP):
        for j, n in enumerate(NU):
            r = D[(D['bp'] == b) & (D['nu_wrec'] == n) & (D['DeltaG'] == DG) & (D['s'] == s)]
            r = r[np.argsort(r['N'])]
            Z[i, j] = q(r[-1] if largest_N else r)
    return Z


def key(fig, lines=True, markers=True, ncol=6):
    """Shared legend: what the lines and the markers on the map panels mean."""
    hs, ls = [], []
    if lines:
        for lab, st in ((r'$h=1/3$', ':'), (r'$h=1/2$', '--'), (r'$h=1$', '-')):
            hs.append(Line2D([], [], color='k', ls=st, lw=0.7))
            ls.append(lab)
    if markers:
        for lab, mk in ((r"Case A, $b_p=1/2$ (D1 $F'=2$)", 'o'), (r"Case A, $b_p=5/6$ (D1 $F'=1$)", 's'), (r'Case B, $b_p=1/2$', '^')):
            hs.append(Line2D([], [], marker=mk, ls='', mfc='white', mec='k', ms=3.5, mew=0.6))
            ls.append(lab)
    fig.legend(hs, ls, loc='outside upper center', ncol=ncol, fontsize=6)


def decorate(ax):
    """h lines and the real configurations."""
    x = np.geomspace(EDGE_NU[0], EDGE_NU[-1], 200)
    for h, ls in ((1 / 3, ':'), (1 / 2, '--'), (1.0, '-')):
        ax.plot(x, ALPHA / (h * x), ls, color='k', lw=0.7)
    for nu, b, mk in ((NU_A, 1 / 2, 'o'), (NU_A, 5 / 6, 's'), (NU_B, 1 / 2, '^')):
        ax.plot(nu, b, mk, mfc='white', mec='k', ms=3.5, mew=0.6)
    ax.set_xscale('log')
    ax.set_xlim(EDGE_NU[0], EDGE_NU[-1])
    ax.set_ylim(EDGE_BP[0], EDGE_BP[-1])
    ax.set_xticks([1, 2, 5, 10], ['1', '2', '5', '10'])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_yticks(BP, ['1/6', '1/3', '1/2', '2/3', '5/6'])
    ax.set_xlabel(r'$\nu_r/\omega_{\rm rec}$')
    ax.set_ylabel(r'$b_p$')


def panel(fig, ax, Z, label, norm=None, cmap='Blues'):
    m = ax.pcolormesh(EDGE_NU, EDGE_BP, Z, cmap=cmap, norm=norm, shading='flat')
    fig.colorbar(m, ax=ax, label=label, fraction=0.05, pad=0.02)
    decorate(ax)


def fig_map():
    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots(2, 4, figsize=(WDOUBLE, 98 * MM), constrained_layout=True)
    for k, DG in enumerate((3.0, 6.0)):
        x = grid(lambda r: 1 - r['P0'], DG, 1.0)
        panel(fig, ax[k, 0], x, r'$1-P_0$', LogNorm(3e-3, 0.45))
        w = grid(lambda r: r['W0'], DG, 1.0)
        w[(grid(lambda r: r['P0'], DG, 1.0) < P0_COLD) | (grid(lambda r: r['n'], DG, 1.0) >= 1)] = np.nan
        panel(fig, ax[k, 1], w, r'$W_0/\nu_r$', LogNorm(0.008, 0.11))
        t = grid(lambda r: r['n'] / (1 - r['P0']), DG, 1.0)
        panel(fig, ax[k, 2], t, r'$\langle n\rangle/(1-P_0)$', LogNorm(1, 60))
        dn = grid(lambda r: 100 * (r['n'][-1] / r['n'][0] - 1) if r.size > 1 else np.nan, DG, 1.0, largest_N=False)
        panel(fig, ax[k, 3], dn, r'$n_{\rm ss}$ change, $N=120\to250$ (%)', plt.Normalize(0, 45))
        ax[k, 0].text(0.0, 1.02, r'   $\Delta=%g\Gamma$, $s=1$' % DG, transform=ax[k, 0].transAxes, va='bottom')
    for a, l in zip(ax.ravel(), 'abcdefgh'):
        letter(a, l)
    key(fig)
    fig.savefig(OUT / 'map.pdf')
    fig.savefig(OUT / 'map.png')
    plt.close(fig)


def tail_ratio_model(h, nx):
    n = np.arange(1, 200000.0)
    w = n**(-1 / h) * np.exp(-n / nx)
    return (w @ n) / w.sum()


def fig_collapse():
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 70 * MM), constrained_layout=True)
    cols = dict(zip(np.round(BP, 3), (CG, CB, CA, CL, CO)))
    mks = dict(zip(np.round(BP, 3), 'vs^Do'))
    for a, DG in zip(ax, (3.0, 6.0)):
        for b in BP:
            xs, ys = [], []
            for n in NU:
                r = D[(D['bp'] == b) & (D['nu_wrec'] == n) & (D['DeltaG'] == DG) & (D['s'] == 1.0)]
                r = r[np.argmax(r['N'])]
                xs.append(r['h'])
                ys.append(r['n'] / (1 - r['P0']))
            k = np.round(b, 3)
            a.plot(xs, ys, mks[k], color=cols[k], ms=3.5, label=r'$b_p=%.2f$' % b)
        h = np.geomspace(0.03, 2.5, 60)
        nc = 2 * DG**2
        a.plot(h, [tail_ratio_model(x, 2 * x * nc) for x in h], '-', color='k', lw=0.8,
               label=r'$P\propto n^{-1/h}e^{-n/2hn_c}$')
        for v, ls in ((1 / 3, ':'), (0.5, '--'), (1.0, '-')):
            a.axvline(v, color=CG, lw=0.5, ls=ls)
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xticks([0.05, 0.1, 0.2, 0.5, 1, 2], ['0.05', '0.1', '0.2', '0.5', '1', '2'])
        a.xaxis.set_minor_locator(plt.NullLocator())
        a.set_yticks([1, 2, 5, 10, 20, 50], ['1', '2', '5', '10', '20', '50'])
        a.set_xlabel(r'$h=\alpha\,\omega_{\rm rec}/(b_p\nu_r)$')
        a.set_ylabel(r'$\langle n\rangle/(1-P_0)$')
        a.set_ylim(0.9, 120)
        a.legend(loc='upper left', ncol=2, fontsize=5.5)
        a.text(0.0, 1.02, r'   $\Delta=%g\Gamma$, $s=1$' % DG, transform=a.transAxes, va='bottom')
    letter(ax[0], 'a')
    letter(ax[1], 'b')
    hs = [Line2D([], [], color=CG, ls=st, lw=0.5) for st in (':', '--', '-')]
    fig.legend(hs, [r'$h=1/3$', r'$h=1/2$', r'$h=1$'], loc='outside upper center', ncol=3, fontsize=6, title=None)
    fig.savefig(OUT / 'collapse.pdf')
    fig.savefig(OUT / 'collapse.png')
    plt.close(fig)


def fig_target():
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    from matplotlib.colors import LogNorm
    for a, tgt, l in zip(ax, (0.95, 0.98, 0.99), 'abc'):
        Z = np.full((BP.size, NU.size), np.nan)
        for i, b in enumerate(BP):
            for j, n in enumerate(NU):
                r = D[(D['bp'] == b) & (D['nu_wrec'] == n) & (D['P0'] >= tgt)]
                if r.size:
                    Z[i, j] = r['W0'].max()
        panel(fig, a, Z, r'best $W_0/\nu_r$', LogNorm(0.008, 0.11))
        a.text(0.0, 1.02, r'   $P_0\geq%g$' % tgt, transform=a.transAxes, va='bottom')
        letter(a, l)
    key(fig)
    fig.savefig(OUT / 'target.pdf')
    fig.savefig(OUT / 'target.png')
    plt.close(fig)


if __name__ == '__main__':
    fig_map()
    fig_collapse()
    fig_target()
