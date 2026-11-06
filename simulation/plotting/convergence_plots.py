"""Plots of the approximation and cutoff checks (studies/convergence.py).

Reads results/convergence/*.dat and writes results/convergence/convergence.{pdf,png}
(convergence checks) and discarded.{pdf,png} (the elements removed by the band
cutoff |n - m| <= B and by the 1e-11 sparsity threshold of operators.function_of_x).
Usage: python convergence_plots.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

OUT = Path(__file__).resolve().parents[1] / 'results' / 'convergence'
style.apply()
S_COL = {0.1: CL, 1.0: CA, 5.0: CO}


def load(name):
    return np.genfromtxt(OUT / name, names=True)


def main():
    fig, ax = plt.subplots(2, 3, figsize=(WDOUBLE, 110 * MM), constrained_layout=True)
    ax = ax.ravel()

    # a, b: convergence in the Fock truncation N (case A), relative to N = 160 banded
    d = load('nconv.dat')
    for i, q in enumerate(('n', 'W')):
        for s, col in S_COL.items():
            r = d[(d['case'] == 0) & (d['s'] == s)]
            ref = r[(r['method'] == 1) & (r['N'] == 160)][q][0]
            for m, mk, ls in ((0, 's', ''), (1, 'o', '-')):
                x = r[r['method'] == m]
                ax[i].plot(x['N'], x[q] / ref, ls + mk, color=col, ms=3,
                           mfc='white' if m == 0 else col,
                           label=('$s=%g$' % s) if m == 1 else None)
        rb = d[(d['case'] == 1) & (d['method'] == 1)]
        ax[i].plot(rb['N'], rb[q] / rb[q][rb['N'] == 160][0], '-^', color=CB, ms=3,
                   label='case B, $s=1$')
        ax[i].axvline(20, color=CG, lw=0.6, ls=':')
        ax[i].set_xlabel('Fock truncation $N$')
        ax[i].set_ylabel(r'$n_{\rm ss}/n_{\rm ss}(N{=}160)$' if q == 'n' else r'$W/W(N{=}160)$')
        ax[i].set_ylim(0.78, 1.06)
        letter(ax[i], 'ab'[i])
    ax[0].legend(loc='lower right', ncol=1)
    ax[1].text(0.97, 0.05, 'filled: banded, eliminated\nopen: full master equation',
               transform=ax[1].transAxes, ha='right', va='bottom')

    # c: Fock-state populations at N = 160 (case A tail, case B)
    t = load('tail.dat')
    for nm, col, lab in (('P_A_s0p1', CL, 'A, $s=0.1$'), ('P_A_s1', CA, 'A, $s=1$'),
                         ('P_A_s5', CO, 'A, $s=5$'), ('P_B_s1', CB, 'B, $s=1$')):
        ax[2].semilogy(t['n'], np.clip(t[nm], 1e-16, None), color=col, label=lab)
    ax[2].axvline(20, color=CG, lw=0.6, ls=':')
    ax[2].text(21, 1e-1, '$N=20$', color=CG)
    ax[2].set_xlim(0, 100)
    ax[2].set_ylim(1e-12, 2)
    ax[2].set_xlabel('Fock state $n$')
    ax[2].set_ylabel('$P(n)$')
    ax[2].legend(loc='upper right')
    letter(ax[2], 'c')

    # d: decay of the steady-state coherences with |n - m| (exact, N = 40)
    c = load('coherences.dat')
    for s, col in S_COL.items():
        r = c[c['s'] == s]
        ax[3].semilogy(r['k'], np.clip(r['max_abs_rho'], 1e-16, None), 'o-', color=col, ms=2.5,
                       label='$s=%g$' % s)
    ax[3].axvline(6, color=CG, lw=0.6, ls=':')
    ax[3].text(6.5, 1e-1, 'band $B=6$', color=CG)
    ax[3].set_xlim(0, 20)
    ax[3].set_ylim(1e-14, 2)
    ax[3].set_xlabel(r'$|n-m|$')
    ax[3].set_ylabel(r'$\max\,|\rho_{n,m}|$')
    ax[3].legend(loc='upper right')
    letter(ax[3], 'd')

    # e: error of the banded solution vs band width B (N = 40)
    b = load('band_check.dat')
    for s, col in S_COL.items():
        r = b[b['s'] == s]
        ax[4].semilogy(r['B'], np.abs(r['n_band'] / r['n_exact'] - 1), 'o-', color=col, ms=3,
                       label=r'$n_{\rm ss}$, $s=%g$' % s)
        ax[4].semilogy(r['B'], np.abs(r['W_band'] / r['W_exact'] - 1), 's--', color=col, ms=3,
                       mfc='white')
    ax[4].axhline(1e-3, color=CG, lw=0.6, ls=':')
    ax[4].set_xlabel('band width $B$')
    ax[4].set_ylabel('relative error')
    ax[4].set_ylim(1e-12, 1e3)                  # top kept empty for the legend
    ax[4].legend(loc='upper right')
    ax[4].text(0.03, 0.05, r'solid: $n_{\rm ss}$, dashed: $W$', transform=ax[4].transAxes)
    letter(ax[4], 'e')

    # f: error of the adiabatic elimination (N = 20)
    e = load('elim_check.dat')
    base = (e['eps'] == 0) & (e['d2'] == 0)
    for cid, DG, col, mk, lab in ((0, 3.0, CA, 'o', r'A, $\Delta=3\Gamma$'),
                                  (0, 1.0, CA, 's', r'A, $\Delta=\Gamma$'),
                                  (1, 3.0, CB, 'o', r'B, $\Delta=3\Gamma$'),
                                  (1, 1.0, CB, 's', r'B, $\Delta=\Gamma$')):
        r = e[base & (e['case'] == cid) & (e['DeltaG'] == DG)]
        ax[5].loglog(r['s'], np.abs(r['n_elim'] / r['n_full'] - 1), mk + '-', color=col, ms=3,
                     mfc='white' if DG == 1.0 else col, label=lab)
    ax[5].axhline(1e-2, color=CG, lw=0.6, ls=':')
    ax[5].set_xlabel('saturation $s$')
    ax[5].set_ylabel(r'$|n_{\rm ss}^{\rm elim}/n_{\rm ss}^{\rm full}-1|$')
    ax[5].set_ylim(1e-5, 1)
    ax[5].legend(loc='upper left', ncol=2)
    letter(ax[5], 'f')

    for ext in ('pdf', 'png'):
        fig.savefig(OUT / f'convergence.{ext}')
    plt.close(fig)


BAND, TOL = 6, 1e-11


def discarded():
    fig, axs = plt.subplots(2, 2, figsize=(WDOUBLE, 150 * MM), constrained_layout=True)
    ax = axs.ravel()
    norm = LogNorm(1e-16, 1)
    # a: steady state, exact eliminated (case A, s = 1, N = 40), band edges
    rho = np.loadtxt(OUT / 'rho_A_s1.dat', skiprows=1)
    N = rho.shape[0] // 2
    im = ax[0].imshow(np.clip(rho / rho.max(), 1e-16, None), norm=norm, cmap='viridis')
    x = np.arange(N)
    for a in range(2):
        for b in range(2):
            for sgn in (1, -1):
                y = x + sgn * BAND
                ok = (y >= 0) & (y < N)
                ax[0].plot(b * N + y[ok], a * N + x[ok], '--', color='w', lw=0.6)
    ax[0].axhline(N - 0.5, color='w', lw=0.5)
    ax[0].axvline(N - 0.5, color='w', lw=0.5)
    ax[0].set_xticks([N / 2, 3 * N / 2], [r'$g_p$', r'$g_c$'])
    ax[0].set_yticks([N / 2, 3 * N / 2], [r'$g_p$', r'$g_c$'])
    n_of = np.tile(x, 2)
    outside = np.abs(np.subtract.outer(n_of, n_of)) > BAND
    ax[0].set_title(r'$|\rho_{(a,n),(b,m)}|$ (case A, $s=1$, $N=%d$)' % N
                    + '\n' + r'discarded ($|n-m|>%d$): max %.0e, sum %.0e'
                    % (BAND, rho[outside].max(), rho[outside].sum()))
    letter(ax[0], 'a', x=-0.2, y=1.0)
    # b-d: operators f(x) without cutoff; white contour at TOL x max (zeroed below)
    for i, (name, title) in enumerate((('op_recoil.dat', r'recoil $e^{i\eta_{\rm sp}u\hat x}$, $|u|=|u|_{\max}$'),
                                       ('op_probe.dat', r'probe profile $f_p(\hat x)$'),
                                       ('op_coupling.dat', r'coupling profile $f_c(\hat x)$')), start=1):
        M = np.loadtxt(OUT / name, skiprows=1)
        mx = M.max()
        ax[i].imshow(np.clip(M / mx, 1e-16, None), norm=norm, cmap='viridis')
        ax[i].contour(M / mx, levels=[TOL], colors='w', linewidths=0.7)
        drop = M < TOL * mx
        ax[i].set_title(title + '\n' + r'zeroed: %d of %d, largest %.0e of max'
                        % (drop.sum(), M.size, (M[drop].max() / mx) if drop.any() else 0))
        ax[i].set_xlabel('$n$')
        ax[i].set_ylabel('$m$')
        letter(ax[i], 'bcd'[i - 1], x=-0.2, y=1.0)
    fig.colorbar(im, ax=ax, shrink=0.6, label='magnitude / max')
    fig.suptitle(r'white dashed: band edges $|n-m|=%d$;  white solid: threshold $10^{-11}\times$max of function_of_x' % BAND)
    for ext in ('pdf', 'png'):
        fig.savefig(OUT / f'discarded.{ext}')
    plt.close(fig)

if __name__ == '__main__':
    main()
    discarded()
