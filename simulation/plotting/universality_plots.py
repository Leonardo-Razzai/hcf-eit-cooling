"""Figures for studies/universality.py (reads the run JSON files only). Writes results/universality/:
  fixedh.pdf/png  a: fitted exponent against Delta for each s at fixed h = 0.4, 0.7 (lines: -1/h), only where
                  the fit window spans >= 0.3 decades; b, c: compensated tails P(n) n^(1/h) against n / n_c
                  (flat = power law with exponent exactly 1/h; the drop = the fast-channel cut-off)
  drift.pdf/png   local exponent gamma(n) = 2 n v / D2 - dlnD2/dln n from drift and diffusion, against n, for each
                  h (dashed: -1/h; markers: n* and n_c); band check B = 4 vs 8
Usage: python universality_plots.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

R = Path(__file__).resolve().parents[1] / 'results' / 'universality'
style.apply()
SCOL = {0.1: CL, 0.3: CA, 1.0: CB, 3.0: CO}
HCOL = {0.3: CL, 0.4: CA, 0.5: CG, 0.7: CB, 1.0: CO}


def load(sub):
    return [json.loads(f.read_text()) for f in sorted((R / sub).glob('*.json'))]


def fixedh(runs):
    if not runs:
        return
    hs = sorted({r['h'] for r in runs})
    fig, ax = plt.subplots(1, 1 + len(hs), figsize=(WDOUBLE, 58 * MM), constrained_layout=True)
    a = ax[0]
    for h in hs:
        for s in sorted({r['s'] for r in runs}):
            q = sorted([r for r in runs if r['h'] == h and r['s'] == s and r['decades'] >= 0.3 and np.isfinite(r['slope'])],
                       key=lambda r: r['DG'])
            if q:
                a.plot([r['DG'] for r in q], [-r['slope'] for r in q], 'o-' if h == hs[0] else 's--', color=SCOL.get(s, CG),
                       ms=3, lw=0.8, label='$h=%g$, $s=%g$' % (h, s))
        a.axhline(1 / h, color='k', lw=0.6, ls=':')
    a.set_xlabel(r'$\Delta/\Gamma$')
    a.set_ylabel('fitted tail exponent')
    a.legend(fontsize=4.5, ncol=2)
    for b, h in zip(ax[1:], hs):
        for r in sorted([r for r in runs if r['h'] == h], key=lambda r: (r['DG'], r['s'])):
            P = np.asarray(r['P'])
            n = np.arange(len(P))
            ok = (n >= 2) & (n < len(P) * 0.6) & (P > 0)
            y = P[ok] * n[ok]**(1 / h)
            y = y / np.median(y[(n[ok] >= r['fit_lo']) & (n[ok] <= max(r['fit_hi'], r['fit_lo'] + 1))]) if r['decades'] > 0 else y / y[0]
            b.loglog(n[ok] / r['nc'], y, color=SCOL.get(r['s'], CG), lw=0.5 + 0.25 * r['DG'] / 3, alpha=0.8)
        b.set_xlabel(r'$n/n_c$')
        b.set_ylabel(r'$P(n)\,n^{1/h}$ (normalised)')
        b.set_title('$h=%g$ (line width: $\\Delta$; colour: $s$)' % h, fontsize=6.5)
        b.set_ylim(1e-2, 1e1)
    for x, l in zip(ax, 'abc'):
        letter(x, l)
    fig.savefig(R / 'fixedh.pdf')
    fig.savefig(R / 'fixedh.png', dpi=200)
    plt.close(fig)


def drift(runs):
    if not runs:
        return
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE * 0.8, 58 * MM), constrained_layout=True)
    a = ax[0]
    for r in sorted([r for r in runs if r['B'] == 4 and r['N'] >= 1000], key=lambda r: (r['DG'], r['h'])):
        c = HCOL.get(r['h'], CG)
        ls = '-' if r['DG'] == 8 else ':'
        a.semilogx(r['n0'], r['gamma'], 'o' + ls, color=c, ms=2.5, lw=0.8, label='$h=%g$, $\\Delta=%g\\Gamma$' % (r['h'], r['DG']))
        a.axhline(-1 / r['h'], color=c, lw=0.5, ls='--')
        a.axvline(r['nc'], color=c, lw=0.4, ls='-.')
    a.set_xlabel('$n$')
    a.set_ylabel(r'local exponent $\gamma(n)$')
    a.set_ylim(-4, 1)
    a.legend(fontsize=5)
    b = ax[1]
    for r in [r for r in runs if r['N'] < 1000]:
        b.semilogx(r['n0'], r['gamma'], 'o-' if r['B'] == 4 else 's--', color=CA if r['B'] == 4 else CB, ms=2.5, lw=0.8,
                   label='$B=%d$' % r['B'])
    b.set_xlabel('$n$')
    b.set_ylabel(r'$\gamma(n)$, band check ($h=0.5$)')
    b.legend(fontsize=5.5)
    letter(a, 'a')
    letter(b, 'b')
    fig.savefig(R / 'drift.pdf')
    fig.savefig(R / 'drift.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    fixedh(load('fixedh'))
    drift(load('drift'))
    print('figures in', R)
