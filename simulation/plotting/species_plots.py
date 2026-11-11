"""Species, ions and tweezers on the cooling map (studies/species.py points).

Writes results/universal/species.{pdf,png}:
  a  overview on (nu_r/omega_rec, b_p) with lines of constant h; ions, tweezers and this fibre
     (case A and B geometry, nu_z = 300 kHz); grey box: range of the simulated map
  b  the simulated map (1 - P0, Delta = 3 Gamma, s = 1) with the fibre points
Usage: python species_plots.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import style
from style import MM, WDOUBLE, CA, CB, CG, CO, letter
import universal_plots as U

R = Path(__file__).resolve().parents[1] / 'results'
OUT = R / 'universal'
style.apply()
ALPHA = 0.397

rows = []
for line in open(R / 'species_points.txt').read().splitlines()[1:]:
    name, cat, x, bp, h = line.split()[:5]
    rows.append((name, cat, float(x), float(bp), float(h)))

LABEL = {  # short labels; None: not labelled
    'Ca+_Roos2000': None, 'Ca+_Xing2025': None, 'Ca+_Lechner2016': 'Ca$^+$ (3 expts)',
    'Yb+_Feng2020': 'Yb$^+$', 'Be+_Lin2013': 'Be$^+$', 'Mg+_Wu2024': 'Mg$^+$',
    'Rb87_Qiu2026_radial': 'Rb, EIT tw.', 'Rb87_Qiu2026_axial': 'Rb tw. axial',
    'Li6_Blodgett2023': 'Li $\\Lambda$GM', 'K39_Angonga2022': 'K GM', 'Cs_nanofibre_Pache2026': 'Cs nanofibre',
    'Sr88_Cooper2018': 'Sr 689', 'Yb171_Jenkins2022': 'Yb RSC',
}
SPEC = {'6Li': 'Li', '23Na': 'Na', '39K': 'K', '87Rb': 'Rb', '133Cs': 'Cs',
        'Yb171_556': 'Yb 556', 'Yb171_399': 'Yb 399', 'Sr87_689': 'Sr 689', 'Sr87_461': 'Sr 461'}
STYLE = {'ion': dict(marker='D', color='k'), 'tweezer': dict(marker='v', color=CO)}


def fibre_label(name):
    base, tag = name.rsplit('_fibre', 1)
    return SPEC[base], tag


def hlines(ax, x0, x1):
    x = np.geomspace(x0, x1, 200)
    for h, ls in ((1 / 3, ':'), (1 / 2, '--'), (1.0, '-')):
        ax.plot(x, ALPHA / (h * x), ls, color='k', lw=0.7, zorder=1)


def fibre_points(ax, labels=True, offs=None):
    offs = offs or {}
    for name, cat, x, bp, h in rows:
        if cat != 'fibre':
            continue
        lab, tag = fibre_label(name)
        c = CA if tag == 'A' else CB
        ax.plot(x, bp, 'o' if tag == 'A' else 's', ms=3.5, mfc=c if tag == 'A' else 'white', mec=c, mew=0.7, zorder=4)
        if labels:
            dx, dy = offs.get(name, (3, 2))
            ax.annotate(lab, (x, bp), xytext=(dx, dy), textcoords='offset points', fontsize=5, color=c,
                        ha='center' if dx == 0 else 'left' if dx > 0 else 'right', va='center')


def main():
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 72 * MM), constrained_layout=True,
                           gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    hlines(a, 0.05, 1500)
    a.add_patch(Rectangle((U.EDGE_NU[0], U.EDGE_BP[0]), U.EDGE_NU[-1] - U.EDGE_NU[0], U.EDGE_BP[-1] - U.EDGE_BP[0],
                          fc='0.9', ec='none', zorder=0))
    offs = {'Ca+_Roos2000': (3, 4), 'Ca+_Lechner2016': (0, -6), 'Ca+_Xing2025': (3, 4), 'Be+_Lin2013': (-3, 4),
            'Mg+_Wu2024': (3, -4), 'Yb+_Feng2020': (-3, -5), 'Rb87_Qiu2026_radial': (3, 0), 'Rb87_Qiu2026_axial': (0, -6),
            'Li6_Blodgett2023': (0, -6), 'K39_Angonga2022': (-3, 4), 'Cs_nanofibre_Pache2026': (3, -4),
            'Sr88_Cooper2018': (-3, 0), 'Yb171_Jenkins2022': (3, 4)}
    for name, cat, x, bp, h in rows:
        if cat == 'fibre':
            continue
        st = STYLE[cat]
        a.plot(x, bp, st['marker'], ms=3.5, color=st['color'], mec='white', mew=0.4, zorder=4)
        dx, dy = offs.get(name, (3, 0))
        if LABEL[name]:
            a.annotate(LABEL[name], (x, bp), xytext=(dx, dy), textcoords='offset points', fontsize=5,
                       ha='center' if dx == 0 else 'left' if dx > 0 else 'right', va='center', color='k')
    fibre_points(a, labels=False)
    a.set_xscale('log')
    a.set_xlim(0.05, 1500)
    a.set_ylim(0, 1.05)
    a.set_xticks([0.1, 1, 10, 100, 1000], ['0.1', '1', '10', '100', '1000'])
    a.set_xlabel(r'$\nu_r/\omega_{\rm rec}$')
    a.set_ylabel(r'$b_p$ (EIT) or $1/\kappa$ (other schemes)')
    a.text(1.0, 0.89, 'simulated map', fontsize=5, color=CG, va='bottom')

    b = ax[1]
    Z = U.grid(lambda r: 1 - r['P0'], 3.0, 1.0)
    m = b.pcolormesh(U.EDGE_NU, U.EDGE_BP, Z, cmap='Blues', norm=LogNorm(3e-3, 0.45), shading='flat', zorder=0)
    fig.colorbar(m, ax=b, label=r'$1-P_0$ ($\Delta=3\Gamma$, $s=1$)', fraction=0.05, pad=0.02)
    hlines(b, 0.05, 30)
    offs_b = {'87Rb_fibreA': (4, 3), '133Cs_fibreA': (4, -3), '39K_fibreA': (0, -7), '23Na_fibreA': (4, 3), '6Li_fibreA': (4, 3),
              'Yb171_556_fibreA': (4, -5), 'Yb171_399_fibreA': (-4, -3), 'Sr87_689_fibreA': (0, -7), 'Sr87_461_fibreA': (-4, 3),
              '87Rb_fibreB': (4, 3), '133Cs_fibreB': (4, -3), '39K_fibreB': (4, 3), '23Na_fibreB': (4, 3), '6Li_fibreB': (4, 3),
              'Yb171_556_fibreB': (4, 3), 'Yb171_399_fibreB': (0, 7), 'Sr87_689_fibreB': (4, -3), 'Sr87_461_fibreB': (0, 7)}
    fibre_points(b, offs=offs_b)
    b.set_xscale('log')
    b.set_xlim(0.08, 25)
    b.set_ylim(0.0, 1.0)
    b.set_xticks([0.1, 0.3, 1, 3, 10], ['0.1', '0.3', '1', '3', '10'])
    b.xaxis.set_minor_locator(plt.NullLocator())
    b.set_xlabel(r'$\nu_r/\omega_{\rm rec}$')
    b.set_ylabel(r'$b_p$')
    letter(a, 'a', x=-0.13)
    letter(b, 'b', x=-0.2)
    hs = [Line2D([], [], color='k', ls=st, lw=0.7) for st in (':', '--', '-')]
    ls = [r'$h=1/3$', r'$h=1/2$', r'$h=1$']
    hs += [Line2D([], [], marker='D', ls='', color='k', ms=3.5), Line2D([], [], marker='v', ls='', color=CO, ms=3.5),
           Line2D([], [], marker='o', ls='', color=CA, ms=3.5),
           Line2D([], [], marker='s', ls='', mfc='white', mec=CB, ms=3.5, mew=0.7)]
    ls += ['ions (EIT)', 'tweezers, nanofibre', 'fibre case A ($a=15\\,\\mu$m)', 'fibre case B ($a=3.5\\,\\mu$m)']
    fig.legend(hs, ls, loc='outside upper center', ncol=7, fontsize=5.5)
    fig.savefig(OUT / 'species.pdf')
    fig.savefig(OUT / 'species.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
