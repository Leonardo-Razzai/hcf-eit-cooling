"""Paper figures: simulation/results/*.dat -> paper/manuscript/figures/*.pdf.

Usage: python make_figures.py [mechanism] [saturation] [robust] [dynamics] [temperature] [lines]
(default: all). Only reads result files, so it runs in seconds.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from scipy.special import jv

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

SIM = Path(__file__).resolve().parents[1]
RESULTS = SIM / 'results'
FIG = SIM.parent / 'paper' / 'manuscript' / 'figures'
sys.path.insert(0, str(SIM))
import eitsim as E  # noqa: E402

style.apply()


def load(name):
    return np.genfromtxt(RESULTS / name, names=True)


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name)
    plt.close(fig)


def fig_mechanism():
    """Fig. 1: mode profiles and Fano absorption spectrum (computed directly)."""
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 58 * MM), constrained_layout=True)
    x = np.linspace(-1, 1, 400)
    ax[0].axvspan(-0.05, 0.05, color=CO, alpha=0.35, lw=0)
    l1, = ax[0].plot(x, jv(0, E.U01 * x), 'k')
    l2, = ax[0].plot(x, E.amp_ratio(1) * jv(1, E.U1[0] * x), color=CA)
    l3, = ax[0].plot(x, E.gfac(1) * x, '--', color=CA, lw=0.7)
    ax[0].annotate('atoms', xy=(-0.04, 1.55), xytext=(-0.62, 1.7), color='#B07800',
                   va='center', arrowprops=dict(arrowstyle='->', color='#B07800', lw=0.6))
    ax[0].set_xlabel('$x/a$')
    ax[0].set_ylabel('Field (LP$_{01}$ peak = 1)')
    ax[0].set_ylim(-2.2, 2.2)
    ax[0].legend([l1, l2, l3], ['LP$_{01}$', 'LP$_{11}$', '$G_{11}x/a$'],
                 loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=3)
    letter(ax[0], 'a')
    p = E.params(E.paper_cases()['A'])
    d = np.linspace(-2.2, 2.2, 2000)
    R = E.fano(p, d)
    ax[1].plot(d, R / R.max(), 'k')
    tr = blended_transform_factory(ax[1].transData, ax[1].transAxes)
    for dd, c, lab in ((-1, CB, 'Blue sb. $A_+$'), (0, CG, 'Carrier'), (1, CA, 'Red sb. $A_-$')):
        ax[1].axvline(dd, color=c, ls='--' if dd else ':', lw=0.7)
        ax[1].text(dd, 1.02, lab, color=c, transform=tr, ha='center', va='bottom')
    ax[1].set_xlabel(r'Two-photon detuning $\delta/\nu_r$')
    ax[1].set_ylabel(r'Absorption $R(\delta)/R_{\max}$')
    ax[1].set_ylim(0, 1.05)
    letter(ax[1], 'b')
    save(fig, 'mechanism.pdf')


def fig_saturation():
    """Fig. 3: rate and steady state versus probe strength."""
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    for key, c in (('A', CA), ('B', CB)):
        d = load(f'sat_{key}.dat')
        ax[0].loglog(d['s'], d['WME'], 'o-', color=c)
        ax[0].loglog(d['s'], d['WSM2'], '--', color=c)
        ax[0].loglog(d['s'], d['WSM1'], ':', color=c)
        ax[1].loglog(d['s'], d['nME'], 'o-', color=c)
        ax[1].loglog(d['s'], d['nSM2'], '--', color=c)
    ax[0].set_xlabel('Saturation parameter $s$')
    ax[0].set_ylabel('Cooling rate $W$ (s$^{-1}$)')
    ax[1].set_xlabel('Saturation parameter $s$')
    ax[1].set_ylabel(r'Steady-state $\langle n_x\rangle$')
    letter(ax[0], 'a')
    letter(ax[1], 'b')
    h = [Line2D([], [], color=CA, lw=1.5), Line2D([], [], color=CB, lw=1.5),
         Line2D([], [], color='k', marker='o'),
         Line2D([], [], color='k', ls='--'), Line2D([], [], color='k', ls=':')]
    t = [r'Case A ($a=15\ \mu$m)', r'Case B ($a=3.5\ \mu$m)', 'Master equation',
         r'Simple model, $\kappa=2$', r'Simple model, $\kappa=1$']
    fig.legend(h, t, loc='outside upper center', ncol=5)
    save(fig, 'saturation.pdf')


def fig_robust():
    """Fig. 4: mode impurity, two-photon detuning, trap-frequency mismatch."""
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    kw = dict(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2,
              columnspacing=0.8, handlelength=1.8)
    d = load('carrier_A.dat')
    e, x = d['kind'] == 0, d['kind'] == 1
    w0 = d['WmaxME'][e][d['epseff'][e] == 0][0]
    ax[0].plot(d['epseff'][e][1:], d['WmaxME'][e][1:] / 1e3, 'o-', color=CA, label=r'Admixture $\epsilon$')
    ax[0].plot(d['epseff'][x], d['WmaxME'][x] / 1e3, 's', color=CB, ms=4, mfc='none', mew=0.9,
               label=r'Offset $x_s$')
    ax[0].plot(d['epseff'][e][1:], d['Wunified'][e][1:] / 1e3, '--', color=CG, label='Simple model')
    ax[0].axhline(w0 / 1e3, color=CA, lw=0.6, ls=':', label=r'$\epsilon=0$')
    ax[0].set_xscale('log')
    ax[0].set_ylim(0, 3.2)
    ax[0].set_xlabel(r'$\epsilon_{\rm eff}$')
    ax[0].set_ylabel(r'Max. rate (10$^3$ s$^{-1}$)')
    ax[0].legend(**kw)
    t = load('twophoton_A.dat')
    ax[1].semilogy(t['d2'], t['n0'], 'o-', color=CA, label=r'ME, $\epsilon=0$')
    ax[1].semilogy(t['d2'], t['n5'], 's-', color=CB, label=r'ME, $\epsilon=0.05$')
    ok = np.isfinite(t['nSM']) & (t['nSM'] > 0)
    ax[1].semilogy(t['d2'][ok], t['nSM'][ok], '--', color=CG, label=r'Simple, $\epsilon=0$')
    ax[1].set_xlabel(r'Two-photon error $\delta_{\rm e}/\nu_r$')
    ax[1].set_ylabel(r'$\langle n_x\rangle$')
    ax[1].legend(**kw)
    for key, c in (('A', CA), ('B', CB)):
        m = load(f'mismatch_{key}.dat')
        ax[2].plot(m['nurel'], m['WME'] / m['WME'].max(), 'o-', color=c, ms=2.5, label=f'{key}, ME')
        ax[2].plot(m['nurel'], m['WSM'] / m['WME'].max(), '--', color=c, lw=0.8, label=f'{key}, simple')
    ax[2].set_xlabel(r'$\nu_r/\nu_r^{\rm design}$')
    ax[2].set_ylabel(r'$W/W_{\max}^{\rm ME}$')
    ax[2].legend(**kw)
    for a, s in zip(ax, 'abc'):
        letter(a, s, x=-0.3, y=1.0)
    save(fig, 'robust.pdf')


def fig_dynamics():
    """Fig. 5: relaxation (ME vs rate models) and capture from 10 uK."""
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 58 * MM), constrained_layout=True)
    d = load('dyn_A.dat')
    ax[0].plot(d['tms'], d['nME'], 'o', color=CA, ms=2.5, label='Master eq.')
    ax[0].plot(d['tms'], d['nGRM'], '-', color=CB, label='GRM, calibrated')
    ax[0].plot(d['tms'], d['nGRM1'], ':', color=CB, label=r'GRM, $c_{\rm sat}=1$')
    ax[0].plot(d['tms'], d['nSM'], '--', color=CG, label='Simple model')
    ax[0].set_yscale('log')
    ax[0].set_ylim(0.015, 8)
    ax[0].set_xlabel('Time (ms)')
    ax[0].set_ylabel(r'$\langle n_x\rangle$')
    ax[0].legend(loc='upper right')
    letter(ax[0], 'a')
    ax[0].text(0.0, 1.02, r'   Case A, $\langle n\rangle_0=2$', transform=ax[0].transAxes, va='bottom')
    labs = (('s1', CA, '-', r'$\Delta=3\Gamma$, $s=1$'), ('s03', CL, '-', r'$\Delta=3\Gamma$, $s=0.3$'),
            ('staged', CB, '-', r'$\Delta=\Gamma\to3\Gamma$'), ('LP12', 'k', ':', r'LP$_{12}$, $s=1$'))
    ylims = {'A': (0.1, 300), 'B': (0.01, 60)}
    for i, key in ((1, 'A'), (2, 'B')):
        c = load(f'capture_{key}.dat')
        for lab, col, ls, txt in labs:
            ax[i].semilogy(c['tms'], c['n_' + lab], ls, color=col, label=txt)
        ax[i].set_ylim(*ylims[key])
        ax[i].set_xlabel('Time (ms)')
        ax[i].set_ylabel(r'$\langle n_x\rangle$')
        ax[i].legend(loc='upper right')
        letter(ax[i], 'bc'[i - 1])
        ax[i].text(0.0, 1.02, r'   Case %s, from 10 $\mu$K' % key, transform=ax[i].transAxes, va='bottom')
    save(fig, 'dynamics.pdf')


def fig_temperature():
    """Fig. 6: cooling versus initial radial temperature k_B T_r / U0 (calibrated GRM)."""
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    pick = (0, 5, 10, 15, 19)
    cols = plt.cm.viridis(np.linspace(0.0, 0.85, len(pick)))[::-1]
    ylims = {'A': (0.02, 300), 'B': (0.01, 60)}
    for i, key in enumerate('AB'):
        tr = load(f'temptrace_{key}.dat')
        m = load(f'temperature_{key}.dat')
        hs = []
        for j, c in zip(pick, cols):
            h, = ax[i].semilogy(tr['tms'], tr['n%02d' % j], color=c)
            ax[i].semilogy(tr['tms'], 1 - tr['p%02d' % j], ':', color=c, lw=0.8)
            hs.append(h)
        ax[i].axhline(m['nss'][0], color=CG, lw=0.6, ls='--')
        ax[i].set_xlim(0, tr['tms'][-1])
        ax[i].set_ylim(*ylims[key])
        ax[i].set_xlabel('Time (ms)')
        ax[i].set_ylabel(r'$\langle n_x\rangle$ (solid), $1-P_0$ (dotted)')
        letter(ax[i], 'ab'[i])
        ax[i].text(0.0, 1.02, '   Case %s' % key, transform=ax[i].transAxes, va='bottom')
    fig.legend(hs, ['%.3g' % x for x in load('temperature_A.dat')['TU0'][list(pick)]],
               title=r'$k_BT_r/U_0$', loc='outside upper center', ncol=len(pick))
    for key, col in (('A', CA), ('B', CB)):
        m = load(f'temperature_{key}.dat')
        ax[2].loglog(m['TU0'], m['t90_ms'], 'o-', color=col, label=r'Case %s, $P_0\geq0.9$' % key)
        ax[2].loglog(m['TU0'], m['te_ms'], '--', color=col, label=r'Case %s, $1/e$' % key)
    ax[2].axvline(10e-6 * E.KB / E.paper_cases()['A']['U0'], color=CG, lw=0.6, ls=':')
    ax[2].set_xlim(1e-3, 1e-1)
    ax[2].set_ylim(0.02, 3e4)             # upper left kept empty for the legend
    ax[2].set_xlabel(r'Initial $k_BT_r/U_0$')
    ax[2].set_ylabel('Cooling time (ms)')
    ax[2].legend(loc='upper left', ncol=2, columnspacing=0.8, handlelength=1.6)
    top = ax[2].secondary_xaxis('top', functions=(lambda x: x * 533.0, lambda t: t / 533.0))
    top.set_xlabel(r'$T_r$ ($\mu$K)')
    letter(ax[2], 'c', y=1.1)
    save(fig, 'temperature.pdf')


# D1/D2 configurations shown in fig_lines: (line, F', g_p) -> colour, label
LINE_CFGS = (((1, 2, 1), CA, r"D1 $F'$=2, $g_p$=1 (paper)"),
             ((1, 1, 2), CL, r"D1 $F'$=1, $g_p$=2"),
             ((2, 2, 1), CB, r"D2 $F'$=2, $g_p$=1"),
             ((2, 1, 1), CO, r"D2 $F'$=1, $g_p$=1"))


def _select(d, line, Fe, gp):
    return d[(d['line'] == line) & (d['Fe'] == Fe) & (d['gp'] == gp)]


def fig_lines():
    """D1 versus D2 (stage 2): all excited hyperfine levels in the ME, two-photon
    detuning re-optimised; three-level model for reference."""
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 62 * MM), constrained_layout=True)
    dl = load('lines_delta_A.dat')
    hs = []
    for cfg, col, lab in LINE_CFGS:
        r = _select(dl, *cfg)
        ok, ok3 = r['okopt'] == 1, r['ok3'] == 1
        for i, (y3, y) in enumerate(((r['W3'], r['Wopt']), (r['n3'], r['nopt']))):
            # rates only where n_ss < 0.2 (gaps, no bridging); n_ss everywhere
            show3 = ok3 if i == 0 else np.ones_like(ok3)
            ax[i].plot(r['DeltaG'], np.where(show3, y3, np.nan), '--', color=col, lw=0.7)
            h, = ax[i].plot(r['DeltaG'], np.where(ok, y, np.nan), 'o-', color=col, ms=3)
            ax[i].plot(r['DeltaG'][~ok], y[~ok], 'o', mfc='white', color=col, ms=3)
        hs.append(h)
    ax[0].set_yscale('log')
    ax[0].set_ylim(80, 2e4)
    ax[0].set_ylabel(r'Cooling rate $W$ (s$^{-1}$)')
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-3, 20)
    ax[1].axhline(0.2, color=CG, lw=0.5, ls=':')
    ax[1].set_ylabel(r'$n_{\rm ss}$')
    for i in (0, 1):
        ax[i].set_xscale('log')
        ax[i].set_xticks([1, 2, 3, 5, 8, 12], ['1', '2', '3', '5', '8', '12'])
        ax[i].minorticks_off()
        ax[i].set_xlabel(r'$\Delta/\Gamma$')
        letter(ax[i], 'ab'[i])
        ax[i].text(0.0, 1.02, '   Case A, $s=1$', transform=ax[i].transAxes, va='bottom')
    # c: rate at Delta = 3 Gamma relative to the paper configuration (three-level)
    for j, (key, mk) in enumerate((('A', 'o'), ('B', 's'))):
        m = load(f'lines_ml_{key}.dat')
        ref = _select(m, 1, 2, 1)['W3'][0]
        for k, (cfg, col, _) in enumerate(LINE_CFGS):
            r = _select(m, *cfg)
            y = len(LINE_CFGS) - 1 - k + (0.15 if key == 'A' else -0.15)
            ax[2].plot(sorted(r['Wopt'] / ref), [y, y], '-', color=col, lw=2.2, solid_capstyle='butt')
            ax[2].plot(r['W3'][0] / ref, y, mk, mfc='white', color=col, ms=3.5)
    ax[2].axvline(1, color=CG, lw=0.5, ls=':')
    ax[2].set_xlim(0, 2.6)
    ax[2].set_ylim(-0.6, len(LINE_CFGS) - 0.4)
    ax[2].set_yticks(range(len(LINE_CFGS)), [lab for _, _, lab in LINE_CFGS][::-1])
    ax[2].set_xlabel(r'$W/W_{\rm paper}$ at $\Delta=3\Gamma$, $s=1$')
    letter(ax[2], 'c', x=-0.75)
    ax[2].text(0.0, 1.02, '   Cases A (upper) and B (lower)', transform=ax[2].transAxes, va='bottom')
    extra = [Line2D([], [], color='k', ls='--', lw=0.7),
             Line2D([], [], marker='o', ls='', mfc='white', color='k', ms=3),
             Line2D([], [], color='k', lw=2.2),
             Line2D([], [], marker='s', ls='', mfc='white', color='k', ms=3.5)]
    fig.legend(hs + extra, [lab for _, _, lab in LINE_CFGS]
               + ['three-level (a, b)', r'$n_{\rm ss}>0.2$: rate ill-defined (a, b)',
                  'multi-level, both signs (c)', 'three-level, A: circle, B: square (c)'],
               loc='outside upper center', ncol=4)
    save(fig, 'lines.pdf')


FIGURES = {'mechanism': fig_mechanism, 'saturation': fig_saturation,
           'robust': fig_robust, 'dynamics': fig_dynamics,
           'temperature': fig_temperature, 'lines': fig_lines}

if __name__ == '__main__':
    for name in sys.argv[1:] or FIGURES:
        FIGURES[name]()
