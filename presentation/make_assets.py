"""Images for the slide deck: schematics, equations and result plots.

Reads simulation/results, the eitsim package and notes/quantum-memory; writes presentation/assets/*.png
(300 dpi). Run from anywhere:  python presentation/make_assets.py
"""
import sys
from pathlib import Path

import numpy as np
import scipy.constants as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, Polygon
from matplotlib.lines import Line2D
from scipy.special import jv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / 'simulation' / 'results'
OUT = HERE / 'assets'
sys.path.insert(0, str(ROOT / 'simulation'))
import eitsim as E  # noqa: E402

# On WSL, register the Windows fonts so that Calibri is used (otherwise DejaVu Sans)
for _f in Path('/mnt/c/Windows/Fonts').glob('*'):
    if _f.name.lower().startswith(('calibri', 'arial', 'cambria')) and _f.suffix.lower() in ('.ttf', '.ttc'):
        try:
            font_manager.fontManager.addfont(str(_f))
        except Exception:
            pass

# deck palette (matches build_deck.js THEME)
NAVY, INK, GREY, LIGHT = '#0F1B2D', '#1A2233', '#5B6475', '#EEF2F7'
PROBE, COUPL, ATOM, GREEN, LATT, SKY = '#0072B2', '#D55E00', '#E69F00', '#009E73', '#8C2D19', '#56B4E9'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Calibri', 'Arial', 'DejaVu Sans'],
    'mathtext.fontset': 'stix', 'font.size': 13, 'axes.labelsize': 14,
    'legend.fontsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'axes.linewidth': 0.9, 'lines.linewidth': 2.0, 'lines.markersize': 5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'text.color': INK, 'legend.frameon': False, 'savefig.dpi': 300,
    'savefig.transparent': False, 'figure.facecolor': 'white',
})


def save(fig, name, transparent=False):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / name, bbox_inches='tight', pad_inches=0.06, transparent=transparent)
    plt.close(fig)


def load(name):
    return np.genfromtxt(RESULTS / name, names=True)


# ------------------------------------------------------------------ equations
EQUATIONS = {
    'eq_eta': r'$\eta_\perp=G_{11}\,\dfrac{x_0}{a},\qquad G_{11}=3.49$',
    'eq_h': r'$h=\dfrac{\alpha}{b_p}\,\dfrac{\omega_{\rm rec}}{\nu_r}\ \propto\ \dfrac{a}{b_p\,\nu_z}$',
}


def equations():
    for name, tex in EQUATIONS.items():
        fig = plt.figure(figsize=(0.1, 0.1))
        fig.text(0, 0, tex, fontsize=26, color=INK)
        save(fig, name + '.png', transparent=True)


# ------------------------------------------------------------- mode maps
def mode_maps():
    """Cross-sections: lattice/coupling in LP01, probe in LP11, atom at centre."""
    n = 401
    X, Y = np.meshgrid(np.linspace(-1.15, 1.15, n), np.linspace(-1.15, 1.15, n))
    R = np.hypot(X, Y)
    inside = R <= 1
    I01 = np.where(inside, jv(0, E.U01 * R)**2, np.nan)
    I11 = np.where(inside, (jv(1, E.U1[0] * R) * X / np.maximum(R, 1e-9))**2, np.nan)
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.9))
    for a, I, cmap, title in ((ax[0], I01, 'Oranges', 'Coupling + lattice: LP$_{01}$'),
                              (ax[1], I11, 'Blues', 'Probe: LP$_{11}$')):
        a.imshow(I, extent=(-1.15, 1.15, -1.15, 1.15), cmap=cmap, origin='lower')
        a.add_patch(Circle((0, 0), 1, fill=False, lw=2, ec=GREY))
        a.add_patch(Circle((0, 0), 0.07, fc=ATOM, ec=INK, lw=1.2, zorder=5))
        a.set_title(title, fontsize=15, color=INK, pad=8)
        a.set_axis_off()
        a.set_aspect('equal')
    for a, lab in ((ax[0], 'atoms'), (ax[1], 'node on axis')):
        a.annotate(lab, xy=(0.06, -0.08), xytext=(0.55, -1.42), fontsize=13, color=INK,
                   ha='center', arrowprops=dict(arrowstyle='->', color=INK, lw=1.2))
        a.set_xlim(-1.2, 1.2)
        a.set_ylim(-1.55, 1.2)
    save(fig, 'mode_maps.png')


def lattice_sketch():
    """Atoms in LP01 lattice pancakes; axial vs radial trap frequencies."""
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    z = np.linspace(0, 6 * np.pi, 800)
    r = np.linspace(-1, 1, 200)
    Z, Rr = np.meshgrid(z, r)
    I = np.cos(Z)**2 * np.exp(-2 * Rr**2 / 0.45)
    ax.imshow(I, extent=(0, z[-1], -1, 1), aspect='auto', cmap='Reds', origin='lower', alpha=0.85)
    ax.plot([0, z[-1]], [1, 1], color=GREY, lw=3)
    ax.plot([0, z[-1]], [-1, -1], color=GREY, lw=3)
    rng = np.random.default_rng(3)
    for k in range(1, 6):         # red-detuned lattice: atoms at the antinodes
        zc = k * np.pi
        ax.scatter(zc + rng.normal(0, 0.06, 9), np.clip(rng.normal(0, 0.2, 9), -0.45, 0.45),
                   s=26, c=ATOM, ec=INK, lw=0.6, zorder=4)
    ax.annotate('', xy=(z[-1] + 1.4, 0), xytext=(z[-1] + 0.2, 0),
                arrowprops=dict(arrowstyle='-|>', color=LATT, lw=2.4))
    ax.text(z[-1] + 0.8, 0.25, '$k$', ha='center', fontsize=16, color=LATT)
    ax.annotate('', xy=(1.5 * np.pi, 0.95), xytext=(1.5 * np.pi, -0.95),
                arrowprops=dict(arrowstyle='<->', color=PROBE, lw=2))
    ax.text(1.5 * np.pi + 0.25, 0.75, r'radial: $\nu_r\approx 7$ kHz', fontsize=14, color=PROBE)
    ax.annotate('', xy=(5 * np.pi, -1.28), xytext=(4 * np.pi, -1.28),
                arrowprops=dict(arrowstyle='<->', color=LATT, lw=2))
    ax.text(4.5 * np.pi, -1.45, r'axial: $\nu_z=300$ kHz', ha='center', va='top',
            fontsize=14, color=LATT)
    ax.set_xlim(0, z[-1] + 1.6)
    ax.set_ylim(-1.95, 1.15)
    ax.set_axis_off()
    save(fig, 'lattice.png')


# ----------------------------------------------------------- level scheme
def levels():
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.2), gridspec_kw=dict(width_ratios=[1, 1.1]))
    a = ax[0]
    lv = dict(lw=3, color=INK, solid_capstyle='butt')
    a.plot([0, 1.6], [0, 0], **lv)
    a.plot([3.4, 5.0], [0.6, 0.6], **lv)
    a.plot([1.7, 3.3], [4.0, 4.0], **lv)
    a.plot([1.7, 3.3], [4.6, 4.6], color=GREY, lw=1.5, ls='--')
    a.text(0.8, -0.45, r'$|g_p\rangle$', ha='center', fontsize=17)
    a.text(4.2, 0.15, r'$|g_c\rangle$', ha='center', fontsize=17)
    a.text(1.55, 4.0, r'$|e\rangle$', ha='right', va='center', fontsize=17)
    a.annotate('', xy=(3.55, 4.6), xytext=(3.55, 4.0), arrowprops=dict(arrowstyle='<->', color=INK, lw=1.3))
    a.text(3.7, 4.3, r'$\Delta$', va='center', fontsize=17)
    a.annotate('', xy=(2.15, 4.6), xytext=(0.55, 0.0), arrowprops=dict(arrowstyle='-|>', color=PROBE, lw=3))
    a.annotate('', xy=(2.85, 4.6), xytext=(4.45, 0.6), arrowprops=dict(arrowstyle='-|>', color=COUPL, lw=3))
    a.text(0.75, 2.6, r'$\Omega_p$' + '\nLP$_{11}$', ha='right', fontsize=16, color=PROBE)
    a.text(4.3, 2.8, r'$\Omega_c$' + '\nLP$_{01}$', ha='left', fontsize=16, color=COUPL)
    for x0, x1, y1 in ((2.35, 1.35, 0.15), (2.65, 3.65, 0.75)):
        t = np.linspace(0, 1, 200)
        xs = x0 + (x1 - x0) * t + 0.06 * np.sin(2 * np.pi * 14 * t)
        a.plot(xs, 3.85 + (y1 - 3.85) * t, color=GREY, lw=1.2)
    a.text(2.5, 1.8, r'$\Gamma$', ha='center', fontsize=16, color=GREY)
    a.set_title('a   bare $\\Lambda$ system ($^{87}$Rb D1)', loc='left', fontsize=15, fontweight='bold')
    a.set_xlim(-0.6, 5.8)
    a.set_ylim(-0.9, 5.0)
    a.set_axis_off()
    b = ax[1]
    for n in range(3):
        b.plot([0, 1.4], [1.3 * n] * 2, **lv)
        b.text(-0.15, 1.3 * n, rf'$|g_p,{n}\rangle$', ha='right', va='center', fontsize=16)
    for m in range(2):
        y = 1.3 * (m + 1)
        b.plot([3.0, 4.4], [y, y], lw=3, color=ATOM)
        b.text(4.55, y, rf'$|+,{m}\rangle$', va='center', fontsize=16, color='#9A6700')
    for n, lab in ((1, r'$\Omega_\perp$'), (2, r'$\Omega_\perp\sqrt{2}$')):
        b.annotate('', xy=(2.9, 1.3 * n), xytext=(1.5, 1.3 * n),
                   arrowprops=dict(arrowstyle='-|>', color=PROBE, lw=2.6))
        b.text(2.2, 1.3 * n + 0.15, lab, ha='center', fontsize=15, color=PROBE)
        t = np.linspace(0, 1, 200)
        xs = 3.2 + (1.3 - 3.2) * t
        ys = 1.3 * n - 0.1 + (1.3 * (n - 1) + 0.1 - (1.3 * n - 0.1)) * t + 0.05 * np.sin(2 * np.pi * 12 * t)
        b.plot(xs, ys, color=GREY, lw=1.2)
    b.text(2.75, 0.35, r'$\Gamma_n$', fontsize=15, color=GREY)
    b.set_title('b   dressed picture, $\\delta_+=\\nu_r$', loc='left', fontsize=15, fontweight='bold')
    b.set_xlim(-1.5, 5.6)
    b.set_ylim(-0.9, 3.4)
    b.set_axis_off()
    save(fig, 'levels.png')


def fano_plot():
    p = E.params(E.paper_cases()['A'])
    d = np.linspace(-2.2, 2.2, 2000)
    R = E.fano(p, d)
    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    ax.plot(d, R / R.max(), color=INK)
    for dd, c, lab in ((-1, COUPL, 'blue sb.'), (0, GREY, 'carrier'), (1, PROBE, 'red sb.')):
        ax.axvline(dd, color=c, ls='--' if dd else ':', lw=1.6)
        ax.text(dd, 1.07, lab, color=c, ha='center', va='bottom', fontsize=14)
    ax.set_xlabel(r'two-photon detuning $\delta/\nu_r$')
    ax.set_ylabel('absorption (norm.)')
    ax.set_ylim(0, 1.05)
    save(fig, 'fano.png')


# ------------------------------------------------------- setup schematic
def box(ax, x, y, w, h, text, fc, ec, tc=INK, fs=13):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.02,rounding_size=0.12',
                                fc=fc, ec=ec, lw=1.6, zorder=3))
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, color=tc, zorder=4)


def beam(ax, pts, color, lw=3.2, alpha=0.9, z=1):
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw, alpha=alpha, solid_capstyle='round', zorder=z)


def mirror(ax, x, y, angle=45, color=GREY, length=0.55):
    t = np.deg2rad(angle)
    dx, dy = 0.5 * length * np.cos(t), 0.5 * length * np.sin(t)
    ax.plot([x - dx, x + dx], [y - dy, y + dy], color=color, lw=4.5, solid_capstyle='butt', zorder=4)


def setup():
    fig, ax = plt.subplots(figsize=(13.0, 7.0))
    ax.set_xlim(0, 16)
    ax.set_ylim(-0.2, 9.0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    fx = 3.2                      # fibre axis
    # fibre (vertical), hollow core, atom strip
    ax.add_patch(Rectangle((fx - 0.45, 1.9), 0.9, 4.9, fc='#D9DEE6', ec=GREY, lw=1.5, zorder=2))
    ax.add_patch(Rectangle((fx - 0.13, 1.9), 0.26, 4.9, fc='white', ec='none', zorder=2))
    for y in np.arange(3.0, 6.3, 0.17):
        ax.add_patch(Circle((fx, y), 0.055, fc=ATOM, ec='none', zorder=3))
    ax.annotate('', xy=(fx + 0.75, 6.25), xytext=(fx + 0.75, 3.0),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=1.2))
    ax.text(fx + 0.9, 4.6, 'atom strip\nin 1D lattice\n(cm long)', fontsize=13, va='center')
    ax.text(fx - 0.6, 2.25, 'hollow-core\nfibre', fontsize=13, ha='right', va='center')
    # MOT above the tip
    ax.add_patch(Ellipse((fx, 7.55), 1.1, 0.75, fc=ATOM, alpha=0.55, ec='none', zorder=3))
    for ang in (20, 160, 200, 340):
        t = np.deg2rad(ang)
        ax.annotate('', xy=(fx + 0.5 * np.cos(t), 7.55 + 0.33 * np.sin(t)),
                    xytext=(fx + 1.25 * np.cos(t), 7.55 + 0.8 * np.sin(t)),
                    arrowprops=dict(arrowstyle='-|>', color=GREY, lw=1.4))
    ax.text(fx - 1.3, 7.55, 'MOT\n$^{87}$Rb', ha='right', va='center', fontsize=13)
    ax.annotate('', xy=(fx, 6.85), xytext=(fx, 7.2), arrowprops=dict(arrowstyle='-|>', color='#9A6700', lw=1.6))
    # top lattice beam (counter-propagating) via dichroic above MOT
    mirror(ax, fx, 8.55, 45, color=LATT)
    beam(ax, [(7.2, 8.55), (fx, 8.55), (fx, 6.85)], LATT, lw=2.6, alpha=0.65)
    box(ax, 8.9, 8.55, 3.4, 0.62, '1064 nm lattice beam 2', '#F6E3DF', LATT)
    # bottom path: combined 795 + 1064 into the fibre
    yb = 0.75
    mirror(ax, fx, yb, 45)
    beam(ax, [(fx, yb), (fx, 1.9)], LATT, lw=6.0, alpha=0.35)
    beam(ax, [(fx, yb), (fx, 1.9)], COUPL, lw=3.0)
    beam(ax, [(fx + 0.09, yb), (fx + 0.09, 1.9)], PROBE, lw=2.0)
    ax.add_patch(Ellipse((fx, 1.45), 0.9, 0.22, fc=SKY, alpha=0.6, ec=PROBE, lw=1.0, zorder=4))
    ax.text(fx - 0.6, 1.45, 'lens', ha='right', va='center', fontsize=12, color=GREY)
    # dichroic: lattice beam 1 from below-right
    xd = 6.0
    mirror(ax, xd, yb, 135, color=LATT)
    beam(ax, [(xd, 2.55), (xd, yb), (fx, yb)], LATT, lw=6.0, alpha=0.35)
    box(ax, xd, 2.95, 3.4, 0.62, '1064 nm lattice beam 1', '#F6E3DF', LATT)
    # PBS combining coupling (straight) and probe (from SLM)
    xp = 9.0
    ax.add_patch(Rectangle((xp - 0.3, yb - 0.3), 0.6, 0.6, fc='#E8F1F8', ec=INK, lw=1.4, zorder=4))
    ax.plot([xp - 0.3, xp + 0.3], [yb - 0.3, yb + 0.3], color=INK, lw=1.2, zorder=5)
    ax.text(xp, yb - 0.6, 'PBS', ha='center', va='top', fontsize=12, color=GREY)
    beam(ax, [(13.3, yb), (xp, yb), (fx, yb)], COUPL, lw=3.0, z=2)
    beam(ax, [(xp + 0.09, 3.95), (xp + 0.09, yb + 0.06), (fx, yb + 0.06)], PROBE, lw=2.0, z=2)
    box(ax, 14.6, yb, 2.6, 0.75, 'coupling 795 nm\nLP$_{01}$, $\\Omega_c$', '#FBE9DE', COUPL)
    # probe arm: laser -> SLM -> mirror down
    box(ax, 14.6, 3.95, 2.6, 0.75, 'probe 795 nm\n(+6.8 GHz, locked)', '#E1EEF7', PROBE)
    ax.add_patch(Rectangle((11.0, 3.65), 0.32, 0.6, fc='#3A4A63', ec=INK, zorder=4))
    for k in range(4):
        ax.add_patch(Rectangle((11.0, 3.65 + 0.15 * k), 0.32, 0.15,
                               fc='#3A4A63' if k % 2 == 0 else '#7D8CA6', ec='none', zorder=5))
    ax.text(11.16, 4.45, 'SLM\n0/$\\pi$ phase', ha='center', va='bottom', fontsize=12)
    beam(ax, [(13.3, 3.95), (11.32, 3.95)], PROBE, lw=2.4)
    beam(ax, [(11.0, 3.95), (xp + 0.09, 3.95)], PROBE, lw=2.4)
    mirror(ax, xp + 0.09, 3.95, 135)
    # mode inset next to the probe arm (LP11 two lobes)
    for sgn in (-1, 1):
        ax.add_patch(Ellipse((10.0 + sgn * 0.22, 3.35), 0.3, 0.42, fc=PROBE, alpha=0.75, ec='none', zorder=4))
    ax.text(10.0, 2.85, 'LP$_{11}$', ha='center', va='top', fontsize=12, color=PROBE)
    # cross-section inset of the fibre core
    ax.add_patch(FancyBboxPatch((11.6, 5.55), 4.0, 2.9, boxstyle='round,pad=0.02,rounding_size=0.15',
                                fc=LIGHT, ec='none', zorder=1))
    ax.text(13.6, 8.15, 'core cross-section', ha='center', fontsize=13, fontweight='bold')
    n = 201
    Xg, Yg = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    Rg = np.hypot(Xg, Yg)
    I01 = np.where(Rg <= 1, jv(0, E.U01 * Rg)**2, np.nan)
    I11 = np.where(Rg <= 1, (jv(1, E.U1[0] * Rg) * Xg / np.maximum(Rg, 1e-9))**2, np.nan)
    for cx, I, cmap, lab in ((12.6, I01, 'Oranges', 'LP$_{01}$'), (14.6, I11, 'Blues', 'LP$_{11}$')):
        ax.imshow(I, extent=(cx - 0.8, cx + 0.8, 5.95, 7.55), cmap=cmap, origin='lower', zorder=2)
        ax.add_patch(Circle((cx, 6.75), 0.8, fill=False, ec=GREY, lw=1.5, zorder=3))
        ax.add_patch(Circle((cx, 6.75), 0.07, fc=ATOM, ec=INK, lw=0.8, zorder=4))
        ax.text(cx, 5.72, lab, ha='center', fontsize=13)
    # legend of beams
    hs = [Line2D([], [], color=COUPL, lw=3), Line2D([], [], color=PROBE, lw=3),
          Line2D([], [], color=LATT, lw=5, alpha=0.5), Line2D([], [], marker='o', ls='', color=ATOM, ms=8)]
    ax.legend(hs, ['coupling (LP$_{01}$)', 'probe (LP$_{11}$)', 'lattice 1064 nm (LP$_{01}$)', 'atoms'],
              loc='lower left', bbox_to_anchor=(0.665, 0.17), fontsize=13, frameon=False)
    save(fig, 'setup.png')


# ------------------------------------------------------------- result plots
def saturation():
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.6))
    for key, c in (('A', PROBE), ('B', COUPL)):
        d = load(f'sat_{key}.dat')
        ax[0].loglog(d['s'], d['WME'], 'o-', color=c)
        ax[0].loglog(d['s'], d['WSM2'], '--', color=c, lw=1.6)
        ax[1].loglog(d['s'], d['nME'], 'o-', color=c)
        ax[1].loglog(d['s'], d['nSM2'], '--', color=c, lw=1.6)
    ax[0].set_ylabel(r'cooling rate $W$ (s$^{-1}$)')
    ax[1].set_ylabel(r'steady state $\langle n_x\rangle$')
    for a in ax:
        a.set_xlabel('saturation parameter $s$')
        a.axvline(1, color=GREY, lw=1, ls=':')
    h = [Line2D([], [], color=PROBE, lw=3), Line2D([], [], color=COUPL, lw=3),
         Line2D([], [], color=INK, marker='o'), Line2D([], [], color=INK, ls='--', lw=1.6)]
    fig.legend(h, [r'A: $a=15\,\mu$m', r'B: $a=3.5\,\mu$m', 'master equation', 'simple model'],
               loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout()
    save(fig, 'saturation.png')


ALPHA = 0.397            # recoil fraction per radial axis


def hscale():
    """Where systems sit on the recoil parameter h, with the thresholds 1/3, 1/2, 1."""
    fig, ax = plt.subplots(figsize=(11.5, 3.3))
    edges = [3e-3, 1 / 3, 0.5, 1, 3]
    cols = ['#E8F1F8', '#FFF4E0', '#FBE3D6', '#EBD5D1']
    labs = ['near-thermal tail', 'heavy tail', 'mean set by cut-offs', 'no ground-state cooling']
    for x0, x1, c in zip(edges[:-1], edges[1:], cols):
        ax.axvspan(x0, x1, color=c, lw=0)
    ax.text(0.04, 1.75, labs[0], ha='center', va='center', fontsize=13, color=GREY)
    ax.text(0.41, 1.5, labs[1], ha='center', va='center', fontsize=12.5, color=GREY, rotation=90)
    ax.text(0.71, 1.38, 'mean set\nby cut-offs', ha='center', va='center', fontsize=12, color=GREY, rotation=90)
    ax.text(1.73, 1.15, 'no ground-\nstate cooling', ha='center', va='center', fontsize=12.5, color=GREY)
    for x in (1 / 3, 0.5, 1):
        ax.axvline(x, color=INK, lw=1.0, ls=':')
    ax.text(1 / 3, -0.55, '1/3', ha='center', fontsize=12)
    ax.text(0.5, -0.55, '1/2', ha='center', fontsize=12)
    ax.text(1, -0.55, '1', ha='center', fontsize=12)
    pts = [(0.008, 'trapped ion\n(Ca$^+$, 3 MHz)', GREY, 0.75), (0.03, 'Rb tweezer\n(100 kHz)', GREY, 0.75),
           (0.094, 'fibre B\n$a=3.5\\,\\mu$m', COUPL, 0.75), (0.24, 'fibre A\n$b_p=5/6$', PROBE, 0.75),
           (0.40, 'fibre A\n$b_p=1/2$', PROBE, 0.75)]
    for x, lab, c, y in pts:
        ax.plot(x, y, 'o', ms=11, color=c, mec=INK, mew=0.8, zorder=5)
        if x > 0.3:
            ax.text(x * 1.07, y - 0.32, lab, ha='left', va='top', fontsize=12, color=INK,
                    bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
        else:
            ax.text(x, y - 0.32, lab, ha='center', va='top', fontsize=12, color=INK)
    ax.set_xscale('log')
    ax.set_xlim(3e-3, 3)
    ax.set_ylim(-0.75, 2.15)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(r'recoil per removed phonon, $h$   ($b_p=1/2$ for ion and tweezer)')
    save(fig, 'hscale.png')


def tail():
    """Steady-state Fock distribution: power law against a thermal state of equal <n>."""
    d = load('convergence/tail.dat')
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    n = d['n']
    for col, c, lab in (('P_A_s0p1', PROBE, 'case A, $s=0.1$'), ('P_A_s1', SKY, 'case A, $s=1$'),
                        ('P_B_s1', COUPL, 'case B, $s=1$')):
        P = d[col]
        ok = P > 0
        ax.loglog(n[ok], P[ok], 'o-', color=c, ms=4, lw=1.6, label=lab)
    P = d['P_A_s0p1']
    nm = np.sum(n * P) / np.sum(P)
    th = (nm / (1 + nm))**n / (1 + nm)
    ax.loglog(n, th, '--', color=GREY, lw=1.6, label=r'thermal, same $\langle n\rangle$')
    h = ALPHA / (0.5 * E.paper_cases()['A']['nu_wrec'])
    x = np.array([5, 80])
    ref = np.interp(5, n, P)
    ax.loglog(x, ref * (x / 5)**(-1 / h) * 2.2, '-', color=INK, lw=1.2, label=r'$n^{-1/h}$, $h=%.2f$ (case A)' % h)
    ax.set_ylim(1e-9, 0.2)
    ax.set_xlim(0.9, 300)
    ax.set_xlabel('Fock state $n$')
    ax.set_ylabel('steady-state $P(n)$')
    ax.legend(loc='upper right', fontsize=12)
    save(fig, 'tail.png')


def regimes():
    d = load('convergence/rates_vs_n.dat')
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    x0, x1 = 3e-3, 400
    ax.axvspan(x0, 1, color='#E8F1F8', lw=0)
    ax.axvspan(1, 18, color='#F2F2F2', lw=0)
    ax.axvspan(18, x1, color='#FBE9DE', lw=0)
    for x, lab in ((0.06, 'resolved\nsideband'), (4.2, 'saturated'), (85, 'fast channel')):
        ax.text(x, 1.3e5, lab, ha='center', va='top', fontsize=13, color=GREY)
    for s, c in ((0.01, GREY), (0.1, SKY), (1.0, PROBE)):
        r = d[d['s'] == s]
        ax.loglog(s * r['n0'], -r['drift'], 'o-', color=c, ms=5, label=f'$s={s:g}$')
    ax.set_xlim(x0, x1)
    ax.set_ylim(20, 2e5)
    ax.set_xlabel(r'$s\,n$  (saturation of the $n\to n-1$ sideband)')
    ax.set_ylabel(r'phonon removal rate (s$^{-1}$), case A')
    ax.legend(loc='lower right', fontsize=12)
    save(fig, 'regimes.png')


def cooling_map():
    D = load('universal_map.dat')
    NU, BP = np.unique(D['nu_wrec']), np.unique(D['bp'])

    def grid(q):
        Z = np.full((BP.size, NU.size), np.nan)
        for i, b in enumerate(BP):
            for j, v in enumerate(NU):
                r = D[(D['bp'] == b) & (D['nu_wrec'] == v) & (D['DeltaG'] == 3) & (D['s'] == 1)]
                r = r[np.argsort(r['N'])][-1]
                Z[i, j] = q(r)
        return Z
    en = np.sqrt(NU[1:] * NU[:-1])
    en = np.r_[NU[0]**2 / en[0], en, NU[-1]**2 / en[-1]]
    eb = (BP[1:] + BP[:-1]) / 2
    eb = np.r_[2 * BP[0] - eb[0], eb, 2 * BP[-1] - eb[-1]]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.4))
    Z1 = grid(lambda r: r['P0'])
    Z2 = grid(lambda r: r['W0'] if (r['P0'] > 0.85 and r['n'] < 1) else np.nan)
    from matplotlib.colors import LogNorm
    m1 = ax[0].pcolormesh(en, eb, Z1, cmap='Blues', vmin=0.5, vmax=1)
    m2 = ax[1].pcolormesh(en, eb, Z2, cmap='Blues', norm=LogNorm(0.008, 0.12))
    fig.colorbar(m1, ax=ax[0], label='ground-state population $P_0$')
    fig.colorbar(m2, ax=ax[1], label=r'cooling rate $W_0/\nu_r$')
    x = np.geomspace(0.7, 20, 200)
    for a in ax:
        for hh, ls in ((1 / 3, ':'), (0.5, '--'), (1, '-')):
            a.plot(x, ALPHA / (hh * x), color=INK, ls=ls, lw=1.4)
        for xv, bv, mk in ((NU[2], 0.5, 'o'), (NU[2], 5 / 6, 's'), (NU[5], 0.5, '^')):
            a.plot(xv, bv, mk, mfc='white', mec=INK, ms=10, mew=1.4)
        a.set_xscale('log')
        a.set_xlim(en[0], en[-1])
        a.set_ylim(eb[0], eb[-1])
        a.set_yticks(BP, ['1/6', '1/3', '1/2', '2/3', '5/6'])
        a.set_xticks([1, 2, 5, 10], ['1', '2', '5', '10'])
        a.set_xlabel(r'trap stiffness $\nu_r/\omega_{\rm rec}$')
        a.set_ylabel('branching to probe ground state $b_p$')
    hs = [Line2D([], [], color=INK, ls=':'), Line2D([], [], color=INK, ls='--'), Line2D([], [], color=INK, ls='-'),
          Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=9), Line2D([], [], marker='s', ls='', mfc='white', mec=INK, ms=9),
          Line2D([], [], marker='^', ls='', mfc='white', mec=INK, ms=9)]
    fig.legend(hs, ['$h=1/3$', '$h=1/2$', '$h=1$', "A, D1 F'=2", "A, D1 F'=1", 'B'], loc='upper center',
               ncol=6, bbox_to_anchor=(0.5, 1.07), fontsize=12)
    fig.tight_layout(w_pad=2)
    save(fig, 'map.png')


def two_axes():
    d = load('two_axes_steady.dat')
    d = d[d['s'] == 1]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for DG, c in ((3, PROBE), (6, COUPL)):
        r = d[d['DeltaG'] == DG]
        o = np.argsort(r['h'])
        ax.semilogy(r['h'][o], 1 - r['P0_1D'][o], 'o', mfc='white', color=c, ms=7, mew=1.5)
        ax.semilogy(r['h'][o], 1 - r['P00_2D'][o], 'o', color=c, ms=7)
    for x in (1 / 3, 0.5):
        ax.axvline(x, color=GREY, ls=':', lw=1.2)
    ax.set_xlabel('recoil parameter $h$ (one axis)')
    ax.set_ylabel('population outside the ground state')
    ax.set_ylim(2e-3, 1.2)
    hs = [Line2D([], [], marker='o', ls='', color=PROBE, ms=7), Line2D([], [], marker='o', ls='', color=COUPL, ms=7),
          Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=7), Line2D([], [], marker='o', ls='', color=INK, ms=7)]
    ax.legend(hs, [r'$\Delta=3\Gamma$', r'$\Delta=6\Gamma$', r'one axis, $1-P_0$', r'both axes, $1-P_{00}$'],
              loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=12, ncol=4, handletextpad=0.3, columnspacing=1.0)
    save(fig, 'two_axes.png')


def capture():
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for key, c in (('A', PROBE), ('B', COUPL)):
        m = load(f'temperature_{key}.dat')
        ax.loglog(m['TuK'], m['t90_ms'], 'o-', color=c, ms=5, label=f'case {key}')
    ax.set_xlabel(r'initial radial temperature ($\mu$K)')
    ax.set_ylabel(r'time to 90 % in $n=0$ (ms)')
    ax.set_ylim(0.1, 1000)
    ax.legend(loc='upper left', fontsize=12)
    save(fig, 'capture.png')


def memory_panel():
    sys.path.insert(0, str(ROOT / 'notes' / 'quantum-memory'))
    import memory_budget as mb
    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.4))
    a = np.geomspace(2.5, 25, 100)
    sig = 3 * mb.LAM_S**2 / (2 * np.pi) * 1e12
    rho = 1e4 / 0.532 / (mb.a_eff(a * 1e-6) * 1e12)          # um^-2 per cm, one atom per site
    for r, c, lab in ((0.25, PROBE, 'clock Λ ($r=1/4$)'), (0.5, COUPL, 'stretched Λ ($r=1/2$)')):
        ax[0].loglog(a, sig * r * rho, color=c, label=lab)
    for od, lab in ((93, r'$\eta=0.90$'), (201, r'$\eta=0.95$')):
        ax[0].axhline(od, color=GREY, ls='--', lw=1.2)
        ax[0].text(24, od * 1.08, lab, ha='right', va='bottom', fontsize=12, color=GREY)
    ax[0].set_xlabel(r'core radius $a$ ($\mu$m)')
    ax[0].set_ylabel('optical depth per cm')
    ax[0].set_xticks([3, 5, 10, 20], ['3', '5', '10', '20'])
    ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].set_ylim(2, 600)
    ax[0].legend(loc='lower left', fontsize=12)
    ax[0].set_title('one atom per site; OD = 0.302 r µm² × ρ$_s$', fontsize=12, loc='left')
    T = np.geomspace(1e-3, 100, 300)
    scen = (('full depth, no echo', dict(u=1.0), GREY), ('full depth, echo', dict(u=1.0, echo=True), SKY),
            (r'depth ×0.01, echo, $\tau_{\rm bg}$ = 2 s', dict(u=.01, echo=True), PROBE),
            (r'depth ×0.01, echo, $\tau_{\rm bg}$ = 100 s', dict(u=.01, echo=True, tau_bg=100.), COUPL))
    for lab, kw, c in scen:
        s = mb.storage_terms(T, **kw)
        ax[1].semilogx(T, s['heat'] * s['inh'] * s['mag'] * s['bg'], color=c, label=lab)
    ax[1].set_xlabel('storage time (s)')
    ax[1].set_ylabel(r'$\eta(T)/\eta(0)$')
    ax[1].set_ylim(0, 1.05)
    ax[1].legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=11)
    ax[1].set_title('1 % lattice-intensity spread along the column', fontsize=12, loc='left')
    fig.tight_layout(w_pad=3)
    save(fig, 'memory.png')


def rb_figure():
    """Rate against ground-state population: hyperfine Lambda configurations (Delta scan) and the
    cost of the Zeeman structure (remixing schemes), case A."""
    om = E.paper_cases()['A']['om']
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 8.2), sharex=True)
    d = load('lines_delta_A.dat')
    cfgs = (((1, 1, 2), SKY, "D1 F'=1, probe F=2"), ((1, 2, 1), PROBE, "D1 F'=2"),
            ((2, 2, 1), COUPL, "D2 F'=2"), ((2, 1, 1), ATOM, "D2 F'=1, probe F=1"))
    for (line, Fe, gp), c, lab in cfgs:
        r = d[(d['line'] == line) & (d['Fe'] == Fe) & (d['gp'] == gp) & (d['okopt'] == 1)]
        r = r[np.argsort(r['DeltaG'])]
        ax[0].loglog(1 - r['P0opt'], r['W0opt'], 'o-', color=c, ms=6, label=lab)
        if (line, Fe) == (1, 1):
            for x, y, D in zip(1 - r['P0opt'], r['W0opt'], r['DeltaG']):
                ax[0].annotate(r'%g$\Gamma$' % D, (x, y), textcoords='offset points', xytext=(-4, 7), fontsize=11, ha='right', color=GREY)
    ax[0].set_ylabel(r'ground-state rate $W_0$ (s$^{-1}$)')
    ax[0].set_title(r'a   all excited hyperfine levels, $\Delta$ = 2–12$\Gamma$', loc='left', fontsize=13)
    ax[0].legend(loc='upper left', fontsize=11)
    z = np.genfromtxt(RESULTS / 'zeeman_refine.dat', names=True)
    z = z[(z['P0'] > 0.85) & (z['n'] < 1) & (z['W0'] > 0) & np.isfinite(z['W0']) & (np.abs(z['nu_wrec'] - 1.97) < 0.01)]
    m = load('universal_map.dat')
    ref = m[(np.abs(m['nu_wrec'] - 1.968) < 0.01) & (np.abs(m['bp'] - 5 / 6) < 0.01) & (m['P0'] > 0.85)]
    groups = ((ref, INK, r"closed $\Lambda$, $b_p$ = 5/6 (no Zeeman)"),
              (z[(z['scheme'] == 0) & (z['kind'] == 2)], SKY, "D1 F'=1 + repumper"),
              (z[(z['scheme'] == 0) & (z['kind'] == 0)], GREEN, "D1 F'=1 + transverse field"),
              (z[z['scheme'] == 1], COUPL, r"Zeeman $\Lambda$ in F=2 + repumper"))
    for g, c, lab in groups:
        x, y = 1 - g['P0'], g['W0'] * om
        o = np.argsort(x)
        ax[1].loglog(x, y, '.', color=c, ms=4, alpha=0.35)
        ax[1].loglog(x[o], np.maximum.accumulate(y[o]), '-', color=c, lw=2.2, label=lab)
    ax[1].axvline(0.03, color=GREY, ls=':', lw=1.2)
    ax[1].set_xlabel(r'population outside the ground state, $1-P_0$')
    ax[1].set_ylabel(r'ground-state rate $W_0$ (s$^{-1}$)')
    ax[1].set_title('b   all 8 ground sublevels: best rate at each $P_0$ (dotted: $P_0$ = 0.97)', loc='left', fontsize=12)
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=10.5)
    ax[1].set_xlim(1.5e-3, 0.16)
    ax[1].set_ylim(50, 1.5e4)
    ax[0].set_ylim(300, 3e4)
    fig.tight_layout(h_pad=1.5)
    save(fig, 'rb.png')


def robust_figure():
    """Rate against probe-mode impurity / trap offset, and against trap-frequency mismatch."""
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 7.6))
    c = load('carrier_A.dat')
    eta = 0.021
    e, x = c['kind'] == 0, c['kind'] == 1
    w0 = c['WmaxME'][e][0]
    ax[0].semilogx(c['epseff'][e][1:] / eta, c['WmaxME'][e][1:] / w0, 'o-', color=PROBE, ms=6, label=r'LP$_{01}$ admixture $\epsilon$')
    ax[0].semilogx(c['epseff'][x] / eta, c['WmaxME'][x] / w0, 's', color=COUPL, ms=9, mfc='none', mew=1.8, label=r'trap offset $x_s$ ($\epsilon = \eta_\perp x_s/x_0$)')
    ax[0].semilogx(c['epseff'][e][1:] / eta, c['Wunified'][e][1:] / w0, '--', color=GREY, lw=1.6, label='simple model')
    ax[0].axhline(0.5, color=GREY, lw=1, ls=':')
    ax[0].set_ylim(0, 1.05)
    ax[0].set_xlabel(r'impurity $\epsilon/\eta_\perp$')
    ax[0].set_ylabel(r'maximum rate / clean node')
    ax[0].set_title('a   probe purity, case A', loc='left', fontsize=13)
    ax[0].legend(loc='lower left', fontsize=11)
    for key, col in (('A', PROBE), ('B', COUPL)):
        m = load(f'mismatch_{key}.dat')
        ax[1].plot(m['nurel'], m['WME'] / m['WME'].max(), 'o-', color=col, ms=5, label=f'case {key}')
    ax[1].axvspan(0.9, 1.1, color=GREEN, alpha=0.12, lw=0)
    ax[1].text(1.0, 0.08, '±10 %', ha='center', fontsize=12, color=GREEN)
    ax[1].set_ylim(0, 1.05)
    ax[1].set_xlabel(r'radial frequency / design value, $\nu_r/\nu_r^{\rm design}$')
    ax[1].set_ylabel(r'rate / maximum')
    ax[1].set_title('b   lasers fixed, trap frequency varies', loc='left', fontsize=13)
    ax[1].legend(loc='upper left', fontsize=11)
    fig.tight_layout(h_pad=2.0)
    save(fig, 'robust.png')


if __name__ == '__main__':
    for f in (equations, mode_maps, lattice_sketch, levels, fano_plot, setup, saturation,
              hscale, tail, regimes, cooling_map, two_axes, capture, memory_panel, rb_figure, robust_figure):
        f()
        print('ok', f.__name__)
