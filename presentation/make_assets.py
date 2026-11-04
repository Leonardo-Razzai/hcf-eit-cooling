"""Images for the slide deck: schematics, equations and result plots.

Reads simulation/results and the eitsim package; writes presentation/assets/*.png
(300 dpi). Run from anywhere:  python presentation/make_assets.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, Polygon
from matplotlib.lines import Line2D
from scipy.special import jv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / 'simulation' / 'results'
OUT = HERE / 'assets'
sys.path.insert(0, str(ROOT / 'simulation'))
import eitsim as E  # noqa: E402

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
    'eq_gradient': r'$\Omega_p(\hat x)\simeq\bar\Omega_p\,G_{11}\,\dfrac{\hat x}{a}=\Omega_\perp\,(b+b^\dagger)$',
    'eq_eta': r'$\eta_\perp=G_{11}\,\dfrac{x_0}{a},\qquad G_{11}=3.49$',
    'eq_fano': r'$R(\delta)=\dfrac{\Omega_p^2\,\Gamma\,\delta^2}{\Gamma^2\delta^2+4\left[(\Delta+\delta)\delta-\Omega_c^2/4\right]^2}'
               r',\qquad \Omega_c^2=4\nu_r(\Delta+\nu_r)$',
    'eq_W': r'$W\simeq(1-h)\,\dfrac{\Gamma_n}{2}\,\dfrac{s}{1+s},\qquad \Gamma_n\simeq\dfrac{\Gamma\,\nu_r}{\Delta}$',
    'eq_n': r'$n_{\rm ss}\simeq\left(\dfrac{\Gamma}{4\Delta}\right)^2\dfrac{1+h}{1-h}\,(1+s),\qquad h=\kappa\,\alpha\,\eta_{\rm sp}^2,\ \ \kappa=1/b_p$',
    'eq_s': r'$s=\dfrac{2\Omega_\perp^2}{\Gamma\,\Gamma_n}=2\left(\dfrac{\Omega_R}{\Gamma_n}\right)^2,\qquad \Omega_R=\dfrac{\Omega_c\,\Omega_\perp}{2\Delta}$',
    'eq_lightshift': r'$\delta_{\rm diff}(r)\simeq\delta_0\left(1-\dfrac{V_\perp(r)}{U_0}\right),\qquad \delta_0\approx\dfrac{U_0\,\omega_{\rm hfs}}{\Delta_{\rm eff}}\approx 0.8\ \mathrm{kHz}$',
    'eq_spread': r'$\sigma_\delta\simeq\delta_0\,\dfrac{\hbar\omega_r\,\langle n_r\rangle}{U_0},\qquad T_2^*\approx\dfrac{1}{\pi\,\sigma_\delta}$',
    'eq_spinwave': r'$\rho_{pc}(z)\propto e^{i\Delta k\,z},\qquad \Delta k=\dfrac{\omega_{\rm hfs}}{c}\ \Rightarrow\ \lambda_{\rm sw}\approx 4.4\ \mathrm{cm}$',
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


def robust():
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
    d = load('carrier_A.dat')
    e, x = d['kind'] == 0, d['kind'] == 1
    w0 = d['WmaxME'][e][d['epseff'][e] == 0][0]
    ax[0].plot(d['epseff'][e][1:], d['WmaxME'][e][1:] / 1e3, 'o-', color=PROBE, label=r'LP$_{01}$ admixture $\epsilon$')
    ax[0].plot(d['epseff'][x], d['WmaxME'][x] / 1e3, 's', color=COUPL, ms=8, mfc='none', mew=1.8, label='trap offset $x_s$')
    ax[0].plot(d['epseff'][e][1:], d['Wunified'][e][1:] / 1e3, '--', color=GREY, lw=1.6, label='simple model')
    ax[0].axhline(w0 / 1e3, color=PROBE, lw=1, ls=':')
    ax[0].set_xscale('log')
    ax[0].set_ylim(0, 3.2)
    ax[0].set_xlabel(r'$\epsilon_{\rm eff}$')
    ax[0].set_ylabel(r'max. rate ($10^3$ s$^{-1}$)')
    ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=11)
    t = load('twophoton_A.dat')
    ax[1].semilogy(t['d2'], t['n0'], 'o-', color=PROBE, label=r'clean node, $\epsilon=0$')
    ax[1].semilogy(t['d2'], t['n5'], 's-', color=COUPL, label=r'$\epsilon=0.05$')
    ax[1].set_xlabel(r'two-photon error $\delta_{\rm e}/\nu_r$')
    ax[1].set_ylabel(r'$\langle n_x\rangle$')
    ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=11)
    for key, c in (('A', PROBE), ('B', COUPL)):
        m = load(f'mismatch_{key}.dat')
        ax[2].plot(m['nurel'], m['WME'] / m['WME'].max(), 'o-', color=c, ms=4, label=f'case {key}')
    ax[2].axvspan(0.9, 1.1, color=GREEN, alpha=0.12, lw=0)
    ax[2].set_xlabel(r'$\nu_r/\nu_r^{\rm design}$')
    ax[2].set_ylabel(r'$W/W_{\max}$')
    ax[2].legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=11)
    fig.tight_layout(w_pad=2.5)
    save(fig, 'robust.png')


def dynamics():
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.6))
    d = load('dyn_A.dat')
    ax[0].plot(d['tms'], d['nME'], 'o', color=PROBE, ms=5, label='master equation')
    ax[0].plot(d['tms'], d['nGRM'], '-', color=COUPL, label='rate model, calibrated')
    ax[0].plot(d['tms'], d['nSM'], '--', color=GREY, lw=1.6, label='simple model')
    ax[0].set_yscale('log')
    ax[0].set_ylim(0.015, 10)
    ax[0].set_xlabel('time (ms)')
    ax[0].set_ylabel(r'$\langle n_x\rangle$')
    ax[0].set_title(r'case A, from $\langle n\rangle=2$', fontsize=14, loc='left')
    ax[0].legend(loc='upper right', fontsize=11)
    for key, c in (('A', PROBE), ('B', COUPL)):
        cp = load(f'capture_{key}.dat')
        ax[1].semilogy(cp['tms'], cp['n_staged'], '-', color=c, label=f'case {key} (staged)')
    ax[1].set_xscale('log')
    ax[1].set_xlim(0.05, 150)
    ax[1].set_ylim(0.01, 100)
    ax[1].set_xlabel('time (ms)')
    ax[1].set_ylabel(r'$\langle n_x\rangle$')
    ax[1].set_title(r'cooling from 10 $\mu$K', fontsize=14, loc='left')
    ax[1].legend(loc='upper right', fontsize=11)
    fig.tight_layout(w_pad=3)
    save(fig, 'dynamics.png')


def temperature():
    """Cooling versus initial radial temperature k_B T_r / U0 (calibrated GRM)."""
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.6))
    tr = load('temptrace_A.dat')
    TU = load('temperature_A.dat')['TU0']
    pick = (0, 5, 10, 15, 19)
    cols = plt.cm.viridis(np.linspace(0.0, 0.85, len(pick)))[::-1]
    for j, c in zip(pick, cols):
        ax[0].semilogy(tr['tms'], tr['n%02d' % j], color=c, lw=2, label='%.3g' % TU[j])
    ax[0].set_xlim(0, 400)
    ax[0].set_ylim(0.1, 1e4)             # top kept empty for the legend
    ax[0].set_xlabel('time (ms)')
    ax[0].set_ylabel(r'$\langle n_x\rangle$')
    ax[0].set_title('case A', fontsize=14, loc='left')
    ax[0].legend(title=r'$k_BT_r/U_0$', loc='upper right', fontsize=10, title_fontsize=10, ncol=3)
    for key, c in (('A', PROBE), ('B', COUPL)):
        m = load(f'temperature_{key}.dat')
        ax[1].loglog(m['TU0'], m['t90_ms'], 'o-', color=c, ms=5, label=f'case {key}')
    ax[1].set_xlim(1e-3, 1e-1)
    ax[1].set_ylim(0.05, 3000)
    ax[1].set_xlabel(r'initial $k_BT_r/U_0$')
    ax[1].set_ylabel(r'time to 90 % in $n=0$ (ms)')
    ax[1].legend(loc='upper left', fontsize=11)
    top = ax[1].secondary_xaxis('top', functions=(lambda x: x * 533.0, lambda t: t / 533.0))
    top.set_xlabel(r'$T_r$ ($\mu$K)')
    fig.tight_layout(w_pad=3)
    save(fig, 'temperature.png')


def memory():
    """Thermal spread of the differential light shift vs radial occupation."""
    h = 6.62607015e-34
    kl = 2 * np.pi / 1064e-9
    U0 = E.M_RB * (2 * np.pi * 300e3)**2 / (2 * kl**2)          # J
    delta0 = 800.0                                               # Hz
    n = np.logspace(-2.3, 2, 200)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    pts = {}                                   # (10 uK, end of staged capture)
    for key in 'AB':
        cp = load(f'capture_{key}.dat')
        pts[key] = (cp['n_staged'][0], cp['n_staged'][-1])
    for key, c in (('A', PROBE), ('B', COUPL)):
        fc = E.paper_cases()[key]
        sig = delta0 * h * fc['nur'] * n / U0
        ax.loglog(n, sig, color=c, label=f'case {key} ($\\nu_r$ = {fc["nur"] / 1e3:.1f} kHz)')
        for nn, mk, lab in ((pts[key][0], 'o', '10 $\\mu$K'), (pts[key][1], '*', 'after EIT')):
            s = delta0 * h * fc['nur'] * nn / U0
            ax.plot(nn, s, mk, color=c, ms=11 if mk == '*' else 8, mec=INK, mew=0.8, zorder=5)
    ax.set_xlabel(r'radial occupation $\langle n_r\rangle$ (each axis)')
    ax.set_ylabel(r'light-shift spread $\sigma_\delta$ (Hz)')
    ax.set_ylim(3e-3, 300)
    sec = ax.secondary_yaxis('right', functions=(lambda s: 1 / (np.pi * s), lambda t: 1 / (np.pi * t)))
    sec.set_ylabel(r'$T_2^*\approx 1/\pi\sigma_\delta$ (s)')
    hs, ls = ax.get_legend_handles_labels()
    hs += [Line2D([], [], marker='o', ls='', color=GREY, mec=INK, ms=8),
           Line2D([], [], marker='*', ls='', color=GREY, mec=INK, ms=11)]
    ls += [r'thermal, 10 $\mu$K', 'after radial EIT cooling']
    ax.legend(hs, ls, loc='upper left', fontsize=12)
    save(fig, 'memory.png')


if __name__ == '__main__':
    for f in (equations, mode_maps, lattice_sketch, levels, fano_plot, setup,
              saturation, robust, dynamics, temperature, memory):
        f()
        print('ok', f.__name__)
