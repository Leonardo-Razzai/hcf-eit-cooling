#!/usr/bin/env python3
"""Analytic and very light numerical budget for a spin-wave memory in a hollow-core fibre.

Self-contained on purpose: it needs numpy, scipy, sympy and matplotlib only, and does not
import the cooling code. Runs in about a minute. Writes budget_output.txt and memory_budget.png
next to this file. Every number in memory_theory.md tagged [script] comes from here.

Conventions (Gorshkov et al., PRA 76, 033805):
  d       resonant intensity attenuation is exp(-2 d), so d = OD / 2
  gamma   polarisation decay, Gamma / 2
  h       int Omega^2 dt / gamma, with Omega half the Rabi frequency
"""
import os
import numpy as np
import scipy.constants as sc
from scipy.linalg import eigh, expm
from scipy.special import i0e, i1e, jv, jn_zeros
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_6j

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = None   # opened in __main__, so that crosstalk.py can import this module


def say(*args):
    line = ' '.join(str(a) for a in args)
    print(line)
    if OUT:
        OUT.write(line + '\n')


# ---------------------------------------------------------------- constants
# Fundamental constants: scipy.constants (CODATA). Rb data: Steck, Rubidium 87 D Line Data rev. 2.3.4,
# except the D1 width, kept at the repository value (5.746 MHz).
HBAR, H, C, KB = sc.hbar, sc.h, sc.c, sc.k
M = 86.909180531 * sc.atomic_mass   # 87Rb
LAM_S = 794.978851156e-9            # D1
LAM_D2 = 780.241209686e-9
GAMMA = 2 * np.pi * 5.746e6
GAMMA_D2 = 2 * np.pi * 6.0666e6
LAM_L = 1064e-9
HFS = 6.834682610904e9
U01 = jn_zeros(0, 1)[0]
U02 = jn_zeros(0, 2)[1]
WREC_L = HBAR * (2 * np.pi / LAM_L)**2 / (2 * M)   # lattice recoil, rad/s
ALPHA = 0.397              # fraction of emission recoil per radial axis (repo value)
G_EARTH = sc.g


def a_eff(a, u=U01):
    """On-axis effective area of an LP0m mode: P / I(0)."""
    return np.pi * a**2 * jv(1, u)**2


# ---------------------------------------------------------------- 1. line strengths
def _raw(Fg, mg, Fe, me):
    I_, J_ = Rational(3, 2), Rational(1, 2)
    q = me - mg
    if abs(q) > 1 or abs(mg) > Fg or abs(me) > Fe:
        return 0.0
    return float((2 * Fe + 1) * wigner_6j(J_, J_, 1, Fe, Fg, I_)**2 * clebsch_gordan(Fe, 1, Fg, me, mg - me, mg)**2)


def branching(Fg, mg, Fe, me):
    """Decay fraction |F'm'> -> |F m> on D1, normalised over all F, m.

    sigma(F m -> F' m') = (3 lam^2 / 2 pi) r for D1 (J = J' = 1/2). Check: D1 F'=1 -> F=2 sums to 5/6.
    """
    tot = sum(_raw(F, m, Fe, me) for F in (1, 2) for m in range(-F, F + 1))
    return _raw(Fg, mg, Fe, me) / tot


def lambda_table():
    say('\n== 1. D1 lambda schemes: signal and control branching ratios r')
    say('check: sum r(F\'=1 m=1 -> F=2) = %.4f (5/6), r(|1\',1> -> |2,2>) = %.4f (1/2)' % (
        sum(branching(2, m, 1, 1) for m in range(-2, 3)), branching(2, 2, 1, 1)))
    say('quantisation axis along the fibre (axial B): only sigma+- are available, no pi')
    say('r = 1 is a closed two-level transition; OD scales with r_signal')
    say('%-20s %-14s %6s %6s %s' % ('g -> s', 'via, pol', 'r_sig', 'r_ctl', 'control on g? (FWM)'))
    rows = []
    for g, s, name in (((1, 0), (2, 0), 'clock'), ((2, 0), (1, 0), 'clock'),
                       ((1, -1), (2, 1), 'magic pair'), ((2, 1), (1, -1), 'magic pair'),
                       ((2, 2), (1, 0), 'non-clock')):
        for Fe in (1, 2):
            for me in range(-Fe, Fe + 1):
                if abs(me - g[1]) > 1 or abs(me - s[1]) > 1:
                    continue
                rs = branching(g[0], g[1], Fe, me)
                rc = branching(s[0], s[1], Fe, me)
                if rs == 0 or rc == 0:
                    continue
                qc = me - s[1]
                # control acting on g reaches level m_g + qc; forbidden means no Stokes channel
                ctl_on_g = sum(branching(g[0], g[1], f, g[1] + qc)
                               for f in (1, 2) if abs(g[1] + qc) <= f)
                say('%-20s F\'=%d m=%+d    %6.3f %6.3f   %.3f' % (
                    '|%d,%+d>->|%d,%+d>' % (g + s), Fe, me, rs, rc, ctl_on_g))
                rows.append((g, s, Fe, me, rs, rc, ctl_on_g))
    return rows


# ---------------------------------------------------------------- 2. retrieval kernel
def kernel(d, n=None):
    """Gauss-Legendre discretised K = sqrt(w) k_r sqrt(w), k_r of Gorshkov Eq. (15)."""
    n = n or int(100 + 25 * np.sqrt(d))   # converged to 1e-5 (tested d = 50-1500)
    x, w = np.polynomial.legendre.leggauss(n)
    z = (x + 1) / 2
    w = w / 2
    zz = d * np.sqrt(np.outer(z, z))
    # i0e(x) = exp(-x) I0(x); write the exponent to avoid overflow at large d
    k = (d / 2) * np.exp(-d * (z[:, None] + z[None, :]) / 2 + zz) * i0e(zz)
    sw = np.sqrt(w)
    return z, w, sw[:, None] * k * sw[None, :]


def retrieval(d):
    """Largest retrieval efficiency eta_r and its spin-wave."""
    z, w, K = kernel(d)
    vals, vecs = eigh(K)
    return vals[-1], z, w, vecs[:, -1], K


def total_efficiency(d, phi=0.0, forward=False):
    """Optimal storage + retrieval.

    Backward at phase mismatch phi = Delta k L, or forward with the spin wave flipped.
    eta = largest eigenvalue of K D K D^dagger (see memory_theory.md, Section 2).
    """
    z, w, K = kernel(d)
    if forward:
        # flip z -> 1 - z: Gauss-Legendre nodes are symmetric, so the flip is an index reversal
        D = np.eye(len(z))[::-1]
    else:
        D = np.diag(np.exp(1j * phi * z))
    M = K @ D @ K @ D.conj().T
    return np.max(np.linalg.eigvals(M).real)


def validate_kernel():
    say('\n== 2. Retrieval kernel checks')
    for d in (10, 50, 200):
        z, w, K = kernel(d)
        s = np.sqrt(w)              # flat spin wave
        eta = s @ K @ s
        exact = 1 - (i0e(d) + i1e(d))   # Gorshkov Eq. (30)
        say('d=%-4d flat spin wave eta %.5f, Eq.(30) %.5f' % (d, eta, exact))
    for d in (20, 50, 100, 200, 400):
        say('d=%-4d (1-eta_r) d = %.2f  (Gorshkov: about 2.9 at large d)' % (d, (1 - retrieval(d)[0]) * d))


def efficiency_table():
    say('\n== 3. Optimal efficiency against d = OD/2')
    say('%6s %6s %8s %10s %10s' % ('OD', 'd', 'eta_r', 'backward', 'forward'))
    tab = []
    for od in (10, 20, 30, 50, 80, 100, 150, 200, 300, 500, 1000):
        d = od / 2
        er = retrieval(d)[0]
        fo = total_efficiency(d, forward=True)
        say('%6d %6.1f %8.4f %10.4f %10.4f' % (od, d, er, er**2, fo))
        tab.append((od, er**2, fo))
    for goal in (0.5, 0.7, 0.8, 0.9, 0.95, 0.99):
        lo, hi = 2.0, 4000.0
        for _ in range(26):
            mid = np.sqrt(lo * hi)
            if retrieval(mid / 2)[0]**2 > goal:
                hi = mid
            else:
                lo = mid
        say('backward eta = %.2f needs OD = %.0f' % (goal, hi))
    return np.array(tab)


def geometry_table(r_sig):
    say('\n== 4. Optical depth of a column, single occupancy (OD = N sigma / A_eff)')
    say('sigma = r (3 lam^2 / 2 pi), r = %.3f; lattice period %.0f nm' % (r_sig, LAM_L / 2 * 1e9))
    sig0 = 3 * LAM_S**2 / (2 * np.pi)
    say('%6s %6s %6s %8s %8s %8s %8s' % ('a(um)', 'L(mm)', 'fill', 'N', 'OD(r=1)', 'OD', 'eta_b'))
    for a, L, f in ((15, 10, 1), (15, 10, .5), (15, 50, .5), (7.5, 10, .5), (5, 10, .5),
                    (3.5, 10, .5), (3.5, 3, .5), (3.5, 1, .5), (3.5, 10, 1)):
        N = f * L * 1e-3 / (LAM_L / 2)
        od1 = N * sig0 / a_eff(a * 1e-6)
        od = od1 * r_sig
        say('%6.1f %6.0f %6.1f %8.0f %8.1f %8.1f %8.3f' % (a, L, f, N, od1, od, retrieval(od / 2)[0]**2))


def mismatch_table():
    dk = 2 * np.pi * HFS / C
    say('\n== 5. Backward retrieval, hyperfine phase mismatch')
    say('Delta k L = 2 (omega_hf / c) L = %.2f rad per cm' % (2 * dk * 0.01))
    say('%6s %6s %8s %8s %8s %8s' % ('OD', 'L(cm)', 'phi', 'back', 'back(0)', 'forward'))
    for od in (50, 100, 300):
        for L in (0.3, 1, 3):
            phi = 2 * dk * L * 0.01
            say('%6d %6.1f %8.2f %8.3f %8.3f %8.3f' % (od, L, phi, total_efficiency(od / 2, phi),
                                                         retrieval(od / 2)[0]**2,
                                                         total_efficiency(od / 2, forward=True)))


# ---------------------------------------------------------------- 3. temporal capacity
def capacity_table():
    say('\n== 6. Temporal modes (spectrum of the kernel K): modes with total efficiency above a level')
    say('Each eigenvalue lambda_j of K is the retrieval efficiency of one axial spin wave;')
    say('with backward retrieval the total efficiency of that mode is lambda_j^2.')
    say('%6s %6s %8s %8s %8s %10s' % ('OD', 'd', 'N(0.5)', 'N(0.7)', 'N(0.9)', 'N(0.7)/sqrt(d)'))
    for od in (30, 100, 300, 1000, 3000):
        d = od / 2
        vals = eigh(kernel(d)[2], eigvals_only=True)
        n = [int(np.sum(vals**2 > t)) for t in (.5, .7, .9)]
        say('%6d %6.0f %8d %8d %8d %10.2f' % (od, d, n[0], n[1], n[2], n[1] / np.sqrt(d)))


# ---------------------------------------------------------------- 4. erasure of spectators
def volterra(kappa, n):
    """Matrix V with (V S)_i = int_0^{z_i} exp(-kappa (z_i - z')) S(z') dz' on a uniform grid.

    Exact for piecewise-linear S (exponential time differencing), so it stays accurate at large kappa.
    """
    dz = 1.0 / (n - 1)
    x = kappa * dz
    E = np.exp(-x)
    if abs(x) < 1e-6:
        a = b = dz / 2
    else:
        b = dz * (x - 1 + E) / x**2      # weight of the right node
        a = dz * (1 - E) / x - b         # weight of the left node
    V = np.zeros((n, n), dtype=complex)
    for i in range(1, n):
        V[i] = E * V[i - 1]
        V[i, i - 1] += a
        V[i, i] += b
    return V


def retrieval_vs_h(d, dt, S0, hmax, nh=200, n=600):
    """Adiabatic retrieval at constant control, detuning dt = Delta / gamma.

    Eliminating P in Gorshkov Eqs (5)-(7) gives, with c = 1 + i dt,
      dS/dh = -S / c + (d / c^2) int_0^z exp(-(d/c)(z - z')) S(z') dz'
      E(1)  = -(sqrt(d) Omega / c) int_0^1 exp(-(d/c)(1 - z')) S(z') dz'
    Returns h and the cumulative retrieved energy.
    """
    c = 1 + 1j * dt
    V = volterra(d / c, n)
    A = -np.eye(n, dtype=complex) / c + (d / c**2) * V
    z = np.linspace(0, 1, n)
    w = np.full(n, 1.0 / (n - 1))
    w[[0, -1]] /= 2
    S = S0(z).astype(complex)
    S /= np.sqrt(np.sum(w * abs(S)**2))
    tail = V[-1]
    dh = hmax / nh
    step = expm(A * dh)
    h, cum, out = 0.0, 0.0, []
    for _ in range(nh):
        f0 = abs(tail @ S)**2
        S = step @ S
        f1 = abs(tail @ S)**2
        cum += (d / abs(c)**2) * 0.5 * (f0 + f1) * dh
        h += dh
        out.append((h, cum))
    return np.array(out)


def erasure_table():
    say('\n== 7. Erasure: a spin wave that is not phase matched only sees the control as a pump')
    say('Spectator survival per access = exp(-2 h / (1 + dt^2)), h for 95% of the attainable retrieval.')
    for dd in (20, 150):
        say('check at dt=0, flat S, d=%d: eta(h large) = %.4f, Eq.(30) %.4f' % (
            dd, retrieval_vs_h(dd, 0, np.ones_like, 4 * dd, 200, 800)[-1, 1], 1 - (i0e(dd) + i1e(dd))))
    say('%6s %8s %10s %10s %10s %12s %12s' % ('OD', 'Delta/Gam', 'dt', 'eta_read', 'h_95', 'survive', 'accesses/2'))
    for od in (30, 100, 300):
        d = od / 2
        for D in (0, 5, 20, 50, 100, 300):
            dt = 2.0 * D
            hmax = max(3 * d, 12 * (1 + dt**2) / d)
            r = retrieval_vs_h(d, dt, np.ones_like, hmax, nh=200, n=800)
            top = r[-1, 1]
            h95 = r[np.argmax(r[:, 1] >= 0.95 * top), 0]
            surv = np.exp(-2 * h95 / (1 + dt**2))
            # number of accesses before spectators drop to 1/2
            acc = np.log(2) / max(1 - surv, 1e-12) if surv < 1 else np.inf
            say('%6d %8d %10.0f %10.3f %10.1f %12.2e %12.1f' % (od, D, dt, top, h95, surv, acc))


def control_photons(r_sig, r_ctl):
    say('\n== 8. Control photons per access: N_c = 2 h A / sigma_c = (h/d) N (r_s / r_c)')
    sig0 = 3 * LAM_S**2 / (2 * np.pi)
    e_ph = H * C / LAM_S
    for a, L, f, hd in ((15, 10, .5, 3), (15, 10, 1, 3), (3.5, 10, .5, 3), (3.5, 1, .5, 3)):
        N = f * L * 1e-3 / (LAM_L / 2)
        A = a_eff(a * 1e-6)
        d = N * sig0 * r_sig / (2 * A)
        h = hd * d
        Nc = 2 * h * A / (sig0 * r_ctl)
        say('a=%-4.1f L=%-3d f=%.1f  d=%-7.1f N=%-7.0f h=%-7.1f N_c=%.2e  energy %.1f fJ; '
            '100 ns pulse: %.0f nW' % (a, L, f, d, N, h, Nc, Nc * e_ph * 1e15, Nc * e_ph / 100e-9 * 1e9))


def noise_table(r_sig):
    say('\n== 9. Noise per read [order-of-magnitude]')
    say('Residual s population p: incoherent photons d p. FWM spin-wave probability: '
        '(h/d) d^2 (Gamma/omega_hf)^2 / 2 r_ratio')
    x = (GAMMA / (2 * np.pi * HFS))**2
    say('(Gamma/omega_hf)^2 = %.2e' % x)
    say('%6s %6s %14s %14s %14s' % ('OD', 'd', 'p=1e-3', 'p=1e-4', 'FWM (h/d=3)'))
    for od in (30, 100, 300, 1000, 3000):
        d = od / 2
        say('%6d %6.0f %14.2e %14.2e %14.2e' % (od, d, d * 1e-3, d * 1e-4, 3 * d**2 * x / 2))


# ---------------------------------------------------------------- 5. trap and storage time
def trap(u, nuz0=300e3, a=15e-6):
    """Parameters at depth fraction u of the depth that gives nuz0."""
    kl = 2 * np.pi / LAM_L
    w = 0.67 * a
    nuz = nuz0 * np.sqrt(u)
    nur = nuz * np.sqrt(2) / (kl * w)
    U = M * (2 * np.pi * nuz)**2 / (2 * kl**2)
    return dict(u=u, nuz=nuz, nur=nur, U=U, EL=U / (HBAR * WREC_L))


def lattice_shift_scatter(U):
    """Differential shift (Hz), scattering (1/s); two-line model at 1064 nm."""
    w = 2 * np.pi * C / LAM_L
    w1 = 2 * np.pi * C / 794.98e-9
    w2 = 2 * np.pi * C / 780.24e-9
    g1, g2 = 2 * np.pi * 5.746e6, 2 * np.pi * 6.065e6
    d1, d2 = w1 - w, w2 - w          # positive: red detuned
    inv = 1 / 3 / d1 + 2 / 3 / d2
    delta_eff = 1 / inv
    ddiff = U / H * (2 * np.pi * HFS) / delta_eff
    gam = (g1 * (1 / 3) / d1 + g2 * (2 / 3) / d2) * delta_eff  # weighted line width
    rs = gam / delta_eff * U / HBAR
    return ddiff, rs


def tunnelling(V):
    """Harper: J / E_r for depth V in E_r; returns J in Hz."""
    return (4 / np.sqrt(np.pi)) * V**0.75 * np.exp(-2 * np.sqrt(V)) * WREC_L / (2 * np.pi)


def axis_coherence(T, eps, up, nmax=25):
    """Coherent amplitude on one motional axis starting in n = 0.

    Recoil ladder: n -> n+1 at up (n+1), n -> n-1 at up n. Rayleigh events keep the hyperfine
    coherence, so the coherence amplitude follows the ladder and rotates at 2 pi eps n.
    """
    n = np.arange(nmax)
    D = np.zeros((nmax, nmax), dtype=complex)
    for i in range(nmax):
        D[i, i] = -up * (2 * i + 1) + 2j * np.pi * eps * i
        if i + 1 < nmax:
            D[i + 1, i] = up * (i + 1)
            D[i, i + 1] = up * (i + 1)
    v = np.zeros(nmax, dtype=complex)
    v[0] = 1
    vals, vecs = np.linalg.eig(D)
    coef = np.linalg.solve(vecs, v)
    return np.array([abs(np.sum(vecs @ (coef * np.exp(vals * t)))) for t in np.atleast_1d(T)])


def storage_terms(T, u, comp=1.0, dI=0.01, dB_hz=0.0, tau_bg=2.0, echo=False, a=15e-6):
    """Amplitude and survival factors. Returns a dict of multiplicative terms in eta(T)."""
    t = trap(u, a=a)
    ddiff, rs = lattice_shift_scatter(t['U'])
    ddiff_res = comp * ddiff
    Ufreq = t['U'] / H
    eps_r = ddiff_res * t['nur'] / (2 * Ufreq)
    eps_z = ddiff_res * t['nuz'] / (2 * Ufreq)
    wr = WREC_L / (2 * np.pi)
    up_r = rs * ALPHA * wr / t['nur']
    up_z = rs * (1 + (1 - 2 * ALPHA)) * wr / t['nuz']
    T = np.atleast_1d(T)
    c_r = axis_coherence(T, eps_r, up_r)
    c_z = axis_coherence(T, eps_z, up_z)
    # static longitudinal intensity spread: Gaussian rms dI, removed by an echo
    sigma = ddiff_res * dI
    inh = np.ones_like(T) if echo else np.exp(-(2 * np.pi * sigma * T)**2)
    mag = np.ones_like(T) if echo else np.exp(-(2 * np.pi * dB_hz * T)**2)
    return dict(heat=(c_r**4) * (c_z**2), inh=inh, mag=mag, bg=np.exp(-T / tau_bg),
                eps_r=eps_r, eps_z=eps_z, up_r=up_r, up_z=up_z, ddiff=ddiff_res, rs=rs, sigma=sigma)


def depth_table():
    say('\n== 10. Lattice at depth fraction u, case A geometry (a = 15 um, nu_z = 300 kHz at u = 1)')
    say('%8s %7s %8s %7s %7s %8s %7s %8s %8s %8s %9s %9s' % (
        'u', 'U(uK)', 'nu_z(kHz)', 'nu_r', 'E_L', 'J(Hz)', 'scat', 'up_r/s', 'up_z/s', 'ddiff', 'eps_r(Hz)', 'eps_z(Hz)'))
    for u in (1, .3, .1, .03, .01, .005, .003, .001):
        t = trap(u)
        s = storage_terms([1.0], u)
        say('%8.3f %7.1f %8.1f %7.0f %7.0f %8.1e %7.2f %8.3f %8.4f %8.1f %9.3f %9.3f' % (
            u, t['U'] / KB * 1e6, t['nuz'] / 1e3, t['nur'], t['EL'], tunnelling(t['EL']),
            s['rs'], s['up_r'], s['up_z'], s['ddiff'], s['eps_r'], s['eps_z']))
    t = trap(0.01)
    say('u=0.01: horizontal fibre sag g/omega_r^2 = %.2f um; m g a = %.2f uK vs depth %.2f uK; vertical m g d / h = %.0f Hz' % (
        G_EARTH / (2 * np.pi * t['nur'])**2 * 1e6, M * G_EARTH * 15e-6 / KB * 1e6, t['U'] / KB * 1e6, M * G_EARTH * LAM_L / 2 / H))


def thermal_vs_ground():
    say('\n== 11. Thermal against ground state: can a 10 uK sample follow the lattice down?')
    say('Adiabatic lowering keeps n, so T scales as sqrt(u) and k_B T / U as 1/sqrt(u).')
    t1 = trap(1)
    for T0 in (10e-6, 50e-6):
        for u in (1, .1, .01):
            t = trap(u)
            Tu = T0 * np.sqrt(u)
            ddiff, _ = lattice_shift_scatter(t['U'])
            say('T0=%2.0f uK u=%-5.2f T=%6.2f uK  k_BT/U=%.3f  thermal radial shift spread %.2f Hz' % (
                T0 * 1e6, u, Tu * 1e6, KB * Tu / t['U'], ddiff * KB * Tu / t['U']))
    x = np.sqrt(KB * 10e-6 / (M * (2 * np.pi * t1['nur'])**2))
    say('10 uK, u=1: x_rms = %.2f um, 1/e^2 field radius w = %.1f um, mean intensity loss 2 x^2 / w^2 * 2 = %.2f %%' % (
        x * 1e6, 0.67 * 15, 4 * x**2 / (0.67 * 15e-6)**2 * 100))


def intensity_beat():
    say('\n== 12. Lattice intensity along the column from an admixture of LP02 at 1064 nm')
    g = (jv(1, U01) / jv(1, U02))**2
    a = 15e-6
    k = 2 * np.pi / LAM_L
    dbeta = (U02**2 - U01**2) / (2 * k * a**2)
    say('on-axis field ratio squared g = %.2f, beat period %.2f mm at a = 15 um' % (g, 2 * np.pi / dbeta * 1e3))
    ddiff, _ = lattice_shift_scatter(trap(1)['U'])
    for p in (1e-2, 1e-3, 1e-4, 1e-5):
        pk = 2 * np.sqrt(p * g)
        say('power admixture p=%.0e: intensity +-%.1f %%, shift +-%.1f Hz at u=1, +-%.3f Hz at u=0.01' % (
            p, pk * 100, pk * ddiff, pk * ddiff * 0.01))


def vacuum_table():
    say('\n== 13. Pressure inside the core: molecular flow with wall outgassing q (open at both ends)')
    say('p = 3 q L^2 / (8 a^2 v), L = fibre length; q and species are assumptions (H2, q = 1e-8 Pa m/s)')
    q, v = 1e-8, 1760.
    tau_ref = 2.0 * 1e-8   # s torr, Xin 2019: about 2 s at <= 1e-8 torr; treat as an order of magnitude
    say('%6s %7s %10s %10s %10s' % ('a(um)', 'L(mm)', 'p (torr)', 'tau_bg(s)', 'OD(r=1/2,f=.5)'))
    sig = 3 * LAM_S**2 / (2 * np.pi) * 0.5
    for a, L, col in ((15, 40, 10), (15, 20, 10), (15, 10, 10), (7.5, 10, 10), (3.5, 10, 10),
                      (3.5, 3, 3), (3.5, 1, 1), (5, 2, 2)):
        p = 3 * q * (L * 1e-3)**2 / (8 * (a * 1e-6)**2 * v)
        torr = p / 133.32
        N = 0.5 * col * 1e-3 / (LAM_L / 2)
        say('%6.1f %7.0f %10.1e %10.2f %10.1f' % (a, L, torr, tau_ref / torr, N * sig / a_eff(a * 1e-6)))


def fibre_motion():
    say('\n== 14. Fibre bending modes shake the trap [hazard estimate]')
    say('ndot = m omega^3 S_x(f) / (4 hbar), S_x one-sided in m^2/Hz (Savard 1997, converted)')
    for nur in (715, 7150, 30600):
        w = 2 * np.pi * nur
        co = M * w**3 / (4 * HBAR)
        say('nu_r = %6d Hz: 0.1 quanta/s needs sqrt(S_x) < %.2f pm/rtHz' % (nur, np.sqrt(0.1 / co) * 1e12))
    E, rho = 72e9, 2200.
    for od, L in ((125e-6, 0.01), (125e-6, 0.03), (250e-6, 0.01), (250e-6, 0.03)):
        ri = od / 2
        f1 = (22.37 / (2 * np.pi * L**2)) * (ri / 2) * np.sqrt(E / rho)
        m = rho * np.pi * ri**2 * L
        x2 = KB * 300 / (m * (2 * np.pi * f1)**2)
        Q = 1e4
        Spk = (2 / np.pi) * x2 * Q / f1
        Soff = 4 * KB * 300 / (m * Q * (2 * np.pi * f1)**3) * 2 * np.pi
        co = M * (2 * np.pi * f1)**3 / (4 * HBAR)
        say('OD %3.0f um clamped span %2.0f mm: f1 = %6.0f Hz, x_rms(300 K) = %.0f pm, width f1/Q (Q=1e4); '
            'if nu_r = f1: %.2g quanta/s; off resonance: %.1e /s' % (
                od * 1e6, L * 1e3, f1, np.sqrt(x2) * 1e12, co * Spk, co * Soff))


def lifetime_table():
    say('\n== 15. eta(T): lifetime and term at T for scenarios (eta_0 = 1)')
    say('heat = motional dephasing from lattice scattering (radial axes^4, axial^2); inh = rms intensity 1 % along the column')
    scen = [('deep, no compensation', dict(u=1.0)),
            ('deep, compensation 1 %', dict(u=1.0, comp=.01)),
            ('deep, echo', dict(u=1.0, echo=True)),
            ('u=.1, no compensation', dict(u=.1)),
            ('u=.01, no compensation', dict(u=.01)),
            ('u=.01, 0.1 % intensity', dict(u=.01, dI=.001)),
            ('u=.01, echo', dict(u=.01, echo=True)),
            ('u=.01, echo, tau_bg=10 s', dict(u=.01, echo=True, tau_bg=10.)),
            ('u=.01, echo, tau_bg=100 s', dict(u=.01, echo=True, tau_bg=100.))]
    T = np.geomspace(1e-3, 300, 500)
    curves = {}
    say('%-28s %9s %9s %9s %9s %9s' % ('scenario', 'T_1/e(s)', 'heat', 'inh', 'bg', 'at T_1/e'))
    for name, kw in scen:
        s = storage_terms(T, **kw)
        eta = s['heat'] * s['inh'] * s['mag'] * s['bg']
        i = np.argmax(eta < np.exp(-1))
        say('%-28s %9.3f %9.3f %9.3f %9.3f' % (name, T[i], s['heat'][i], s['inh'][i], s['bg'][i]))
        curves[name] = eta
    return T, curves


def echo_budget():
    say('\n== 16. Ground-state purity ceiling, P0 per axis, without echo at long T')
    for P in (.9, .95, .97, .99, .999):
        say('P0 = %.3f per radial axis: long-T ceiling P0^4 = %.3f' % (P, P**4))


def make_figure(tab, T, curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    oi = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000']
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    od = np.geomspace(5, 3000, 60)
    eta = np.array([retrieval(x / 2)[0]**2 for x in od])
    ax[0].semilogx(od, 1 - eta, color=oi[0], lw=1.8, label='Storage + backward retrieval')
    fw = [total_efficiency(x / 2, forward=True) for x in od[od < 1200]]
    ax[0].semilogx(od[od < 1200], 1 - np.array(fw), color=oi[1], lw=1.8, ls='--', label='Storage + forward retrieval')
    ax[0].semilogx(od, 11.6 / od, color='0.5', lw=1, ls=':', label='11.6 / OD')
    ax[0].set_yscale('log')
    ax[0].set_ylim(5e-3, 1.2)
    ax[0].set_xlabel('Optical depth OD = 2d')
    ax[0].set_ylabel('Error, 1 − η')
    ax[0].legend(frameon=False, loc='lower left')
    ax[0].set_title('a  Optimal efficiency', loc='left', fontweight='bold', fontsize=9)
    show = ['deep, no compensation', 'deep, echo', 'u=.01, no compensation', 'u=.01, echo',
            'u=.01, echo, tau_bg=10 s', 'u=.01, echo, tau_bg=100 s']
    labels = ['Deep, no compensation', 'Deep, echo', 'Depth ×0.01, 1 % intensity spread', 'Depth ×0.01, echo, τ_bg = 2 s',
              'Depth ×0.01, echo, τ_bg = 10 s', 'Depth ×0.01, echo, τ_bg = 100 s']
    for i, (n, l) in enumerate(zip(show, labels)):
        ax[1].semilogx(T, curves[n], color=oi[i], lw=1.6, label=l)
    ax[1].set_xlabel('Storage time T (s)')
    ax[1].set_ylabel('η(T) / η(0)')
    ax[1].set_ylim(0, 1.03)
    ax[1].set_xlim(1e-3, 300)
    ax[1].legend(frameon=False, loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=7)
    ax[1].set_title('b  Storage model, a = 15 µm', loc='left', fontweight='bold', fontsize=9)
    fig.savefig(os.path.join(HERE, 'memory_budget.png'), dpi=200)


if __name__ == '__main__':
    OUT = open(os.path.join(HERE, 'budget_output.txt'), 'w')
    rows = lambda_table()
    # reference scheme: clock pair, g = |1,0>, s = |2,0>, via F'=2, m=+1 (sigma+ for both beams)
    ref = [r for r in rows if r[0] == (1, 0) and r[1] == (2, 0) and r[2] == 2 and r[3] == 1][0]
    r_s, r_c = ref[4], ref[5]
    validate_kernel()
    tab = efficiency_table()
    geometry_table(r_s)
    mismatch_table()
    capacity_table()
    erasure_table()
    control_photons(r_s, r_c)
    noise_table(r_s)
    depth_table()
    thermal_vs_ground()
    intensity_beat()
    vacuum_table()
    fibre_motion()
    T, curves = lifetime_table()
    echo_budget()
    make_figure(tab, T, curves)
    OUT.close()
