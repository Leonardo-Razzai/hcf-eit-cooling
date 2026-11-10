"""Physical constants, hollow-core fibre modes and lattice trap parameters.

Capillary (Marcatili-Schmeltzer) model: LP_lm modes J_l(u_lm r/a) cos(l phi)
with J_l(u_lm) = 0. The LP01 lattice is approximated by a Gaussian with waist
w = 0.67 a for the radial trap frequency.
"""
import numpy as np
import scipy.constants as sc
from scipy.special import jv, jn_zeros

from .atoms import LINES

# physical constants (SI, CODATA via scipy.constants)
HBAR = sc.hbar
KB = sc.k

# 87Rb mass and D1 line (Steck, Rubidium 87 D line data; all lines in atoms.py)
# Steck's tabulated kg value (86.909180531 u converted with an older CODATA u; 1.7e-7 above the
# CODATA 2022 conversion). Kept so that stored results reproduce; the difference is physically irrelevant.
M_RB = 1.443160648e-25
LAM_D1 = LINES['D1']['lam']
GAMMA = LINES['D1']['Gamma']

# mode constants
U01 = jn_zeros(0, 1)[0]          # LP01: J0(u) = 0  -> 2.405
U1 = jn_zeros(1, 10)             # LP1m: J1(u) = 0  -> 3.832, 7.016, ...


def amp_ratio(m=1):
    """Peak-amplitude ratio A_1m / A_01 at equal guided power."""
    return np.sqrt(2) * abs(jv(1, U01)) / abs(jv(2, U1[m - 1]))


def gfac(m=1):
    """On-axis gradient factor G_1m: field = A_01 G_1m x / a near the axis."""
    return U1[m - 1] / 2 * amp_ratio(m)


def fibre_case(a, nuz=300e3, lam_lat=1064e-9, line='D1'):
    """Radial trap parameters of a 1D LP01 lattice in a core of radius a (m).

    nuz is the axial trap frequency (Hz), line the cooling line ('D1', 'D2').
    Returns a dict with the radial frequency nur (Hz), om = 2 pi nur,
    ground-state size x0, x0/a, the spontaneous-emission Lamb-Dicke parameter
    eta_sp = k_line x0, the lattice depth U0 (J) that gives nuz, and the line
    data (line, Gamma, lam).
    """
    w = 0.67 * a
    kl = 2 * np.pi / lam_lat
    nur = nuz * np.sqrt(2) / (kl * w)
    om = 2 * np.pi * nur
    x0 = np.sqrt(HBAR / (2 * M_RB * om))
    U0 = M_RB * (2 * np.pi * nuz)**2 / (2 * kl**2)
    lam = LINES[line]['lam']
    eta_sp = 2 * np.pi / lam * x0
    return dict(a=a, nur=nur, om=om, x0=x0, eta_sp=eta_sp, x0a=x0 / a, U0=U0,
                line=line, Gamma=LINES[line]['Gamma'], lam=lam, nu_wrec=1 / eta_sp**2)


def reduced_case(nu_wrec, Gamma_nu=800.0, x0a=0.006):
    """Atom- and fibre-independent case in internal units (om = nu_r = 1).

    nu_wrec = nu_r / omega_rec = 1 / eta_sp^2 is the only physical parameter of
    the cold-atom problem (closed Lambda, node geometry, Gamma >> nu_r): with
    Gamma/nu_r and x0/a varied by factors 16 and 4 the master-equation rate and
    n_ss change by < 0.3 %. Gamma_nu and x0a are therefore nuisance parameters;
    x0a matters only for hot states (Bessel non-linearity) and the probe power.
    Rates computed with this case are in units of nu_r.
    """
    return dict(om=1.0, Gamma=Gamma_nu, eta_sp=nu_wrec**-0.5, x0a=x0a,
                nu_wrec=nu_wrec, line=None)


def lp01_power(fc, Om, S=1.0):
    """Guided power (W) whose LP01-equivalent peak Rabi frequency is Om (rad/s)
    on a hyperfine transition of strength S (atoms.strength).

    I_peak = 2 I_sat (Om / Gamma)^2 / S, and the LP01 mode of peak intensity
    I_peak carries I_peak pi a^2 J1(u01)^2.
    """
    ln = LINES[fc['line']]
    I = 2 * ln['isat'] * (Om / ln['Gamma'])**2 / S
    return I * np.pi * fc['a']**2 * jv(1, U01)**2


# the two cases of the paper
CASE_A = dict(a=15e-6)           # large-core antiresonant fibre
CASE_B = dict(a=3.5e-6)          # small-core photonic-bandgap fibre


def paper_cases(line='D1'):
    """{'A': fibre_case(15 um), 'B': fibre_case(3.5 um)} on the given line."""
    return {'A': fibre_case(**CASE_A, line=line), 'B': fibre_case(**CASE_B, line=line)}
