"""Physical constants, hollow-core fibre modes and lattice trap parameters.

Capillary (Marcatili-Schmeltzer) model: LP_lm modes J_l(u_lm r/a) cos(l phi)
with J_l(u_lm) = 0. The LP01 lattice is approximated by a Gaussian with waist
w = 0.67 a for the radial trap frequency.
"""
import numpy as np
from scipy.special import jv, jn_zeros

# physical constants (SI)
HBAR = 1.054571817e-34
KB = 1.380649e-23

# 87Rb D1 line (Steck, Rubidium 87 D line data)
M_RB = 1.443160648e-25
LAM_D1 = 794.979e-9
GAMMA = 2 * np.pi * 5.746e6

# mode constants
U01 = jn_zeros(0, 1)[0]          # LP01: J0(u) = 0  -> 2.405
U1 = jn_zeros(1, 10)             # LP1m: J1(u) = 0  -> 3.832, 7.016, ...


def amp_ratio(m=1):
    """Peak-amplitude ratio A_1m / A_01 at equal guided power."""
    return np.sqrt(2) * abs(jv(1, U01)) / abs(jv(2, U1[m - 1]))


def gfac(m=1):
    """On-axis gradient factor G_1m: field = A_01 G_1m x / a near the axis."""
    return U1[m - 1] / 2 * amp_ratio(m)


def fibre_case(a, nuz=300e3, lam_lat=1064e-9):
    """Radial trap parameters of a 1D LP01 lattice in a core of radius a (m).

    nuz is the axial trap frequency (Hz). Returns a dict with the radial
    frequency nur (Hz), om = 2 pi nur, ground-state size x0, x0/a and the
    spontaneous-emission Lamb-Dicke parameter eta_sp = k_795 x0 and the
    lattice depth U0 (J) that gives nuz.
    """
    w = 0.67 * a
    kl = 2 * np.pi / lam_lat
    nur = nuz * np.sqrt(2) / (kl * w)
    om = 2 * np.pi * nur
    x0 = np.sqrt(HBAR / (2 * M_RB * om))
    U0 = M_RB * (2 * np.pi * nuz)**2 / (2 * kl**2)
    return dict(a=a, nur=nur, om=om, x0=x0, eta_sp=2 * np.pi / LAM_D1 * x0,
                x0a=x0 / a, U0=U0)


# the two cases of the paper
CASE_A = dict(a=15e-6)           # large-core antiresonant fibre
CASE_B = dict(a=3.5e-6)          # small-core photonic-bandgap fibre


def paper_cases():
    """{'A': fibre_case(15 um), 'B': fibre_case(3.5 um)}."""
    return {'A': fibre_case(**CASE_A), 'B': fibre_case(**CASE_B)}
