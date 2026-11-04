"""87Rb D-line data and hyperfine line strengths.

Data from Steck, "Rubidium 87 D Line Data", rev. 2.3.4 (2025), Tables 3-5, 7
and the hyperfine level diagrams. Ground state 5S1/2 (J = 1/2, F = 1, 2),
nuclear spin I = 3/2; excited state 5P1/2 (D1, F' = 1, 2) or 5P3/2 (D2,
F' = 0..3). Line strengths, branching ratios and Zeeman couplings follow from
Wigner 3j/6j symbols (Steck Eqs. 36-42), not from typed tables.
"""
from functools import lru_cache

import numpy as np
from sympy import Rational, sqrt
from sympy.physics.wigner import clebsch_gordan, wigner_6j

I_NUC = Rational(3, 2)
J_G = Rational(1, 2)
GROUND_F = (1, 2)
HFS_GROUND = 6.834682610904e9          # Hz, F = 2 above F = 1

# lam: vacuum wavelength (m); Gamma: decay rate (s^-1); Je: excited J;
# isat: far-detuned saturation intensity for pi light (W/m^2), i.e. for the
# effective dipole <J||er||J'>/sqrt(3); hfs: excited hyperfine energies (Hz)
# relative to the lowest F'.
# D1 keeps lam = 794.979 nm and Gamma = 2 pi x 5.746 MHz, the values used for
# all paper results; rev. 2.3.4 gives 794.978851 nm and 2 pi x 5.7500(56) MHz
# (0.7 sigma, 0.07 % in Gamma).
LINES = {
    'D1': dict(lam=794.979e-9, Gamma=2 * np.pi * 5.746e6, Je=Rational(1, 2),
               isat=44.876, hfs={1: 0.0, 2: 814.5e6}),
    'D2': dict(lam=780.241209686e-9, Gamma=2 * np.pi * 6.0666e6, Je=Rational(3, 2),
               isat=25.0399, hfs={0: 0.0, 1: 72.2180e6, 2: 72.2180e6 + 156.9470e6,
                                  3: 72.2180e6 + 156.9470e6 + 266.6500e6}),
}


def excited_levels(line):
    """Hyperfine quantum numbers F' of the excited state."""
    return tuple(LINES[line]['hfs'])


@lru_cache(maxsize=None)
def _six_j_sq(line, F, Fe):
    Je = LINES[line]['Je']
    return wigner_6j(J_G, Je, 1, Fe, F, I_NUC)**2


def strength(line, F, Fe):
    """Relative hyperfine transition strength S_FF' = (2F'+1)(2J+1){J J' 1; F' F I}^2.

    Sum over F' is 1. For pi light and a uniform m_F population the effective
    dipole of F -> F' is S_FF' times the far-detuned one, so the saturation
    intensity of the resolved transition is isat / S_FF'.
    """
    return float((2 * Fe + 1) * (2 * J_G + 1) * _six_j_sq(line, F, Fe))


def branching(line, Fe, F):
    """Probability that F' decays into the ground level F: (2F+1)(2J'+1){J J' 1; F' F I}^2."""
    Je = LINES[line]['Je']
    return float((2 * F + 1) * (2 * Je + 1) * _six_j_sq(line, F, Fe))


@lru_cache(maxsize=None)
def zeeman_coupling(line, F, mF, Fe, mFe, q):
    """<F mF| e r_q |F' mF'> in units of <J||er||J'> (Steck Eqs. 36-37).

    Nonzero only for mF = mFe + q; q = 0 is pi, q = +-1 sigma.
    """
    Je = LINES[line]['Je']
    if mF != mFe + q:
        return 0.0
    red = ((-1)**(Fe + J_G + 1 + I_NUC) * sqrt((2 * Fe + 1) * (2 * J_G + 1))
           * wigner_6j(J_G, Je, 1, Fe, F, I_NUC))
    cg = clebsch_gordan(Fe, 1, F, mFe, q, mF)
    return float(red * cg)


def lambda_config(line, Fe, gp):
    """Closed Lambda via the excited level F', probe on the ground level F = gp.

    Returns the branching ratios and line strengths that the three-level model
    needs: bp (decay into g_p), Sp, Sc (S_FF' of the probe and coupling legs)
    and the detunings (Hz) of the other excited levels from F'.
    """
    if gp not in GROUND_F or Fe not in LINES[line]['hfs']:
        raise ValueError(f'no level F={gp} or F\'={Fe} on {line}')
    gc = 3 - gp
    Sp, Sc = strength(line, gp, Fe), strength(line, gc, Fe)
    if Sp == 0 or Sc == 0:
        raise ValueError(f'F\'={Fe} on {line} does not couple to both F = 1 and 2')
    hfs = LINES[line]['hfs']
    others = {f: hfs[f] - hfs[Fe] for f in hfs if f != Fe}
    return dict(line=line, Fe=Fe, gp=gp, gc=gc, bp=branching(line, Fe, gp),
                bc=branching(line, Fe, gc), Sp=Sp, Sc=Sc, others=others)


def neighbour_shift(line, Fe, gc, Delta):
    """Light shift of g_c by the other excited levels relative to that of F'.

    Far-off-resonant estimate sum_F'' (S_gc,F''/S_gc,F') Delta/Delta'', with
    Delta (rad/s) the blue detuning from F' and Delta'' = Delta - (E_F'' - E_F')
    the detuning from F''. A preview of the multi-level effects of stage 2.
    """
    hfs = LINES[line]['hfs']
    S0 = strength(line, gc, Fe)
    tot = 0.0
    for f, e in hfs.items():
        if f == Fe or strength(line, gc, f) == 0:
            continue
        tot += strength(line, gc, f) / S0 * Delta / (Delta - 2 * np.pi * (e - hfs[Fe]))
    return tot
