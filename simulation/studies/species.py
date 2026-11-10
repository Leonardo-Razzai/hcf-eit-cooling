"""Other atoms, ions and tweezers on the cooling map (nu_r / omega_rec, b_p), with h = alpha omega_rec / (b_p nu_r).

  points  literature experiments and the fibre geometry for Li, Na, K, Rb, Cs (D1), 171Yb and 87Sr
          -> results/species_points.txt (name, category, nu_wrec, bp, h, note)
  sims    master-equation checks where the stored map does not apply: narrow-line Yb (Gamma / nu_r = 36,
          below the validated range) and K (nu_r / omega_rec = 0.83, below the grid)
          -> results/species_sims.dat

Conventions. omega_rec = hbar k^2 / 2m of the cooling (probe) line. For EIT schemes b_p is the branching
ratio of the excited level into the probe ground level (Clebsch-Gordan, D1 hyperfine or Zeeman Lambda).
For schemes without a dark resonance (gray molasses, Raman or narrow-line sideband cooling) the map
coordinate is 1 / kappa, kappa = photons per removed phonon (estimate), marked 'approx'.
Literature numbers are taken from the papers listed with each point (full texts read, Nov 2026).
"""
import sys

import numpy as np
import scipy.constants as sc
from sympy import Rational, S
from sympy.physics.wigner import wigner_6j

from _common import E, save, log, RESULTS

ALPHA = E.ALPHA
U = sc.atomic_mass


def wrec(mass_u, lam):
    """omega_rec / 2 pi in Hz."""
    return sc.h / (2 * mass_u * U * lam**2)


def d1_branching(I, Fe, F):
    """Branching of an excited hyperfine level F' (J' = 1/2) into ground F (J = 1/2), nuclear spin I."""
    J = Rational(1, 2)
    tot = sum((2 * f + 1) * wigner_6j(J, J, 1, Fe, f, I)**2 for f in (I - J, I + J))
    return float((2 * F + 1) * wigner_6j(J, J, 1, Fe, F, I)**2 / tot)


def best_d1_bp(I):
    """Largest and smallest b_p over the four D1 hyperfine Lambda configurations."""
    J = Rational(1, 2)
    vals = [d1_branching(I, Fe, F) for Fe in (I - J, I + J) for F in (I - J, I + J)]
    return max(vals), min(vals)


# fibre geometry: radial frequency is species independent at fixed nu_z for a given lattice wavelength
def fibre_nur(a, nuz=300e3, lam_lat=1064e-9):
    """Radial trap frequency (Hz) of the LP01 lattice in a fibre of core radius a (eitsim.fibre_case)."""
    return E.fibre_case(a, nuz=nuz, lam_lat=lam_lat)['nur']


ALKALI = {  # name: (mass u, D1 wavelength, nuclear spin)
    '6Li': (6.0151, 670.992e-9, S(1)),
    '23Na': (22.9898, 589.756e-9, Rational(3, 2)),
    '39K': (38.9637, 770.108e-9, Rational(3, 2)),
    '87Rb': (86.9092, 794.979e-9, Rational(3, 2)),
    '133Cs': (132.9055, 894.593e-9, Rational(7, 2)),
}


def points():
    rows = []

    def add(name, cat, nu, w, bp, note):
        x = nu / w
        h = ALPHA / (bp * x)
        rows.append((name, cat, x, bp, h, note))
        log('%-28s %-8s nu/wrec %8.2f  b_p %.3f  h %.3g  %s' % (name, cat, x, bp, h, note))

    # --- trapped ions (EIT), J = 1/2 -> 1/2 Zeeman Lambda: pi probe gives b_p = 1/3
    add('Ca+_Roos2000', 'ion', 1.62e6, wrec(39.9626, 396.96e-9), 1 / 3, 'radial 1.62 MHz, 397 nm; n ~ 0.1 from 16 in < 1 ms; Roos et al. PRL 85, 5547 (2000)')
    add('Ca+_Xing2025', 'ion', 3.0e6, wrec(39.9626, 396.96e-9), 1 / 3, 'standing-wave EIT, radial mode, n ~ 0.05 in 150 us; Xing et al. arXiv:2512.10900')
    add('Yb+_Feng2020', 'ion', 4.45e6, wrec(170.9363, 369.53e-9), 1 / 3, 'tripod, transverse 4.45 MHz, n = 0.04-0.08; Feng et al. arXiv:2004.05190 (b_p approx)')
    add('Ca+_Lechner2016', 'ion', 2.08e6, wrec(39.9626, 396.96e-9), 1 / 3, 'ion strings up to 18, radial 2.08 MHz, n = 0.01-0.02; Lechner et al. arXiv:1603.05568')
    add('Be+_Lin2013', 'ion', 2.97e6, wrec(9.0122, 313.1e-9), 1 / 3, 'sympathetic, 2.97 MHz, n = 0.04-0.08; Lin et al. arXiv:1211.6647 (b_p approx)')
    add('Mg+_Wu2024', 'ion', 1.5e6, wrec(24.9858, 279.6e-9), 1 / 3, 'tripod, 1.5 MHz, n = 0.26; Wu et al. arXiv:2408.13407 (b_p approx)')
    # --- tweezers and nanophotonic traps
    add('Rb87_Qiu2026_radial', 'tweezer', 45e3, wrec(86.9092, 794.979e-9), 1 / 6, 'EIT loading/imaging, Zeeman Lambda in F=2 (b_p = 1/6), radial 45 kHz, n_r ~ 3.7 while imaging; Qiu et al. arXiv:2509.12124')
    add('Rb87_Qiu2026_axial', 'tweezer', 4e3, wrec(86.9092, 794.979e-9), 1 / 6, 'same, axial 4 kHz (outside Lamb-Dicke)')
    add('Li6_Blodgett2023', 'tweezer', 700e3, wrec(6.0151, 670.992e-9), 0.5, 'Lambda-GM D1, radial 700 kHz, ~11 % radial ground state while imaging; Blodgett et al. arXiv:2305.02405 (approx)')
    add('K39_Angonga2022', 'tweezer', 180e3, wrec(38.9637, 770.108e-9), 0.5, 'gray molasses D1, radial 180 kHz; Ang\'ong\'a et al. arXiv:2108.11895 (approx)')
    add('Cs_nanofibre_Pache2026', 'tweezer', 46e3, wrec(132.9055, 852.347e-9), 0.5, 'Lambda-GM (D1 F=4) and guided-light EIT (D2) in a nanofibre trap, 46 kHz; Pache et al. arXiv:2605.13387 (approx)')
    add('Sr88_Cooper2018', 'tweezer', 211e3, wrec(87.9056, 689.449e-9), 1.0, 'narrow-line sideband, radial 211 kHz, ground fraction 0.69-0.96; Cooper et al. arXiv:1810.06537 (kappa ~ 1)')
    add('Yb171_Jenkins2022', 'tweezer', 139e3, wrec(170.9363, 555.8e-9), 0.5, 'Raman sideband, radial 139 kHz, 556 nm pumping; Jenkins et al. arXiv:2112.06732 (kappa ~ 2)')
    # --- this fibre: alkali D1, best and worst Lambda
    for a, tag in ((15e-6, 'A'), (3.5e-6, 'B')):
        nur = fibre_nur(a)
        for name, (m, lam, I) in ALKALI.items():
            hi, lo = best_d1_bp(I)
            add('%s_fibre%s' % (name, tag), 'fibre', nur, wrec(m, lam), hi, 'D1, best Lambda (b_p %.3f; worst %.3f), a = %.1f um' % (hi, lo, a * 1e6))
        # 171Yb: 1S0 F = 1/2 Zeeman Lambda via 3P1 (556 nm) or 1P1 (399 nm), magic lattice 759 nm
        nurY = fibre_nur(a, lam_lat=759.35e-9)
        add('Yb171_556_fibre%s' % tag, 'fibre', nurY, wrec(170.9363, 555.8e-9), 2 / 3, '3P1 (Gamma = 183 kHz), sigma probe b_p = 2/3 (pi probe 1/3), lattice 759 nm')
        add('Yb171_399_fibre%s' % tag, 'fibre', nurY, wrec(170.9363, 398.9e-9), 2 / 3, '1P1 (29 MHz), b_p = 2/3')
        # 87Sr: Zeeman Lambda in F = 9/2; 689 nm narrow line or 461 nm broad line, magic lattice 813 nm
        nurS = fibre_nur(a, lam_lat=813.43e-9)
        add('Sr87_689_fibre%s' % tag, 'fibre', nurS, wrec(86.9089, 689.449e-9), 0.5, '3P1, Gamma = 7.5 kHz ~ nu_r: resolved sideband, not EIT (b_p approx)')
        add('Sr87_461_fibre%s' % tag, 'fibre', nurS, wrec(86.9089, 460.862e-9), 0.5, '1P1 (30 MHz), b_p approx')
    out = RESULTS / 'species_points.txt'
    with open(out, 'w') as f:
        f.write('name category nu_wrec bp h note\n')
        for name, cat, x, bp, h, note in rows:
            f.write('%s %s %.6g %.6g %.6g "%s"\n' % (name, cat, x, bp, h, note))
    log('wrote', out)


def sims():
    """Master-equation checks: narrow-line Yb (Gamma/nu_r = 36 against the default 800) and K in the fibre."""
    rows = []
    nurY = fibre_nur(15e-6, lam_lat=759.35e-9)
    xY = nurY / wrec(170.9363, 555.8e-9)
    gY = 183e3 / nurY
    xK = fibre_nur(15e-6) / wrec(38.9637, 770.108e-9)
    for label, nw, bp, gnu in (('Yb556 b_p=2/3, Gamma/nu=%.0f' % gY, xY, 2 / 3, gY), ('Yb556 b_p=2/3, Gamma/nu=800', xY, 2 / 3, 800.0),
                               ('Yb556 b_p=1/3, Gamma/nu=%.0f' % gY, xY, 1 / 3, gY), ('Yb556 b_p=1/3, Gamma/nu=800', xY, 1 / 3, 800.0),
                               ('K39 D1 b_p=5/6', xK, 5 / 6, 800.0), ('K39 D1 b_p=1/2', xK, 0.5, 800.0)):
        for s in (0.3, 1.0):
            r = E.steady(E.reduced_case(nw, Gamma_nu=gnu), s, cut=(250, 4), Delta_G=3.0, bp=bp)
            h = ALPHA / (bp * nw)
            rows.append([nw, bp, gnu, s, h, r['P0'], r['n'], r['W0']])
            log('sims %-30s nu/wrec %.3f h %.3f s %g: P0 %.4f n %.4f W0/nu_r %.4f' % (label, nw, h, s, r['P0'], r['n'], r['W0']))
    save('species_sims.dat', ['nu_wrec', 'bp', 'Gamma_nu', 's', 'h', 'P0', 'n', 'W0'], rows)


STUDIES = {'points': points, 'sims': sims}

if __name__ == '__main__':
    for name in sys.argv[1:]:
        STUDIES[name]()
