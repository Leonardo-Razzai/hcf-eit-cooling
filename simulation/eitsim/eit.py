"""EIT Lambda system: laser settings, emission pattern and Fano profile.

All quantities are dimensionless, in units of the angular radial trap
frequency (nu_r = 1). The parameter dict returned by params() is the common
input of the master-equation and rate models.
"""
import numpy as np

from .fibre import gfac

GP, GC, EX = 0, 1, 2             # atomic level indices: |g_p>, |g_c>, |e>


def emission_nodes(K=10, seed=1):
    """Quantile quadrature for u = projection of the emission direction on x.

    Dipole along the fibre axis z (pattern ~ sin^2 theta), <u^2> = 2/5.
    Returns K equal-weight nodes, symmetrised.
    """
    rng = np.random.default_rng(seed)
    n = 2_000_000
    c = rng.uniform(-1, 1, 3 * n)
    keep = rng.uniform(0, 1, 3 * n) < 1 - c**2
    c = c[keep][:n]
    phi = rng.uniform(0, 2 * np.pi, c.size)
    u = np.sort(np.sqrt(1 - c**2) * np.cos(phi))
    nodes = np.array([b.mean() for b in np.array_split(u, K)])
    nodes = 0.5 * (nodes - nodes[::-1])
    return nodes, np.full(K, 1.0 / K)


NODES, WEIGHTS = emission_nodes()
ALPHA = float(np.sum(WEIGHTS * NODES**2))      # <u_x^2> = 0.397


def params(fc, s=0.5, Delta_G=3.0, m=1, eps=0.0, xs=0.0, d2=0.0,
           nu_rel=1.0, bp=0.5, N=30):
    """Dimensionless parameter set (nu_r = 1) for fibre case fc.

    s       saturation parameter 2 Om_perp^2 / (Gamma Gamma_n)
    Delta_G single-photon detuning in units of Gamma
    m       radial order of the LP1m probe
    eps     LP01 amplitude admixture in the probe (relative to peak-equivalent)
    xs      transverse offset of the trap centre from the node, units of x0
    d2      two-photon detuning error, units of nu_r
    nu_rel  actual nu_r / design nu_r (lasers set for the design value)
    bp      branching ratio of |e> into |g_p>
    N       Fock-space truncation
    """
    # design values (units of the design trap frequency)
    G = fc['Gamma'] / fc['om']
    D = Delta_G * G
    Wc = np.sqrt(4 * (D + 1))                       # narrow peak at delta = nu
    Gn = G * 0.5 * (1 - D / np.sqrt(D**2 + Wc**2))  # width of the narrow state
    etap = gfac(m) * fc['x0a']
    Wp = np.sqrt(s * G * Gn / 2) / etap
    # same lasers, actual trap frequency nu_rel * nu: rescale to actual units
    r = 1.0 / nu_rel
    return dict(G=G * r, D=D * r, Wc=Wc * r, Gn=Gn * r, Wp=Wp * r,
                etap=etap * np.sqrt(r), x0a=fc['x0a'] * np.sqrt(r),
                eta_sp=fc['eta_sp'] * np.sqrt(r), m=m, eps=eps, xs=xs, d2=d2,
                bp=bp, N=N, s=s, nu_rel=nu_rel)


def fano(p, delta):
    """Weak-probe absorption rate per unit Om_p^2 at two-photon detuning delta."""
    delta = np.asarray(delta, float)
    Dd = (p['D'] + delta) * delta - p['Wc']**2 / 4
    return p['G'] * delta**2 / (p['G']**2 * delta**2 + 4 * Dd**2)
