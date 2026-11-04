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
           nu_rel=1.0, bp=0.5, N=30, levels=None, sign=1.0):
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
    levels  excited hyperfine manifold from atoms.excited_manifold (None:
            single excited level, the closed three-level Lambda); the first
            entry is the level the lasers are tuned to, and its bp replaces bp
    sign    factor on the probe amplitudes of the other levels (+1: signs of
            the reduced hyperfine elements, -1: reversed Raman interference)

    p['levels'] is a list of (offset, cp, cc, bp) per excited level, offset
    in units of nu_r. Delta, Wp, Wc and Gn refer to the first level.
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
    if levels is None:
        lv = [(0.0, 1.0, 1.0, bp)]
    else:
        lv = [(2 * np.pi * l['eps'] / fc['om'] * r, l['cp'] * (sign if k else 1.0),
               l['cc'], l['bp']) for k, l in enumerate(levels)]
        bp = lv[0][3]
    return dict(G=G * r, D=D * r, Wc=Wc * r, Gn=Gn * r, Wp=Wp * r, levels=lv,
                etap=etap * np.sqrt(r), x0a=fc['x0a'] * np.sqrt(r),
                eta_sp=fc['eta_sp'] * np.sqrt(r), m=m, eps=eps, xs=xs, d2=d2,
                bp=bp, N=N, s=s, nu_rel=nu_rel)


def fano(p, delta):
    """Weak-probe absorption rate per unit Om_p^2 at two-photon detuning delta."""
    delta = np.asarray(delta, float)
    Dd = (p['D'] + delta) * delta - p['Wc']**2 / 4
    return p['G'] * delta**2 / (p['G']**2 * delta**2 + 4 * Dd**2)



def _manifold(p):
    """Non-Hermitian single-excitation Hamiltonian on {g_c, e_1..e_K} and the
    probe source vector (per unit Om_p), in the frame of params()."""
    lv = p['levels']
    K = len(lv)
    H = np.zeros((K + 1, K + 1), complex)            # index 0: g_c, k: e_k
    src = np.zeros(K + 1, complex)
    for k, (off, cp, cc, _) in enumerate(lv, start=1):
        H[k, k] = -(p['D'] - off) - 0.5j * p['G']
        H[0, k] = H[k, 0] = 0.5 * p['Wc'] * cc
        src[k] = 0.5 * cp
    return H, src


def _amplitudes(p, d):
    H, src = _manifold(p)
    return np.linalg.solve(H - d * np.eye(len(src)), -src)[1:]


def absorption(p, delta):
    """Weak-probe absorption rate per unit Om_p^2 for the excited manifold
    p['levels'], at two-photon detuning delta (scalar or array).

    The single-excitation amplitudes v of {g_c, e_k} driven from g_p solve
    (H_eff - delta) v = -s, with
    H_eff = -sum_k (Delta - eps_k + i Gamma/2)|e_k><e_k| + (Om_c/2) sum_k cc_k (|e_k><g_c| + h.c.)
    and s = (Om_p/2) sum_k cp_k |e_k>; the rate is Gamma sum_k |v_k|^2.
    Identical to fano() for a single level, which is used then.
    """
    if len(p['levels']) == 1:
        return fano(p, delta)
    delta = np.asarray(delta, float)
    out = [p['G'] * np.sum(np.abs(_amplitudes(p, d))**2) for d in delta.ravel()]
    return np.array(out).reshape(delta.shape)


def excited_weights(p, delta=1.0):
    """Fraction of the probe excitation in each excited level at detuning delta
    (default: the red sideband); weights the branching ratio in rate models."""
    if len(p['levels']) == 1:
        return np.array([1.0])
    v = np.abs(_amplitudes(p, delta))**2
    return v / v.sum()
