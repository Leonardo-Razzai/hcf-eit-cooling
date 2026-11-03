"""Rate models: the simple sideband model and the generalised rate model (GRM).

simple():   weak-probe sideband rate equation (Morigi et al. 2000) with
            kappa photons per cycle and saturation of the red sideband only.
unified():  saturation-optimised rate with an LP01 admixture eps.
grm_matrix(): Markov chain on Fock populations of |g_p, n> with exact
            Franck-Condon factors of the Bessel profile and of the recoil
            kernel; its saturation constant is calibrated against the master
            equation (studies/calibration.py).
"""
import numpy as np

from .eit import NODES, WEIGHTS, ALPHA, fano
from .operators import oscillator, profile_ops

# GRM saturation constant fitted to the master-equation phonon drift
# (Fock states n0 = 1..15, s = 0.2 and 1, both fibre cases; 17 % rms).
CSAT = 2.25


def simple(p, photons=1.0, sat=True):
    """Simple model: rate W, steady state n and sideband rates (units of nu_r)."""
    Wperp2 = (p['Wp'] * p['etap'])**2
    Am = Wperp2 * fano(p, 1.0 + p['d2'])
    Ap = Wperp2 * fano(p, -1.0 + p['d2'])
    if sat:                       # only the resonant red sideband saturates
        Am = Am / (1 + p['s'])
    h = photons * ALPHA * p['eta_sp']**2
    W = Am * (1 - h) - Ap * (1 + h)
    return dict(W=W, n=Ap * (1 + h) / W if W > 0 else np.inf, Am=Am, Ap=Ap)


def unified(p, photons=2.0):
    """Saturation-optimised simple-model rate with LP01 admixture p['eps']."""
    r = np.inf if p['eps'] == 0 else p['etap'] / p['eps']
    s = min(1.0, 8 * (p['D'] / (3 * p['G']))**2 * r**2)
    return (1 - photons * ALPHA * p['eta_sp']**2) * p['Gn'] / 2 * s / (1 + s)


def grm_matrix(p, N=300, sat=True, csat=CSAT):
    """GRM rate matrix Q (units of nu_r): dp/dt = Q p for populations of |g_p, n>."""
    b, w, V = oscillator(N)
    fp, _ = profile_ops(dict(p, N=N), w, V)
    fpd = fp.toarray()
    F2 = np.abs(fpd)**2                                # |<m|f_p|n>|^2
    n = np.arange(N)
    dn = n[None, :] - n[:, None]                       # phonons removed n -> m
    kabs = p['Wp']**2 * F2 * fano(p, dn + p['d2'])
    if sat:
        s_n = 2 * p['Wp']**2 * np.abs(np.diag(fpd, 1))**2 / (p['G'] * p['Gn'])
        s_n = np.concatenate([[0.0], s_n])             # sideband n -> n-1
        kabs = kabs / (1 + csat * s_n[None, :])
    P = np.zeros((N, N))                               # one-photon recoil kernel
    for u, wk in zip(NODES, WEIGHTS):
        Dk = (V * np.exp(1j * p['eta_sp'] * u * w)) @ V.conj().T
        P += wk * np.abs(Dk)**2
    P /= P.sum(axis=0, keepdims=True)
    T = p['bp'] * P @ np.linalg.inv(np.eye(N) - (1 - p['bp']) * P)
    Q = T @ kabs
    Q -= np.diag(kabs.sum(axis=0))
    return Q


def grm_steady_state(Q):
    w, v = np.linalg.eig(Q)
    ss = np.real(v[:, np.argmin(np.abs(w))])
    return ss / ss.sum()


def thermal(nbar, N):
    """Thermal Fock distribution with mean nbar, truncated at N."""
    if nbar <= 0:
        v = np.zeros(N)
        v[0] = 1
        return v
    q = nbar / (nbar + 1)
    v = (1 - q) * q**np.arange(N)
    return v / v.sum()
