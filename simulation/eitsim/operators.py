"""Harmonic-oscillator operators and field profiles as functions of x.

Functions of the position operator are evaluated in the eigenbasis of the
truncated x = x0 (b + b^dagger), which is accurate for states well below the
truncation edge.
"""
import numpy as np
import scipy.sparse as sp
from scipy.special import jv

from .fibre import U01, U1, amp_ratio


def oscillator(N):
    """Annihilation operator b (sparse), eigenvalues w and eigenvectors V of b + b^dagger."""
    b = sp.diags(np.sqrt(np.arange(1, N)), 1, format='csr')
    w, V = np.linalg.eigh((b + b.T).toarray())
    return b, w, V


def function_of_x(vals, V, tol=1e-11):
    """Sparse matrix of f(x) given f at the eigenvalues of b + b^dagger."""
    M = (V * vals) @ V.conj().T
    M[np.abs(M) < tol * np.abs(M).max()] = 0
    return sp.csr_matrix(M)


def projector(i, j, dim=3):
    """Atomic operator |i><j|."""
    return sp.csr_matrix(([1.0], ([i], [j])), shape=(dim, dim))


def profile_ops(p, w, V):
    """Probe (LP1m + eps LP01) and coupling (LP01) fields as operators of x."""
    xa = p['x0a'] * (w + p['xs'])
    u1 = U1[p['m'] - 1]
    fp = amp_ratio(p['m']) * jv(1, u1 * xa) + p['eps'] * jv(0, U01 * xa)
    fcpl = jv(0, U01 * xa)
    return function_of_x(fp, V), function_of_x(fcpl, V)


def recoil_op(eta, u, w, V):
    """Momentum kick exp(i eta u (b + b^dagger)) of a photon emitted along u."""
    return function_of_x(np.exp(1j * eta * u * w), V)
