"""Lindblad master equation: Lambda atom x one radial oscillator.

Column-stacking vectorisation, vec(A rho B) = (B^T kron A) vec(rho).
Solves for the steady state, the cooling rate (inverse integrated relaxation
time of <n>) and time evolution (fixed-step BDF2).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .eit import GP, GC, EX, NODES, WEIGHTS
from .operators import oscillator, projector, profile_ops, recoil_op


def liouvillian(p):
    """Return (L, d, b): sparse Liouvillian, Hilbert dimension 3N, annihilation op."""
    N = p['N']
    b, w, V = oscillator(N)
    Im = sp.identity(N, format='csr')
    fp, fcpl = profile_ops(p, w, V)
    H = sp.kron(sp.identity(3), b.T @ b)
    H = H - p['D'] * sp.kron(projector(EX, EX), Im) - p['d2'] * sp.kron(projector(GC, GC), Im)
    Hp = 0.5 * p['Wp'] * sp.kron(projector(EX, GP), fp)
    Hc = 0.5 * p['Wc'] * sp.kron(projector(EX, GC), fcpl)
    H = (H + Hp + Hp.getH() + Hc + Hc.getH()).tocsr()
    d = 3 * N
    I = sp.identity(d, format='csr')
    L = -1j * (sp.kron(I, H) - sp.kron(H.T, I))
    branch = {GP: p['bp'], GC: 1 - p['bp']}
    for u, wk in zip(NODES, WEIGHTS):
        Dk = recoil_op(p['eta_sp'], u, w, V)
        for j, bj in branch.items():
            C = np.sqrt(p['G'] * bj * wk) * sp.kron(projector(j, EX), Dk)
            L = L + sp.kron(C.conj(), C)
    CdC = p['G'] * sp.kron(projector(EX, EX), Im)
    L = L - 0.5 * (sp.kron(I, CdC) + sp.kron(CdC.T, I))
    return L.tocsc(), d, b


def number_op(b):
    """Phonon number operator on the full (atom x oscillator) space."""
    return sp.kron(sp.identity(3), b.T @ b).tocsr()


def expect(op, rho):
    return float(np.real(np.sum(op.multiply(rho.T))))


def factorize(L, d):
    """Sparse LU of L with the first row replaced by the trace condition."""
    tr = sp.csr_matrix((np.ones(d), (np.zeros(d, int), np.arange(d) * (d + 1))),
                       shape=(1, d * d))
    M = sp.vstack([tr, L.tocsr()[1:, :]]).tocsc()
    return spla.splu(M, permc_spec='MMD_AT_PLUS_A')


def steady_state(L, d, lu=None):
    lu = lu or factorize(L, d)
    rhs = np.zeros(d * d, complex)
    rhs[0] = 1
    rho = lu.solve(rhs).reshape((d, d), order='F')
    return 0.5 * (rho + rho.conj().T)


def relax_rate(lu, d, rho_ss, nop, N, n0=1):
    """Cooling rate W = (n0 - n_ss) / int_0^inf (<n>(t) - n_ss) dt after
    preparing |g_p, n0>; one extra solve with the same factorisation.
    Unlike the slowest Liouvillian eigenvalue it is not affected by slow
    modes at the truncation edge."""
    rho0 = np.zeros((d, d), complex)
    rho0[GP * N + n0, GP * N + n0] = 1
    c = (rho0 - rho_ss).flatten(order='F')
    c[0] = 0
    y = lu.solve(c).reshape((d, d), order='F')
    tau = -expect(nop, y) / (n0 - expect(nop, rho_ss))
    return 1.0 / tau


def solve_me(p, rate=True):
    """Steady state of the master equation.

    Returns dict(n=<n>_ss, W=cooling rate (units of nu_r), pe, pgc, rho).
    """
    L, d, b = liouvillian(p)
    lu = factorize(L, d)
    rho = steady_state(L, d, lu)
    nop = number_op(b)
    N = p['N']
    pop = [float(np.real(np.trace(rho[i * N:(i + 1) * N, i * N:(i + 1) * N])))
           for i in range(3)]
    W = relax_rate(lu, d, rho, nop, N) if rate else np.nan
    return dict(n=expect(nop, rho), W=W, pe=pop[EX], pgc=pop[GC], rho=rho)


def evolve(L, rho0, tmax, observable, h=0.02, n_out=60):
    """Integrate d rho/dt = L rho with fixed-step BDF2 (L-stable for the fast
    optical modes); h in units of 1/nu_r. Returns times and <observable>."""
    d = rho0.shape[0]
    nsteps = int(np.ceil(tmax / h))
    lu2 = spla.splu((3 * sp.identity(d * d, format='csc') - 2 * h * L).tocsc())
    lu1 = spla.splu((sp.identity(d * d, format='csc') - h * L).tocsc())
    y0 = rho0.flatten(order='F')
    every = max(1, nsteps // n_out)
    ts, vals = [0.0], [expect(observable, rho0)]
    ym, y = y0, lu1.solve(y0)                     # implicit Euler start
    for i in range(2, nsteps + 1):
        ym, y = y, lu2.solve(4 * y - ym)
        if i % every == 0:
            ts.append(i * h)
            vals.append(expect(observable, y.reshape((d, d), order='F')))
    return np.array(ts), np.array(vals)


def fock_drift(p, n0s, h=0.05, t1=100, t2=400):
    """Phonon drift d<n>/dt (phonons per unit time 1/nu_r) after preparing
    |g_p, n0>, measured between implicit-Euler steps t1 and t2 (after the
    fast atomic transient)."""
    L, d, b = liouvillian(p)
    nop = number_op(b)
    lu = spla.splu((sp.identity(d * d, format='csc') - h * L).tocsc(),
                   permc_spec='MMD_AT_PLUS_A')
    out = []
    for n0 in n0s:
        r0 = np.zeros((d, d), complex)
        r0[n0, n0] = 1
        y = r0.flatten(order='F')
        ns = []
        for k in range(t2):
            y = lu.solve(y)
            if k in (t1 - 1, t2 - 1):
                ns.append(expect(nop, y.reshape((d, d), order='F')))
        out.append((ns[1] - ns[0]) / ((t2 - t1) * h))
    return np.array(out)
