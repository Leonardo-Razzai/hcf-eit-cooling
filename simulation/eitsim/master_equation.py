"""Lindblad master equation: Lambda atom x one radial oscillator.

Atomic levels: |g_p>, |g_c> and K excited hyperfine levels (p['levels'];
K = 1 is the closed three-level Lambda).

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
    """Return (L, d, b): sparse Liouvillian, Hilbert dimension (2+K)N, annihilation op.

    Excited level k (index EX + k) has energy offset eps_k from the level the
    lasers are tuned to, relative couplings cp_k, cc_k and branching bp_k into
    g_p. Decays of different excited levels are independent (secular
    approximation, valid for hyperfine splittings >> Gamma).
    """
    N = p['N']
    levels = p.get('levels', [(0.0, 1.0, 1.0, p['bp'])])
    na = 2 + len(levels)
    b, w, V = oscillator(N)
    Im = sp.identity(N, format='csr')
    fp, fcpl = profile_ops(p, w, V)
    H = sp.kron(sp.identity(na), b.T @ b)
    H = H - p['d2'] * sp.kron(projector(GC, GC, na), Im)
    for k, (off, cp, cc, _) in enumerate(levels):
        e = EX + k
        H = H - (p['D'] - off) * sp.kron(projector(e, e, na), Im)
        if cp:
            Hp = 0.5 * p['Wp'] * cp * sp.kron(projector(e, GP, na), fp)
            H = H + Hp + Hp.getH()
        if cc:
            Hc = 0.5 * p['Wc'] * cc * sp.kron(projector(e, GC, na), fcpl)
            H = H + Hc + Hc.getH()
    H = H.tocsr()
    d = na * N
    I = sp.identity(d, format='csr')
    L = -1j * (sp.kron(I, H) - sp.kron(H.T, I))
    for u, wk in zip(NODES, WEIGHTS):
        Dk = recoil_op(p['eta_sp'], u, w, V)
        for k, (_, _, _, bpk) in enumerate(levels):
            for j, bj in ((GP, bpk), (GC, 1 - bpk)):
                if bj:
                    C = np.sqrt(p['G'] * bj * wk) * sp.kron(projector(j, EX + k, na), Dk)
                    L = L + sp.kron(C.conj(), C)
    Pe = sum(projector(EX + k, EX + k, na) for k in range(len(levels)))
    CdC = p['G'] * sp.kron(Pe, Im)
    L = L - 0.5 * (sp.kron(I, CdC) + sp.kron(CdC.T, I))
    return L.tocsc(), d, b


def effective_ops(p):
    """Ground-manifold operators of the master equation with the excited
    manifold adiabatically eliminated (effective-operator formalism, Reiter &
    Sorensen, PRA 85, 032111 (2012)). Returns (H, Cs, b), H and the jump
    operators Cs as sparse (2N x 2N) matrices.

    Valid for |Delta - eps_k + i Gamma/2| >> Om_c, nu_r and the gradient
    coupling. With c_k = 1 / (-(Delta - eps_k) - i Gamma/2) and the excitation
    operator A_k = (Wp cp_k / 2) f_p <g_p| + (Wc cc_k / 2) f_c <g_c| (ground
    manifold -> motion in level k):
        H_eff = H_g - sum_k Re(c_k) A_k^dag A_k,
        C_kjl = sqrt(Gamma b_kj w_l) c_k |g_j> D_l A_k.
    The dark state, annihilated by every A_k, stays dark; the motional energy
    is neglected in the denominators (relative error ~ nu_r / Delta).
    """
    N = p['N']
    levels = p.get('levels', [(0.0, 1.0, 1.0, p['bp'])])
    b, w, V = oscillator(N)
    Im = sp.identity(N, format='csr')
    fp, fcpl = profile_ops(p, w, V)
    H = sp.kron(sp.identity(2), b.T @ b) - p['d2'] * sp.kron(projector(GC, GC, 2), Im)
    A, c = [], []
    for off, cp, cc, _ in levels:
        ck = 1.0 / (-(p['D'] - off) - 0.5j * p['G'])
        Ak = sp.hstack([0.5 * p['Wp'] * cp * fp, 0.5 * p['Wc'] * cc * fcpl]).tocsr()  # N x 2N
        H = H - ck.real * (Ak.getH() @ Ak)
        A.append(Ak)
        c.append(ck)
    d = 2 * N
    Z = sp.csr_matrix((N, d))
    Cs = []
    for u, wk in zip(NODES, WEIGHTS):
        Dk = recoil_op(p['eta_sp'], u, w, V)
        for (_, _, _, bpk), Ak, ck in zip(levels, A, c):
            DA = (Dk @ Ak).tocsr()
            for j, bj in ((GP, bpk), (GC, 1 - bpk)):
                if bj:
                    blocks = [DA, Z] if j == GP else [Z, DA]
                    Cs.append(np.sqrt(p['G'] * bj * wk) * ck * sp.vstack(blocks).tocsr())
    return H.tocsr(), Cs, b


def liouvillian_eliminated(p):
    """Return (L, d, b): sparse Liouvillian of the ground-state master equation
    (excited manifold eliminated, see effective_ops); Hilbert dimension 2N."""
    H, Cs, b = effective_ops(p)
    d = H.shape[0]
    I = sp.identity(d, format='csr')
    L = -1j * (sp.kron(I, H) - sp.kron(H.T, I))
    CdC = sp.csr_matrix((d, d), dtype=complex)
    for C in Cs:
        L = L + sp.kron(C.conj(), C)
        CdC = CdC + C.getH() @ C
    L = L - 0.5 * (sp.kron(I, CdC) + sp.kron(CdC.T, I))
    return L.tocsc(), d, b


def solve_banded(p, B=6, rate=True):
    """Steady state and cooling rate of the eliminated master equation keeping
    only the coherences rho_(a,n),(b,m) with |n - m| <= B.

    The steady state is almost diagonal in n (case A, N = 40: |rho_n,n+k| < 4e-5
    for k >= 4), while the nearly dense recoil operator makes the full
    Liouvillian grow as N^4. Every term X rho Y restricted to the kept set S is
    the element-wise product X[r_s, r_s'] * Y^T[c_s, c_s'], so the restricted
    Liouvillian (dimension ~ 4 N (2B + 1)) is assembled from dense 2N x 2N
    operators and solved with a dense LU. Returns dict(n, W, pgc, P (Fock
    populations), dim, rho (2N x 2N, zero outside the band)).
    """
    from scipy.linalg import lu_factor, lu_solve
    H, Cs, b = effective_ops(p)
    N = p['N']
    n_of = np.tile(np.arange(N), 2)                     # Fock number of each ground index
    r, c = np.nonzero(np.abs(n_of[:, None] - n_of[None, :]) <= B)   # kept (row, col)
    S = r.size
    Hd = H.toarray()
    K = sum((C.getH() @ C) for C in Cs).toarray()
    G = -1j * Hd - 0.5 * K                              # rho -> G rho + rho G^dag + jumps
    dr = (r[:, None] == r[None, :])
    dc = (c[:, None] == c[None, :])
    M = G[np.ix_(r, r)] * dc + dr * G.conj()[np.ix_(c, c)]
    del dr, dc
    for C in Cs:
        Cd = C.toarray()
        M += Cd[np.ix_(r, r)] * Cd.conj()[np.ix_(c, c)]
    diag = np.nonzero(r == c)[0]
    i0 = int(np.nonzero((r == 0) & (c == 0))[0][0])     # |g_p,0><g_p,0|
    M[i0, :] = 0.0
    M[i0, diag] = 1.0                                   # trace condition
    lu = lu_factor(M)
    rhs = np.zeros(S, complex)
    rhs[i0] = 1.0
    x = lu_solve(lu, rhs)
    nvec = np.zeros(S)
    nvec[diag] = n_of[r[diag]]
    nss = float(np.real(nvec @ x))
    W = np.nan
    if rate:
        x0 = np.zeros(S, complex)
        x0[int(np.nonzero((r == 1) & (c == 1))[0][0])] = 1.0   # |g_p,1><g_p,1|
        src = x0 - x
        src[i0] = 0.0
        y = lu_solve(lu, src)
        W = 1.0 / (-float(np.real(nvec @ y)) / (1.0 - nss))
    rho = np.zeros((2 * N, 2 * N), complex)
    rho[r, c] = x
    P = np.real(np.diag(rho)[:N] + np.diag(rho)[N:])
    return dict(n=nss, W=W, pgc=float(np.real(np.trace(rho[N:, N:]))), P=P, dim=S, rho=rho)


def number_op(b, na=3):
    """Phonon number operator on the full (atom x oscillator) space, na atomic levels."""
    return sp.kron(sp.identity(na), b.T @ b).tocsr()


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


def solve_me(p, rate=True, eliminate=False):
    """Steady state of the master equation.

    eliminate=True uses the ground-state master equation with the excited
    manifold adiabatically eliminated (liouvillian_eliminated); pe is then 0.
    Returns dict(n=<n>_ss, W=cooling rate (units of nu_r), pe (total excited
    population), pgc, rho).
    """
    L, d, b = (liouvillian_eliminated if eliminate else liouvillian)(p)
    lu = factorize(L, d)
    rho = steady_state(L, d, lu)
    N = p['N']
    na = d // N
    nop = number_op(b, na)
    pop = [float(np.real(np.trace(rho[i * N:(i + 1) * N, i * N:(i + 1) * N])))
           for i in range(na)]
    W = relax_rate(lu, d, rho, nop, N) if rate else np.nan
    return dict(n=expect(nop, rho), W=W, pe=sum(pop[EX:]), pgc=pop[GC], rho=rho)


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


def fock_drift(p, n0s, h=0.05, t1=100, t2=400, eliminate=False):
    """Phonon drift d<n>/dt (phonons per unit time 1/nu_r) after preparing
    |g_p, n0>, measured between implicit-Euler steps t1 and t2 (after the
    fast atomic transient). eliminate=True: excited manifold eliminated."""
    L, d, b = (liouvillian_eliminated if eliminate else liouvillian)(p)
    nop = number_op(b, d // p['N'])
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
