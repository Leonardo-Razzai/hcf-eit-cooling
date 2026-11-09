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

from .config import me_cutoffs
from .eit import GP, GC, EX, NODES, WEIGHTS, params
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


def banded_assemble(H, Cs, na, N, B):
    """Generator on the coherences |n - m| <= B for H and jumps Cs (dense or sparse) on
    (na internal levels) x (N Fock states). Every term X rho Y restricted to the kept set S is the
    element-wise product X[r_s, r_s'] * Y^T[c_s, c_s'], so M (dimension ~ na^2 N (2B + 1)) is
    assembled in row blocks to bound the temporaries. Returns (M, r, c, n_of): kept (row, col)
    indices and the Fock number of each index."""
    n_of = np.tile(np.arange(N), na)                    # Fock number of each index
    r, c = np.nonzero(np.abs(n_of[:, None] - n_of[None, :]) <= B)   # kept (row, col)
    S = r.size
    H = H.toarray() if sp.issparse(H) else H
    Cd = [C.toarray() if sp.issparse(C) else C for C in Cs]
    K = sum(C.conj().T @ C for C in Cd)
    G = -1j * H - 0.5 * K                               # rho -> G rho + rho G^dag + jumps
    Gc = G.conj()
    Cc = [x.conj() for x in Cd]
    M = np.empty((S, S), complex)
    step = max(1, int(2e7 // S))                        # rows per block (bounded temporaries)
    for i0_ in range(0, S, step):
        sl = slice(i0_, min(S, i0_ + step))
        rs, cs = r[sl], c[sl]
        blk = G[np.ix_(rs, r)] * (cs[:, None] == c[None, :])
        blk += (rs[:, None] == r[None, :]) * Gc[np.ix_(cs, c)]
        for Cj, Ccj in zip(Cd, Cc):
            blk += Cj[np.ix_(rs, r)] * Ccj[np.ix_(cs, c)]
        M[sl] = blk
    return M, r, c, n_of


def banded_generator(p, B=6):
    """Generator of the eliminated master equation restricted to the coherences
    rho_(a,n),(b,m) with |n - m| <= B: d x/dt = M x for the kept elements x.
    Returns (M, r, c, n_of), see banded_assemble."""
    H, Cs, b = effective_ops(p)
    return banded_assemble(H, Cs, 2, p['N'], B)


def solve_banded(p, B=6, rate=True):
    """Steady state and cooling rate of the eliminated master equation keeping
    only the coherences |n - m| <= B (banded_generator), solved with a dense LU.

    The steady state is almost diagonal in n (case A, N = 40: |rho_n,n+k| < 4e-5
    for k >= 4), while the nearly dense recoil operator makes the full
    Liouvillian grow as N^4. Both rates start from |g_p, 1>: W is the integrated
    relaxation of <n>, W0 that of the ground-state population,
    1/W0 = int (P0_ss - P0(t)) dt / P0_ss. In a heavy recoil tail <n> undershoots
    n_ss and rises slowly, which biases W upwards and can make it diverge or
    change sign; the tail enters W0 through its weight, not its energy.
    Returns dict(n, W, W0, P0, pgc, P (Fock populations), dim,
    rho (2N x 2N, zero outside the band)).
    """
    from scipy.linalg import lu_factor, lu_solve
    N = p['N']
    M, r, c, n_of = banded_generator(p, B)
    S = r.size
    diag = np.nonzero(r == c)[0]
    i0 = int(np.nonzero((r == 0) & (c == 0))[0][0])     # |g_p,0><g_p,0|
    M[i0, :] = 0.0
    M[i0, diag] = 1.0                                   # trace condition
    lu = lu_factor(M, overwrite_a=True, check_finite=False)
    rhs = np.zeros(S, complex)
    rhs[i0] = 1.0
    x = lu_solve(lu, rhs)
    nvec = np.zeros(S)
    nvec[diag] = n_of[r[diag]]
    nss = float(np.real(nvec @ x))
    p0vec = np.zeros(S)
    p0vec[diag[n_of[r[diag]] == 0]] = 1.0               # projector on n = 0, both ground levels
    P0 = float(np.real(p0vec @ x))
    W = W0 = np.nan
    if rate:
        x0 = np.zeros(S, complex)
        x0[int(np.nonzero((r == 1) & (c == 1))[0][0])] = 1.0   # |g_p,1><g_p,1|
        src = x0 - x
        src[i0] = 0.0
        y = lu_solve(lu, src)
        W = 1.0 / (-float(np.real(nvec @ y)) / (1.0 - nss))
        W0 = P0 / float(np.real(p0vec @ y))
    rho = np.zeros((2 * N, 2 * N), complex)
    rho[r, c] = x
    P = np.real(np.diag(rho)[:N] + np.diag(rho)[N:])
    return dict(n=nss, W=W, W0=W0, P0=P0, pgc=float(np.real(np.trace(rho[N:, N:]))), P=P,
                dim=S, rho=rho)


def banded_evolve_states(p, pops_list, tmax, B=4, h=0.25, n_out=60):
    """Time evolution of <n> and P0 from the states sum_n pops[n] |g_p,n><g_p,n|, for
    each pops in pops_list (length <= N), in the banded eliminated master equation.
    Fixed-step BDF2, (3 - 2hM) y_k+1 = 4 y_k - y_k-1, with one dense in-place LU
    (memory: one S x S matrix); the start y_-1 = y_0 - h M y_0 keeps second order
    overall (y_-1 = y_0 would add a constant time lag, a global first-order error).
    h in units of 1/nu_r: BDF2 is L-stable, the fast internal rate is only
    Gamma_n ~ nu_r/3 after the elimination. Returns times, <n>, P0, the last two
    of shape (states, times). P0 sums n = 0 over both ground levels."""
    from scipy.linalg import lu_factor, lu_solve
    M, r, c, n_of = banded_generator(p, B)
    S = r.size
    diag = np.nonzero(r == c)[0]
    nvec = np.zeros(S)
    nvec[diag] = n_of[r[diag]]
    p0vec = np.zeros(S)
    p0vec[diag[n_of[r[diag]] == 0]] = 1.0
    y = np.zeros((S, len(pops_list)), complex)
    for m, pops in enumerate(pops_list):
        for n, pn in enumerate(pops):
            y[int(np.nonzero((r == n) & (c == n))[0][0]), m] = pn
    ym = y - h * (M @ y)
    M *= -2 * h
    M[np.diag_indices(S)] += 3.0
    lu = lu_factor(M, overwrite_a=True, check_finite=False)
    nsteps = int(np.ceil(tmax / h))
    every = max(1, nsteps // n_out)
    ts, ns, ps = [0.0], [np.real(nvec @ y)], [np.real(p0vec @ y)]
    for k in range(1, nsteps + 1):
        ym, y = y, lu_solve(lu, 4 * y - ym)
        if k % every == 0:
            ts.append(k * h)
            ns.append(np.real(nvec @ y))
            ps.append(np.real(p0vec @ y))
    return np.array(ts), np.array(ns).T, np.array(ps).T


def banded_evolve(p, pops, tmax, B=4, h=0.25, n_out=60):
    """<n>(t) from one initial state: see banded_evolve_states. Returns times, <n>."""
    t, n, _ = banded_evolve_states(p, [pops], tmax, B=B, h=h, n_out=n_out)
    return t, n[0]


def banded_fock_drift(p, n0s, B=4, h=0.05, t1=100, t2=400):
    """Phonon drift d<n>/dt and energy diffusion d Var(n)/dt (per unit time
    1/nu_r) after preparing |g_p, n0>, from implicit-Euler steps of the banded
    generator (one dense LU); measured between steps t1 and t2, after the
    internal transient. Needs N well above max(n0s)."""
    from scipy.linalg import lu_factor, lu_solve
    M, r, c, n_of = banded_generator(p, B)
    S = r.size
    A = -h * M
    A[np.diag_indices(S)] += 1.0
    del M
    lu = lu_factor(A, overwrite_a=True, check_finite=False)
    diag = np.nonzero(r == c)[0]
    n1 = np.zeros(S)
    n1[diag] = n_of[r[diag]]
    n2 = n1**2
    out = []
    for n0 in n0s:
        y = np.zeros(S, complex)
        y[int(np.nonzero((r == n0) & (c == n0))[0][0])] = 1.0
        mom = []
        for k in range(t2):
            y = lu_solve(lu, y)
            if k in (t1 - 1, t2 - 1):
                m1, m2 = float(np.real(n1 @ y)), float(np.real(n2 @ y))
                mom.append((m1, m2 - m1**2))
        T = (t2 - t1) * h
        out.append(((mom[1][0] - mom[0][0]) / T, (mom[1][1] - mom[0][1]) / T))
    return np.array(out)


def steady(fc, s, rate=True, cut=None, **kw):
    """Converged steady state for the fibre case fc at saturation s: builds
    params(fc, s, **kw), takes the Fock truncation N and the band from the
    first matching rule of numerics.toml (config.me_cutoffs, using the actual
    eta_sp, bp and Delta) and calls solve_banded. cut=(N, band) overrides the
    rule (e.g. a cheap search before the final converged solve). Returns the
    solve_banded dict plus N, band, status."""
    p = params(fc, s=s, **kw)
    if cut is None:
        rule = me_cutoffs(p['eta_sp'], s, p['bp'], kw.get('Delta_G', 3.0))
    else:
        rule = dict(N=cut[0], band=cut[1], status='override')
    p['N'] = rule['N']
    out = solve_banded(p, B=rule['band'], rate=rate)
    out.update(N=rule['N'], band=rule['band'], status=rule['status'])
    return out


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
