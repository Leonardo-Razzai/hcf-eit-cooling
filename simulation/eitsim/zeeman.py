"""Zeeman-resolved cooling: eliminated, banded master equation on all ground sublevels.

One radial axis x. A structure (structure_lambda, or a hand-made dict for tests) lists the ground
levels (with their two-photon frame), the excited levels (offset, couplings to the probe and coupling
beams) and the decay groups. A decay group is one distinguishable photon: a matrix D[j, k] of
amplitudes (ground j, excited k), with its own emission pattern along x. Amplitudes of one group
add coherently, which keeps Zeeman coherences correct; different F' are separate groups (their
photons differ by the hyperfine splitting).

Extra fields (repumpers) have their own frequency, so they do not interfere with the cooling beams:
their excitation operators enter H and the jumps as a separate incoherent set.

Internal units: nu_r = 1, as everywhere. Fields are normalised so that the mean |g|^2 over the
sublevels of F for one polarisation component is 1: the Rabi frequencies of eit.params are then the
rms values of the hyperfine-averaged Lambda.
"""
import numpy as np

from .atoms import GROUND_F, LINES, zeeman_coupling
from .eit import NODES, WEIGHTS, params
from .master_equation import banded_assemble
from .operators import oscillator, profile_ops

G_F = {1: -0.5, 2: 0.5}      # Lande factors of 87Rb 5S1/2 (sign and magnitude, nuclear part neglected)


def pattern_nodes(q, axis='perp', K=6, seed=1):
    """Equal-weight nodes of u = k_x/k for a photon of polarisation q (quantisation axis zQ).

    axis='perp': x perpendicular to zQ; 'par': x along zQ. q = 0 has pattern sin^2 th, |q| = 1
    has (1 + cos^2 th)/2, with th measured from zQ. Returns nodes, weights.
    """
    rng = np.random.default_rng(seed)
    n = 1_000_000
    c = rng.uniform(-1, 1, 3 * n)
    w = (1 - c**2) if q == 0 else 0.5 * (1 + c**2)
    c = c[rng.uniform(0, 1, 3 * n) < w][:n]
    if axis == 'par':
        u = c
    else:
        u = np.sqrt(1 - c**2) * np.cos(rng.uniform(0, 2 * np.pi, c.size))
    u = np.sort(u)
    nodes = np.array([b.mean() for b in np.array_split(u, K)])
    return 0.5 * (nodes - nodes[::-1]), np.full(K, 1.0 / K)


def field_matrix(line, F, excited, pol):
    """Couplings g[i, k] from the sublevels of F to the excited levels, for polarisation
    pol = {q: amplitude} (q = me - m on absorption), normalised as stated in the module docstring."""
    ms = list(range(-F, F + 1))
    g = np.zeros((len(ms), len(excited)), complex)
    for i, m in enumerate(ms):
        for k, (Fe, me) in enumerate(excited):
            q = me - m
            if q in pol and abs(me) <= Fe:
                g[i, k] = pol[q] * zeeman_coupling(line, F, m, Fe, me, m - me)
    return g


def _norm(line, F, Fe):
    """Mean over m of |d|^2 for one polarisation component."""
    vals = []
    for q in (-1, 0, 1):
        v = [zeeman_coupling(line, F, m, Fe, m + q, -q)**2 for m in range(-F, F + 1) if abs(m + q) <= Fe]
        vals.append(sum(v) / (2 * F + 1))
    return np.mean(vals)


def structure_lambda(line, Fe, gp, pol_p, pol_c, axis='perp', K=6, fes=None, offsets=None):
    """87Rb Lambda with Zeeman sublevels. Ground order: probe manifold F = gp, then F = 3 - gp.
    fes: excited hyperfine levels included (default only Fe); the lasers are tuned to Fe."""
    gc = 3 - gp
    fes = fes or [Fe]
    hfs = LINES[line]['hfs']
    excited = [(f, m) for f in fes for m in range(-f, f + 1)]
    ground = [(gp, m) for m in range(-gp, gp + 1)] + [(gc, m) for m in range(-gc, gc + 1)]
    nP = 2 * gp + 1
    gP = field_matrix(line, gp, excited, pol_p) / np.sqrt(_norm(line, gp, Fe))
    gC = field_matrix(line, gc, excited, pol_c) / np.sqrt(_norm(line, gc, Fe))
    # eps in Hz; converted by the caller through om
    eps = np.array([hfs[f] - hfs[Fe] for f, _ in excited], float)
    groups = []
    for f in fes:
        idx = [k for k, e in enumerate(excited) if e[0] == f]
        norm = sum(zeeman_coupling(line, F, m, f, 0 if f else 0, m)**2 for F in GROUND_F for m in range(-F, F + 1) if abs(m) <= 1)
        for q in (-1, 0, 1):
            D = np.zeros((len(ground), len(excited)))
            for j, (F, m) in enumerate(ground):
                for k in idx:
                    Fk, me = excited[k]
                    if m == me + q:
                        D[j, k] = zeeman_coupling(line, F, m, Fk, me, q) / np.sqrt(norm)
            u, w = pattern_nodes(q, axis, K)
            groups.append((D, u, w))
    return dict(ground=ground, excited=excited, nP=nP, gP=gP, gC=gC, eps=eps, groups=groups, line=line, Fe=Fe, gp=gp,
                gF=np.array([G_F[F] for F, _ in ground]), mF=np.array([m for _, m in ground]), FF=np.array([F for F, _ in ground]))


def structure_zeeman_lambda(line='D1', Fe=2, K=6):
    """Lambda inside F = 2, after Qiu et al., PRX Quantum 7, 020358 (2026), with guided light and
    the field along the fibre: |2,2> -(probe sigma-)-> |F',1> <-(coupling sigma+)- |2,0>. |2,2> is dark
    to the coupling beam and is the target. F = 1 is populated by decay and needs a repumper.
    Frame: the two beams differ by the Zeeman splitting, which removes the linear Zeeman shift, so
    only the two-photon detuning remains, H = -d2 (2 - m)/2 on F = 2."""
    st = structure_lambda(line, Fe, 2, {-1: 1.0}, {1: 1.0}, K=K)
    na, nE = len(st['ground']), len(st['excited'])
    nrm = np.sqrt(_norm(line, 2, Fe))
    gP, gC = np.zeros((na, nE), complex), np.zeros((na, nE), complex)
    gP[:5] = field_matrix(line, 2, st['excited'], {-1: 1.0}) / nrm
    gC[:5] = field_matrix(line, 2, st['excited'], {1: 1.0}) / nrm
    st.update(gPf=gP, gCf=gC, frame='zeeman')
    return st


def structure_three_level(bp):
    """The closed three-level Lambda written as a structure (tests)."""
    return dict(ground=[('p', 0), ('c', 0)], excited=[('e', 0)], nP=1, gP=np.ones((1, 1), complex), gC=np.ones((1, 1), complex),
                eps=np.zeros(1), groups=[(np.array([[np.sqrt(bp)], [0.0]]), NODES, WEIGHTS),
                                         (np.array([[0.0], [np.sqrt(1 - bp)]]), NODES, WEIGHTS)],
                gF=np.zeros(2), mF=np.zeros(2), FF=np.zeros(2))


def _final_sets(st):
    """Ground indices reached by one photon: a photon of given polarisation and direction leaves a
    coherent superposition within one ground hyperfine level (photons ending in different F differ
    by the hyperfine splitting, so they are distinguishable). The three-level test case has one
    level per set."""
    sets = {}
    for j, (F, _) in enumerate(st['ground']):
        sets.setdefault(F, []).append(j)
    return list(sets.values())


def ops(st, p, b_vec=(0.0, 0.0), repump=()):
    """H and jump operators on (ground sublevels) x (Fock), excited levels eliminated.

    b_vec = (L_z, L_x): Larmor frequencies along zQ and across it, H_Z = g_F (2 L) . F, so that
    F = 2 precesses at L and F = 1 at -L, in units of nu_r.
    repump: list of dict(gm=complex matrix [ground, excited], W=rms Rabi, detuning=from excited level):
    a separate frequency, with the coupling profile (LP01)."""
    N = p['N']
    b, w, V = oscillator(N)
    fp, fcpl = profile_ops(p, w, V)
    fp, fcpl = fp.toarray(), fcpl.toarray()
    na = len(st['ground'])
    d = na * N
    In = np.eye(N)
    H = np.kron(np.eye(na), np.diag(np.arange(N, dtype=float))).astype(complex)
    if st.get('frame') == 'zeeman':
        for i, (Fi, mi) in enumerate(st['ground']):
            if Fi == 2:
                H[i * N:(i + 1) * N, i * N:(i + 1) * N] -= 0.5 * p['d2'] * (2 - mi) * In
    else:
        for i in range(st['nP'], na):
            H[i * N:(i + 1) * N, i * N:(i + 1) * N] -= p['d2'] * In
    Lz, Lx = b_vec
    if Lz or Lx:
        for i, (Fi, mi) in enumerate(st['ground']):
            H[i * N:(i + 1) * N, i * N:(i + 1) * N] += 2 * Lz * st['gF'][i] * mi * In
            for j, (Fj, mj) in enumerate(st['ground']):
                if Fi == Fj and abs(mi - mj) == 1 and Lx:
                    fx = 0.5 * np.sqrt(Fi * (Fi + 1) - mi * mj)
                    H[i * N:(i + 1) * N, j * N:(j + 1) * N] += 2 * Lx * st['gF'][i] * fx * In
    nE = len(st['excited'])
    off = 2 * np.pi * st['eps'] / p.get('om', 1.0)
    c = 1.0 / (-(p['D'] - off) - 0.5j * p['G'])
    Ak = []
    for k in range(nE):
        A = np.zeros((N, d), complex)
        for i in range(na):
            if 'gPf' in st:
                A[:, i * N:(i + 1) * N] = (0.5 * p['Wp'] * np.conj(st['gPf'][i, k]) * fp
                                           + 0.5 * p['Wc'] * np.conj(st['gCf'][i, k]) * fcpl)
            elif i < st['nP']:
                A[:, i * N:(i + 1) * N] = 0.5 * p['Wp'] * np.conj(st['gP'][i, k]) * fp
            else:
                A[:, i * N:(i + 1) * N] = 0.5 * p['Wc'] * np.conj(st['gC'][i - st['nP'], k]) * fcpl
        Ak.append(A)
        H = H - c[k].real * (A.conj().T @ A)
    Cs = []
    for D, u, wg in st['groups']:
        for uk, wk in zip(u, wg):
            Dk = (V * np.exp(1j * p['eta_sp'] * uk * w)) @ V.conj().T
            for jj in _final_sets(st):
                C = np.zeros((d, d), complex)
                used = False
                for j in jj:
                    if np.any(D[j]):
                        C[j * N:(j + 1) * N, :] = np.sqrt(p['G'] * wk) * Dk @ sum(D[j, k] * c[k] * Ak[k] for k in range(nE) if D[j, k])
                        used = True
                if used:
                    Cs.append(C)
    for rp in repump:
        cr = 1.0 / (-(rp['detuning']) - 0.5j * p['G'])
        ARs = []
        for k in range(nE):
            A = np.zeros((N, d), complex)
            for i in range(na):
                if rp['gm'][i, k]:
                    A[:, i * N:(i + 1) * N] = 0.5 * rp['W'] * np.conj(rp['gm'][i, k]) * fcpl
            ARs.append(A)
            H = H - cr.real * (A.conj().T @ A)
        for D, u, wg in st['groups']:
            for uk, wk in zip(u, wg):
                Dk = (V * np.exp(1j * p['eta_sp'] * uk * w)) @ V.conj().T
                for jj in _final_sets(st):
                    C = np.zeros((d, d), complex)
                    used = False
                    for j in jj:
                        if np.any(D[j]):
                            C[j * N:(j + 1) * N, :] = np.sqrt(p['G'] * wk) * cr * Dk @ sum(D[j, k] * ARs[k] for k in range(nE) if D[j, k])
                            used = True
                    if used:
                        Cs.append(C)
    return H, Cs


def solve(st, p, B=3, b_vec=(0.0, 0.0), repump=(), rate=True):
    """Steady state and rates. Rates start from n = 1 with the steady internal state.
    Returns n, P0, W, W0, P (Fock), pint (population per ground sublevel), cond (largest
    relative residual, large if there is no unique steady state)."""
    from scipy.linalg import lu_factor, lu_solve
    N = p['N']
    na = len(st['ground'])
    H, Cs = ops(st, p, b_vec, repump)
    M, r, c, n_of = banded_assemble(H, Cs, na, N, B)
    S = r.size
    diag = np.nonzero(r == c)[0]
    Morig_diag = None
    i0 = int(np.nonzero((r == 0) & (c == 0))[0][0])
    M[i0, :] = 0.0
    M[i0, diag] = 1.0
    lu = lu_factor(M, overwrite_a=True, check_finite=False)
    rhs = np.zeros(S, complex)
    rhs[i0] = 1.0
    x = lu_solve(lu, rhs, check_finite=False)
    pops = np.real(x[diag])
    if not np.all(np.isfinite(x)) or abs(pops.sum() - 1) > 1e-6 or pops.min() < -1e-3:
        # an unmixed dark state: no unique steady state
        return dict(n=np.nan, P0=np.nan, P=None, pint=np.full(na, np.nan), S=S, W=np.nan, W0=np.nan)
    idx = r[diag]
    nvec = n_of[idx].astype(float)
    Pf = np.array([pops[n_of[idx] == m].sum() for m in range(N)])
    lev = idx // N
    pint = np.array([pops[lev == a].sum() for a in range(na)])
    out = dict(n=float(nvec @ pops), P0=float(Pf[0]), P=Pf, pint=pint, S=S, W=np.nan, W0=np.nan)
    if rate:
        x0 = np.zeros(S, complex)
        for a in range(na):
            x0[int(np.nonzero((r == a * N + 1) & (c == a * N + 1))[0][0])] = pint[a]
        src = x0 - x
        src[i0] = 0.0
        y = lu_solve(lu, src, check_finite=False)
        ny = float(np.real(nvec @ y[diag]))
        p0y = float(np.real(y[diag][n_of[idx] == 0].sum()))
        out['W'] = -(1 - out['n']) / ny
        out['W0'] = out['P0'] / p0y
    return out


def setup(fc, s, Delta_G, d2, N, wc_scale=1.0, line=None):
    """Parameters for the structures, with the design laser values of eit.params."""
    p = params(fc, s=s, Delta_G=Delta_G, d2=d2, N=N)
    p['Wc'] = p['Wc'] * wc_scale
    p['om'] = fc['om']
    return p
