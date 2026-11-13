"""Both radial axes with the Zeeman structure of 87Rb: eliminated master equation on
(8 ground sublevels) x (n_x, n_y).

Fields. Every laser field has a frequency offset `delta` (from the zero-field transition of the
manifold it drives, in units of nu_r), the ground manifold F it drives, and one or more components
(Cartesian polarisation in the fibre frame, transverse profile): 'lp01' = f_c(x) f_c(y), 'x' = f_p(x),
'y' = f_p(y), with f_p the LP11 radial profile and f_c the LP01 one (operators.profile_ops; near the
axis LP11 and TE01 are linear in x, y, as in studies/two_axes.py). A TE01 probe is the component pair
(x-polarised, -f_p(y)) + (y-polarised, f_p(x)). Polarisations are converted to spherical components
about the field direction B ('z' = fibre axis, 'x' = transverse).

Coherence. The excited states are eliminated (|Delta| >> everything else, as in zeeman.py). An
excitation path (ground g, field f, excited e) oscillates at nu = E_g + delta_f, with E_g = g_F m Z the
linear Zeeman energy (Z = mu_B B / hbar in units of nu_r). Paths with equal nu (within `tol`) interfere
and form a class; terms between different classes oscillate at their frequency difference and are
dropped (rotating-wave approximation; needs the Zeeman splittings >> tol >> two-photon detunings).
Fields marked incoherent (repumpers) never interfere with the others. One rotating frame theta_g per
ground sublevel makes every kept term static: theta_g + delta_f = Omega_C for all paths of class C
(least squares; the residual is reported and is zero for the schemes in `schemes`). The ground
Hamiltonian is then (E_g - theta_g) + n_x + n_y. Coherences created by a jump between final sublevels
whose frames differ by more than `tol` are dropped (secular).

Decay. Each class, decay polarisation q and final hyperfine level F gives a jump operator, with the
dipole pattern of q about B (sin^2 for pi, (1 + cos^2)/2 for sigma) and the recoil exp(i eta (u_x x + u_y y))
on both axes (quadrature: Gauss-Legendre in cos theta x uniform azimuth). Rabi frequencies are
normalised as in zeeman.py (rms over the sublevels of the hyperfine-averaged Lambda).

Solvers. `evolve`: RK4 on the dense density matrix (dimension 8 N_x N_y), with optional stages
(alternating field sets, e.g. LP11 orientation switching). `steady_direct`: dense Liouvillian, only
for small dimensions (tests). Internal units nu_r = 1. Intended for the cluster at N >= 6 per axis
(see simulation/cluster/README.md); `cost` estimates the work.
"""
import time

import numpy as np

from .atoms import LINES, zeeman_coupling
from .operators import oscillator, profile_ops
from .zeeman import G_F, _norm

GROUND = [(1, m) for m in (-1, 0, 1)] + [(2, m) for m in range(-2, 3)]
S2 = 2**-0.5


def quant_basis(bdir):
    """Rows: x_Q, y_Q, z_Q in fibre coordinates (z_Q along B)."""
    if bdir == 'z':
        return np.eye(3)
    if bdir == 'x':
        return np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    raise ValueError(bdir)


def spherical(vec, bdir):
    """Spherical components E_q = e_q^* . E about B, q = m_e - m_g on absorption,
    e_{+1} = -(x_Q + i y_Q)/sqrt 2, e_{-1} = (x_Q - i y_Q)/sqrt 2, e_0 = z_Q."""
    ex, ey, ez = quant_basis(bdir) @ np.asarray(vec, complex)
    return {1: -(ex - 1j * ey) * S2, -1: (ex + 1j * ey) * S2, 0: ez}


def emission_nodes(q, bdir, n_theta=4, n_phi=8):
    """Nodes (u_x, u_y) of the emission direction projected on the fibre transverse axes, and weights,
    for a photon of polarisation q about B."""
    c, wc = np.polynomial.legendre.leggauss(n_theta)
    pat = (1 - c**2) if q == 0 else 0.5 * (1 + c**2)
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    C, P = np.meshgrid(c, ph, indexing='ij')
    W = (wc * pat)[:, None] * np.ones_like(P)
    s = np.sqrt(1 - C**2)
    nQ = np.stack([s * np.cos(P), s * np.sin(P), C], axis=-1)            # in the B frame
    nF = nQ @ quant_basis(bdir)                                          # fibre frame
    W = W / W.sum()
    return nF[..., 0].ravel(), nF[..., 1].ravel(), W.ravel()


def field(F, comps, amp, delta=0.0, stage=None, incoherent=False, name=''):
    """A laser field: comps = [(polarisation vector in the fibre frame, profile, coefficient)]."""
    return dict(F=F, comps=comps, amp=amp, delta=delta, stage=stage, incoherent=incoherent, name=name)


class Model:
    def __init__(self, fields, p, Nx, Ny, bdir='z', zeeman=0.0, fes=(2,), Fe_ref=None, line='D1', tol=2.0,
                 nodes=(4, 8), stages=1, exc_shift=None):
        self.p, self.Nx, self.Ny, self.bdir, self.Z, self.line = p, Nx, Ny, bdir, zeeman, line
        self.fes = tuple(fes)
        self.Fe_ref = Fe_ref if Fe_ref is not None else self.fes[0]
        self.tol, self.stages = tol, stages
        self.exc_shift = exc_shift or {}          # {(Fe, me): Hz}, e.g. the tensor light shift of the lattice
        self.na = len(GROUND)
        self.M = Nx * Ny
        self.d = self.na * self.M
        self.excited = [(f, m) for f in self.fes for m in range(-f, f + 1)]
        self.E = np.array([G_F[F] * m * zeeman for F, m in GROUND])
        # motional operators
        _, wx, Vx = oscillator(Nx)
        _, wy, Vy = oscillator(Ny)
        fpx, fcx = (o.toarray() for o in profile_ops(p, wx, Vx))
        fpy, fcy = (o.toarray() for o in profile_ops(p, wy, Vy))
        Ix, Iy = np.eye(Nx), np.eye(Ny)
        self.prof = {'x': np.kron(fpx, Iy), 'y': np.kron(Ix, fpy), 'lp01': np.kron(fcx, fcy)}
        self.nvec_x = np.kron(np.arange(Nx, dtype=float), np.ones(Ny))
        self.nvec_y = np.kron(np.ones(Nx), np.arange(Ny, dtype=float))
        self.wx, self.Vx, self.wy, self.Vy = wx, Vx, wy, Vy
        self.nodes = nodes
        self.fields = fields
        self._paths()
        self._classes()
        self._frame()
        self.ops = [self._build(s) for s in range(stages)]

    # ------------------------------------------------------------ paths, classes, frame
    def _paths(self):
        nrm = {F: np.sqrt(_norm(self.line, F, self.Fe_ref)) for F in (1, 2)}
        self.paths = []
        for fi, f in enumerate(self.fields):
            comps = [(spherical(v, self.bdir), prof, coef) for v, prof, coef in f['comps']]
            for gi, (F, m) in enumerate(GROUND):
                if F != f['F']:
                    continue
                for ei, (Fe, me) in enumerate(self.excited):
                    q = me - m
                    if abs(q) > 1:
                        continue
                    dip = zeeman_coupling(self.line, F, m, Fe, me, m - me) / nrm[F]
                    if dip == 0:
                        continue
                    op = sum(coef * Eq[q] * self.prof[prof] for Eq, prof, coef in comps if Eq[q] != 0)
                    if np.isscalar(op) or not np.any(op):
                        continue
                    self.paths.append(dict(g=gi, f=fi, e=ei, op=0.5 * f['amp'] * dip * op,
                                           nu=self.E[gi] + f['delta'], inc=f['incoherent']))

    def _classes(self):
        """Cluster paths by frequency; incoherent fields cluster only among their own paths."""
        groups = {}
        for k, pa in enumerate(self.paths):
            groups.setdefault(pa['f'] if pa['inc'] else -1, []).append(k)
        self.classes = []
        for key, ks in groups.items():
            ks = sorted(ks, key=lambda k: self.paths[k]['nu'])
            cur = [ks[0]]
            for k in ks[1:]:
                if self.paths[k]['nu'] - self.paths[cur[-1]]['nu'] <= self.tol:
                    cur.append(k)
                else:
                    self.classes.append(cur)
                    cur = [k]
            self.classes.append(cur)
        for ci, ks in enumerate(self.classes):
            for k in ks:
                self.paths[k]['c'] = ci
        spread = max(np.ptp([self.paths[k]['nu'] for k in ks]) for ks in self.classes)
        self.class_spread = float(spread)

    def _frame(self):
        """theta_g + delta_f - Omega_C = 0 for every path (least squares), gauge pinned weakly to
        theta_g = E_g."""
        nC = len(self.classes)
        rows, rhs = [], []
        for pa in self.paths:
            r = np.zeros(self.na + nC)
            r[pa['g']] = 1.0
            r[self.na + pa['c']] = -1.0
            rows.append(r)
            rhs.append(-self.fields[pa['f']]['delta'])
        eps = 1e-6
        for g in range(self.na):
            r = np.zeros(self.na + nC)
            r[g] = eps
            rows.append(r)
            rhs.append(eps * self.E[g])
        A, b = np.array(rows), np.array(rhs)
        x = np.linalg.lstsq(A, b, rcond=None)[0]
        self.theta = x[:self.na]
        res = A[:len(self.paths)] @ x - b[:len(self.paths)]
        self.frame_mismatch = float(np.abs(res).max()) if len(res) else 0.0
        self.Delta = self.E - self.theta

    # ------------------------------------------------------------ operators
    def _build(self, stage):
        p, M, na, d = self.p, self.M, self.na, self.d
        active = [fi for fi, f in enumerate(self.fields) if f['stage'] is None or f['stage'] == stage]
        hfs = LINES[self.line]['hfs']
        om = p.get('om', 1.0)
        off = np.array([2 * np.pi * (hfs[Fe] - hfs[self.Fe_ref] + self.exc_shift.get((Fe, me), 0.0)) / om for Fe, me in self.excited])
        cexc = 1.0 / (-(p['D'] - off) - 0.5j * p['G'])
        H = np.kron(np.diag(self.Delta), np.eye(M)).astype(complex)
        H += np.kron(np.eye(na), np.diag(self.nvec_x + self.nvec_y))
        # decay amplitudes, normalised to the total decay of one excited sublevel
        dec = {}
        for Fe in self.fes:
            dec[Fe] = sum(zeeman_coupling(self.line, F, m, Fe, 0, m)**2 for F, m in GROUND if abs(m) <= 1)
        jumps = {q: [] for q in (-1, 0, 1)}
        K = np.zeros((d, d), complex)
        for ci, ks in enumerate(self.classes):
            ks = [k for k in ks if self.paths[k]['f'] in active]
            if not ks:
                continue
            cols = sorted({self.paths[k]['g'] for k in ks})
            cidx = {g: i for i, g in enumerate(cols)}
            A = {}
            for k in ks:
                pa = self.paths[k]
                blk = A.setdefault(pa['e'], np.zeros((M, len(cols) * M), complex))
                j = cidx[pa['g']]
                blk[:, j * M:(j + 1) * M] += pa['op']
            colsel = np.concatenate([np.arange(g * M, (g + 1) * M) for g in cols])
            for e, Ae in A.items():
                Hk = -cexc[e].real * (Ae.conj().T @ Ae)
                H[np.ix_(colsel, colsel)] += Hk
            for q, F, Fg in [(q, F, Fg) for q in (-1, 0, 1) for F in (1, 2) for Fg in self.fes]:
                if True:
                    # photons from different excited hyperfine levels differ by >> Gamma: separate decay channels
                    finals = {}
                    for e, Ae in A.items():
                        Fe, me = self.excited[e]
                        if Fe != Fg:
                            continue
                        m = me + q
                        if abs(m) > F:
                            continue
                        j = GROUND.index((F, m))
                        amp = zeeman_coupling(self.line, F, m, Fe, me, q) / np.sqrt(dec[Fe])
                        if amp == 0:
                            continue
                        finals.setdefault(j, []).append(np.sqrt(p['G']) * amp * cexc[e] * Ae)
                    if not finals:
                        continue
                    # secular grouping of the final sublevels by their frame
                    js = sorted(finals, key=lambda j: self.theta[j])
                    grp = [[js[0]]]
                    for j in js[1:]:
                        if abs(self.theta[j] - self.theta[grp[-1][-1]]) <= self.tol:
                            grp[-1].append(j)
                        else:
                            grp.append([j])
                    for gj in grp:
                        J = np.vstack([sum(finals[j]) for j in gj])               # (len(gj) M) x (len(cols) M)
                        rowsel = np.concatenate([np.arange(j * M, (j + 1) * M) for j in gj])
                        jumps[q].append((rowsel, colsel, J))
                        K[np.ix_(colsel, colsel)] += J.conj().T @ J
        G = -1j * H - 0.5 * K
        rec = {}
        for q in (-1, 0, 1):
            ux, uy, w = emission_nodes(q, self.bdir, *self.nodes)
            eta = p['eta_sp']
            Rs = []
            for a, b in zip(ux, uy):
                Dx = (self.Vx * np.exp(1j * eta * a * self.wx)) @ self.Vx.conj().T
                Dy = (self.Vy * np.exp(1j * eta * b * self.wy)) @ self.Vy.conj().T
                Rs.append((Dx, Dy))
            RR = np.array([np.kron(Dx, Dy) for Dx, Dy in Rs])                    # (K, M, M)
            rec[q] = (Rs, w, RR.reshape(-1, self.M), (w[:, None, None] * RR.conj().transpose(0, 2, 1)).reshape(-1, self.M))
        return dict(H=H, K=K, G=G, Gd=G.conj().T, jumps=jumps, rec=rec)

    # ------------------------------------------------------------ dynamics
    def rhs(self, r, st):
        out = st['G'] @ r + r @ st['Gd']
        na, M = self.na, self.M
        for q in (-1, 0, 1):
            if not st['jumps'][q]:
                continue
            Y = np.zeros((self.d, self.d), complex)
            for rowsel, colsel, J in st['jumps'][q]:
                Y[np.ix_(rowsel, rowsel)] += J @ r[np.ix_(colsel, colsel)] @ J.conj().T
            out += self.recoil(Y, st['rec'][q])
        return out

    def recoil(self, Y, rec):
        """sum_u w_u R_u Y R_u^dag with R_u = 1_internal x Dx(u_x) x Dy(u_y), as two large matrix products:
        left, the stacked R_u (K M x M) times Y; right, sum_u (R_u Y) (w_u R_u^dag) = [R_1 Y ... R_K Y] times
        the stacked w_u R_u^dag."""
        na, M = self.na, self.M
        RL, RH = rec[2], rec[3]
        K = RL.shape[0] // M
        Yr = Y.reshape(na, M, na, M).transpose(1, 0, 2, 3).reshape(M, na * na * M)     # (i, (a, b, j))
        T = (RL @ Yr).reshape(K, M, na * na, M)                                          # (k, i', ab, j)
        T = T.transpose(1, 2, 0, 3).reshape(M * na * na, K * M)                          # ((i', ab), (k, j))
        acc = (T @ RH).reshape(M, na, na, M)                                             # (i', a, b, j')
        return acc.transpose(1, 0, 2, 3).reshape(self.d, self.d)

    def observables(self, r):
        pop = np.real(np.diag(r)).reshape(self.na, self.M)
        pm = pop.sum(axis=0)
        return dict(P00=float(pm[0]), nx=float(pm @ self.nvec_x), ny=float(pm @ self.nvec_y), pint=pop.sum(axis=1),
                    P0x=float(pm[self.nvec_x == 0].sum()), P0y=float(pm[self.nvec_y == 0].sum()), trace=float(pop.sum()))

    def initial(self, nx=1, ny=1, internal=None):
        """|nx, ny> times a diagonal internal state (default: uniform over F = 2)."""
        if internal is None:
            internal = np.array([0.0 if F == 1 else 1.0 for F, _ in GROUND])
        internal = np.asarray(internal, float) / np.sum(internal)
        r = np.zeros((self.d, self.d), complex)
        k = nx * self.Ny + ny
        for g in range(self.na):
            r[g * self.M + k, g * self.M + k] = internal[g]
        return r

    def evolve(self, r0, tmax, dt=None, tau=None, record_every=1, log=None):
        """RK4. tau: stage length (stages cycle 0, 1, ...). Returns times, P00(t), final state."""
        dt = dt or min(0.1, 1.2 / (self.Nx + self.Ny))
        nsteps = int(round(tmax / dt))
        per = None if tau is None else max(1, int(round(tau / dt)))
        r = r0.copy()
        ts, P00 = [], []
        t0 = time.time()
        for k in range(nsteps):
            st = self.ops[0] if per is None else self.ops[(k // per) % self.stages]
            k1 = self.rhs(r, st)
            k2 = self.rhs(r + 0.5 * dt * k1, st)
            k3 = self.rhs(r + 0.5 * dt * k2, st)
            k4 = self.rhs(r + dt * k3, st)
            r = r + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if k % record_every == 0:
                ts.append((k + 1) * dt)
                P00.append(sum(np.real(r[g * self.M, g * self.M]) for g in range(self.na)))
            if log and k % max(1, nsteps // 20) == 0:
                log('  step %d / %d, t %.0f, P00 %.4f (%.0f s)' % (k, nsteps, (k + 1) * dt, P00[-1], time.time() - t0))
        return np.array(ts), np.array(P00), r

    # ------------------------------------------------------------ small systems
    def liouvillian(self, stage=0):
        st = self.ops[stage]
        d = self.d
        I = np.eye(d)
        L = np.kron(I, st['G']) + np.kron(st['G'].conj(), I)              # vec(G r + r G^dag), column stacking
        for q in (-1, 0, 1):
            Rs, w = st['rec'][q][:2]
            for rowsel, colsel, J in st['jumps'][q]:
                Jf = np.zeros((d, d), complex)
                Jf[np.ix_(rowsel, colsel)] = J
                for (Dx, Dy), wk in zip(Rs, w):
                    RJ = np.kron(np.eye(self.na), np.kron(Dx, Dy)) @ Jf
                    L += wk * np.kron(RJ.conj(), RJ)
        return L

    def steady_direct(self, stage=0, rate=True):
        """Steady state and W00 (from |1,1> with the steady internal populations), dense; d^2 <= 1e4."""
        d = self.d
        if d * d > 10000:
            raise ValueError('steady_direct is for small systems (d = %d)' % d)
        L = self.liouvillian(stage)
        tr = np.zeros(d * d)
        tr[np.arange(d) * (d + 1)] = 1.0
        A = L.copy()
        A[0, :] = tr
        b = np.zeros(d * d, complex)
        b[0] = 1.0
        rho = np.linalg.solve(A, b).reshape(d, d, order='F')
        out = self.observables(rho)
        out['rho'] = rho
        if rate:
            r0 = self.initial(1, 1 if self.Ny > 1 else 0, out['pint'])
            src = (r0 - rho).reshape(-1, order='F')
            src[0] = 0.0
            y = np.linalg.solve(A, src).reshape(d, d, order='F')        # y = L^-1 (r0 - rho): int (rho(t) - rho) dt = -y
            integ = -self.observables(y)['P00']
            out['W00'] = out['P00'] / (-integ) if integ else np.nan
        return out

    def cost(self):
        """Rough flop count per RK4 step and memory of the operators."""
        st = self.ops[0]
        fl = 2 * 2 * self.d**3
        mem = 0
        for q in (-1, 0, 1):
            for rowsel, colsel, J in st['jumps'][q]:
                fl += 2 * 2 * J.shape[0] * J.shape[1] * (J.shape[0] + J.shape[1]) / 2
                mem += J.nbytes
            Rs = st['rec'][q][0]
            if st['jumps'][q]:
                fl += len(Rs) * self.na**2 * 8 * self.M**2 * (self.Nx + self.Ny) / 2
        return dict(d=self.d, flops_per_step=4 * fl, n_jumps=sum(len(v) for v in st['jumps'].values()),
                    n_classes=len(self.classes), op_memory_MB=mem / 1e6, frame_mismatch=self.frame_mismatch,
                    class_spread=self.class_spread)


# ------------------------------------------------------------------ schemes
X, Y, ZV = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
SIGP_Z = (S2, 1j * S2, 0.0)          # sigma+ about the fibre axis: E_{+1} = -1
SIGM_Z = (S2, -1j * S2, 0.0)         # sigma- about the fibre axis: E_{-1} = 1


def repump_amp(p, gamma):
    """Rabi frequency (units of nu_r) that gives a scattering rate gamma per sublevel at the probe detuning."""
    return np.sqrt(4 * gamma * (p['D']**2 + p['G']**2 / 4) / p['G'])


def scheme(name, p, zeeman=40.0, d1=0.0, d2=0.0, gamma=0.0, amps=(1.0, 1.0), wc=1.0, te01_norm=True, rep_pol='lin'):
    """87Rb D1 configurations (fields, B direction, excited levels, stages).

    lp11_hf    switched LP11, hyperfine Lambda D1 F'=1, probe pi on F=2 (polarised along B, B transverse
               along x, no field needed: zeeman = 0), coupling linear perpendicular on F=1, repumpers on
               both legs with the complementary polarisation (rate gamma): the best one-axis scheme of
               notes/zeeman_results.md. Stage 0: lobes along x, stage 1: along y. d1 = two-photon detuning.
    lp11_zl    switched LP11, Zeeman Lambda inside F=2 after Qiu et al. (B along the fibre): probe sigma-
               (circular LP11) |2,2> -> |2',1>, coupling sigma+ (LP01) from |2,0>, F=1 repumper sigma+ (gamma).
    te01_bz    TE01 probe, B along the fibre, coupling sigma+ on F=1, F'=2, target |2,1>; probe frequency 1
               for |2,1> -> |1,1> (Raman shift -Z), frequency 2 for |2,1> -> |1,-1> (0); d1, d2 their
               two-photon detunings; optional F=2 repumper (lin, gamma).
    te01_bx    TE01 probe, B transverse (x), coupling pi (x-polarised LP01) on F=1, F'=2, target |2,1>:
               pi channel to |1,1> (shift -Z), sigma- channel to |1,0> (-Z/2).
    rep_pol (TE01 schemes): polarisation of the F=2 repumper, 'lin' (default: x for B along z, y for B
    along x), 'sm' / 'sp' (sigma-/sigma+ about the fibre axis, B along z: drains |2,2> / |2,-2>), 'pi'
    (polarised along B, B along x).
    Saturation s refers to the LP11 power: a TE01 component carries 1/sqrt 2 of the amplitude when
    te01_norm (equal total power); amps scale the two probe frequencies."""
    Wp, Wc = p['Wp'], p['Wc'] * wc
    n2 = S2 if te01_norm else 1.0
    if name == 'lp11_hf':
        f = [field(2, [(X, 'x', 1.0)], Wp, d1, stage=0, name='probe x-lobes'),
             field(2, [(X, 'y', 1.0)], Wp, d1, stage=1, name='probe y-lobes'),
             field(1, [(Y, 'lp01', 1.0)], Wc, 0.0, name='coupling')]
        if gamma:
            f += [field(2, [(Y, 'lp01', 1.0)], repump_amp(p, gamma), 0.0, incoherent=True, name='repump F=2'),
                  field(1, [(X, 'lp01', 1.0)], repump_amp(p, gamma), 0.0, incoherent=True, name='repump F=1')]
        return dict(fields=f, bdir='x', fes=(1,), zeeman=0.0, stages=2, target=GROUND.index((2, 0)))
    if name == 'lp11_zl':
        # E(2,m) = m Z / 2: coupling resonant from |2,0> when delta_c = E(2,2) - E(2,0) + delta_p - d1
        f = [field(2, [(SIGM_Z, 'x', 1.0)], Wp, 0.0, stage=0, name='probe x-lobes'),
             field(2, [(SIGM_Z, 'y', 1.0)], Wp, 0.0, stage=1, name='probe y-lobes'),
             field(2, [(SIGP_Z, 'lp01', 1.0)], Wc, zeeman - d1, name='coupling')]
        if gamma:
            f += [field(1, [(SIGP_Z, 'lp01', 1.0)], repump_amp(p, gamma), 0.0, incoherent=True, name='repump F=1')]
        return dict(fields=f, bdir='z', fes=(2,), zeeman=zeeman, stages=2, target=GROUND.index((2, 2)))
    te = [(X, 'y', -n2), (Y, 'x', n2)]                          # TE01: E = (-f(y), f(x))
    if name == 'te01_bz':
        E = lambda F, m: G_F[F] * m * zeeman
        f = [field(2, te, Wp * amps[0], E(1, 1) - E(2, 1) + d1, name='probe w1'),
             field(2, te, Wp * amps[1], E(1, -1) - E(2, 1) + d2, name='probe w2'),
             field(1, [(SIGP_Z, 'lp01', 1.0)], Wc, 0.0, name='coupling')]
        if gamma:
            pol = {'lin': X, 'sm': SIGM_Z, 'sp': SIGP_Z}[rep_pol]
            f += [field(2, [(pol, 'lp01', 1.0)], repump_amp(p, gamma), 0.0, incoherent=True, name='repump F=2')]
        return dict(fields=f, bdir='z', fes=(2,), zeeman=zeeman, stages=1, target=GROUND.index((2, 1)))
    if name == 'te01_bx':
        E = lambda F, m: G_F[F] * m * zeeman
        f = [field(2, te, Wp * amps[0], E(1, 1) - E(2, 1) + d1, name='probe w1 (pi)'),
             field(2, te, Wp * amps[1], E(1, 0) - E(2, 1) + d2, name='probe w2 (sigma-)'),
             field(1, [(X, 'lp01', 1.0)], Wc, 0.0, name='coupling')]
        if gamma:
            pol = {'lin': Y, 'pi': X}[rep_pol]
            f += [field(2, [(pol, 'lp01', 1.0)], repump_amp(p, gamma), 0.0, incoherent=True, name='repump F=2')]
        return dict(fields=f, bdir='x', fes=(2,), zeeman=zeeman, stages=1, target=GROUND.index((2, 1)))
    raise ValueError(name)


def build(name, p, Nx, Ny, tol=2.0, nodes=(4, 8), **kw):
    sc = scheme(name, p, **kw)
    m = Model(sc['fields'], p, Nx, Ny, bdir=sc['bdir'], zeeman=sc['zeeman'], fes=sc['fes'], tol=tol,
              nodes=nodes, stages=sc['stages'])
    m.target = sc['target']
    return m
