"""Both radial axes with Zeeman structure, by an effective rate chain (no 2D Zeeman master equation).

The analytic chain of two_axes.py is calibrated per scheme on the 1D Zeeman master equation: a rate scale
and a blue-sideband scale fitted to P0 and W0, with the photons per cycle, the emission pattern and the
background scattering at n = 0 read from short-time response. Reading rates level by level was tried first and
failed the closed-Lambda benchmark (decay of a Fock state is as fast as the internal dynamics).
Tests: P(n) of the master equation, not used in the fit; the exact 2D result of the closed Lambda.

Usage: python two_axes_zeeman.py [response] [chain]
       response -> results/two_axes_zeeman_response.dat (~1 h);  chain -> results/two_axes_zeeman.dat
"""
import pickle
import sys

import numpy as np

from _common import E, RESULTS, save, log
from eitsim import zeeman as Z
from eitsim.operators import oscillator
from zeeman_scan import POL, repumpers

N0S = np.array([0, 1, 2, 3, 5, 8, 12, 20])
AA = E.paper_cases()['A']['eta_sp']**-2


def schemes():
    """name, nu_wrec, kind, settings, emission axis of the other radial direction."""
    out = []
    for nw in (AA, 5.0):
        lab = 'A' if nw == AA else '5'
        out.append(('closed5/6_' + lab, nw, 'closed', dict(bp=5 / 6, d2=0.0, s=1.0, rep=0.0), 'perp'))
        out.append(('closed1/2_' + lab, nw, 'closed', dict(bp=1 / 2, d2=0.0, s=1.0, rep=0.0), 'perp'))
    out += [('hf97_A', AA, 'hf', dict(d2=0.60, s=1.0, rep=0.10), 'par'), ('hf99_A', AA, 'hf', dict(d2=0.75, s=0.3, rep=0.03), 'par'),
            ('hf97_5', 5.0, 'hf', dict(d2=0.60, s=1.0, rep=0.10), 'par'),
            ('z97_A', AA, 'z', dict(d2=0.80, s=1.0, rep=1.0), 'perp'), ('z99_A', AA, 'z', dict(d2=0.80, s=0.3, rep=1.0), 'perp'),
            ('z97_5', 5.0, 'z', dict(d2=0.40, s=1.0, rep=1.0), 'perp'), ('z99_5', 5.0, 'z', dict(d2=0.40, s=0.3, rep=1.0), 'perp')]
    return out


def build(kind, nw, c, N):
    fc = E.reduced_case(nw)
    p = Z.setup(fc, c['s'], 3.0, c['d2'], N)
    if kind == 'closed':
        p = dict(E.params(fc, s=c['s'], N=N, bp=c['bp'], d2=c['d2']), om=fc['om'])
        return Z.structure_three_level(c['bp']), p, ()
    if kind == 'hf':
        st = Z.structure_lambda('D1', 1, 2, POL['lin'], POL['pi'])
        return st, p, repumpers(st, 'lin', 'pi', p, c['rep'])
    st = Z.structure_zeeman_lambda()
    na, nE = len(st['ground']), len(st['excited'])
    g = Z.field_matrix('D1', 1, st['excited'], {1: 1.0}) / np.sqrt(Z._norm('D1', 1, 2))
    gm = np.zeros((na, nE), complex)
    gm[5:] = g
    return st, p, [dict(gm=gm, W=np.sqrt(4 * c['rep'] * (p['D']**2 + p['G']**2 / 4) / p['G']), detuning=p['D'])]


def response(only=None):
    out = {}
    for name, nw, kind, c, ax in schemes():
        if only and not name.startswith(only):
            continue
        st, p, rp = build(kind, nw, c, 30)
        B = 2 if kind != 'closed' else 4
        r = Z.fock_response(st, p, N0S, B=B, repump=rp)
        s1 = Z.solve(st, p, B=B, repump=rp)
        out[name] = dict(nw=nw, kind=kind, c=c, ax=ax, W0=s1['W0'], W=s1['W'], **{k: r[k] for k in ('drift', 'rate', 'P', 'n', 'P0', 'ss_rate')})
        log('response', name, 'P0 %.4f n %.4f | G(0) %.4f G(1) %.4f G(20) %.4f | q fractions at n=1: %s | drift' % (r['P0'], r['n'], r['rate'][0].sum(), r['rate'][1].sum(), r['rate'][-1].sum(), np.round(r['rate'][1] / r['rate'][1].sum(), 2)), np.round(r['drift'], 3))
        with open(RESULTS / 'two_axes_zeeman_response.pkl', 'wb') as f:
            pickle.dump(out, f)


def one_photon(nw, frac, axis, N):
    b, w, V = oscillator(N)
    eta = nw**-0.5
    P = np.zeros((N, N))
    for q, fq in zip((-1, 0, 1), frac):
        if fq > 0:
            u, wg = Z.pattern_nodes(q, axis, 8)
            for uk, wk in zip(u, wg):
                P += fq * wk * np.abs((V * np.exp(1j * eta * uk * w)) @ V.conj().T)**2
    return P / P.sum(axis=0, keepdims=True)


def alpha(frac, axis):
    a = {('perp', 0): 0.4, ('perp', 1): 0.3, ('par', 0): 0.2, ('par', 1): 0.4}
    return sum(fq * a[(axis, abs(q))] for q, fq in zip((-1, 0, 1), frac))


C_H = 1.5


def kappa_of(r):
    """Photons per removed phonon, from the 1D response at low n, where the two rates are read
    in the same window; net of recoil heating."""
    tot = r['rate'].sum(axis=1)
    frac = r['rate'][1:].sum(axis=0) / r['rate'][1:].sum()
    a = alpha(frac, 'perp')
    eta2 = 1.0 / r['nw']
    k = [tot[i] / (-r['drift'][i] + tot[i] * a * eta2) for i in (1, 2, 3)]
    return np.array(k), frac, a


def cycle_kernel(P, kappa):
    b = 1.0 / kappa
    return b * P @ np.linalg.inv(np.eye(P.shape[0]) - (1 - b) * P)


def model(r, kappa, A, bf, N):
    """Analytic chain of two_axes.py with a rate scale A, a blue-sideband scale bf, kappa photons
    per cycle, and the measured background photon rate. Returns the matrices."""
    nw, s = r['nw'], r['c']['s']
    Gn = 1.0 / 3.0
    frac = r['rate'][1:].sum(axis=0) / r['rate'][1:].sum()
    a_s = alpha(frac, 'perp')
    h = kappa * a_s / nw
    c = 1.125 * kappa
    n = np.arange(N, dtype=float)
    R = A * s * n * Gn / 2 / (1 + c * s * n)
    B = bf * (Gn / 4)**2 * (n + 1) * s * Gn / 2
    B[-1] = 0.0
    K = np.diag(R[1:], 1) + np.diag(B[:-1], -1)
    d = R + B
    nc = 2 / (Gn**2 * s)
    Hr = C_H * A * max(1 - h, 0.0) * Gn / (2 * c) * np.maximum(n - 0.8 * nc, 0) / nc
    Hm = np.diag(Hr[1:], 1) - np.diag(Hr)
    Ps, Po = one_photon(nw, frac, 'perp', N), one_photon(nw, frac, r['ax'], N)
    G0 = max(r['rate'][0].sum() - kappa * B[0], 0.0)
    return dict(K=K, d=d, Hm=Hm, Ps=Ps, Po=Po, Ts=cycle_kernel(Ps, kappa), To=cycle_kernel(Po, kappa), G0=G0, h=h, a_s=a_s)


def steady_1d(m, N):
    Q = m['Ts'] @ m['K'] - np.diag(m['d']) + m['Hm'] + m['G0'] * (m['Ps'] - np.eye(N))
    A = Q.copy()
    A[0] = 1.0
    rhs = np.zeros(N)
    rhs[0] = 1.0
    p = np.linalg.solve(A, rhs)
    src = -p.copy()
    src[1] += 1.0
    src[0] = 0.0
    y = np.linalg.solve(A, src)
    return p, p[0] / y[0]


def calibrate(r, kappa, N=60):
    from scipy.optimize import least_squares

    def err(v):
        m = model(r, kappa, np.exp(v[0]), np.exp(v[1]), N)
        p, W0 = steady_1d(m, N)
        return [np.log((1 - p[0]) / (1 - r['P0'])), np.log(max(W0, 1e-9) / r['W0'])]
    s = least_squares(err, [0.0, 0.0], diff_step=0.02)
    return np.exp(s.x[0]), np.exp(s.x[1]), np.abs(s.fun).max()


def evolve2(m, N, tmax=20000.0):
    p, _ = steady_1d(m, N)
    P = np.outer(p, p)
    MK = m['Ts'] @ m['K']
    To, Ps, Po, Hm, d, G0 = m['To'], m['Ps'], m['Po'], m['Hm'], m['d'], m['G0']
    nn = np.arange(N, dtype=float)
    dt = 0.9 / (d.max() + abs(Hm).max() + G0)
    t, last, k = 0.0, None, 0
    every = max(1, int(50.0 / dt))
    while t < tmax:
        x = MK @ P @ To.T - d[:, None] * P + Hm @ P + G0 * (Ps @ P @ Po.T - P)
        y = To @ P @ MK.T - P * d[None, :] + P @ Hm.T + G0 * (Po @ P @ Ps.T - P)
        P = P + dt * 0.5 * (x + y)
        t += dt
        k += 1
        if k % every == 0:
            mx = P.sum(axis=1)
            v = (nn @ mx, P[0, 0])
            if last and abs(v[0] - last[0]) < 1e-4 * v[0] and abs(v[1] - last[1]) < 1e-5:
                break
            last = v
    return P, t


def chain_all():
    with open(RESULTS / 'two_axes_zeeman_response.pkl', 'rb') as f:
        data = pickle.load(f)
    KNOWN = {'closed5/6': 1.2, 'closed1/2': 2.0}
    rows = []
    for k, (name, r) in enumerate(data.items()):
        kr, frac, a = kappa_of(r)
        kappa = KNOWN.get(name.split('_')[0], float(np.median(kr)))
        A, bf, fit = calibrate(r, kappa)
        N = 120
        m = model(r, kappa, A, bf, N)
        p, W0 = steady_1d(m, N)
        # independent test: P(n) of the 1D ME, not used in the fit (n = 1..12)
        pme = r['P']
        ns = np.arange(1, 13)
        rms = np.sqrt(np.mean(np.log(p[1:13] / pme[1:13])**2))
        P, t = evolve2(m, N)
        mx = P.sum(axis=1)
        nn = np.arange(N)
        log('chain', name, 'kappa read %s used %.2f A %.2f bf %.2f fit %.1e G0 %.4f h %.2f | P(n) rms log %.2f, n: ME %.3f chain %.3f | 2D P00 %.4f P0x %.4f nx %.3f t %.0f edge %.0e'
            % (np.round(kr, 2), kappa, A, bf, fit, m['G0'], m['h'], rms, r['n'], nn @ p, P[0, 0], mx[0], nn @ mx, t, mx[-10:].sum()))
        rows.append([k, r['nw'], kappa, A, bf, m['G0'], m['h'], r['P0'], r['W0'], rms, r['n'], nn @ p, P[0, 0], mx[0], nn @ mx, t])
    save('two_axes_zeeman.dat', ['scheme', 'nu_wrec', 'kappa', 'A', 'bf', 'G0', 'h', 'P0_1D', 'W0_1D', 'Pn_rms_log', 'n_me', 'n_chain', 'P00', 'P0x', 'nx', 't_end'], rows)
    log('schemes', list(data))


if __name__ == '__main__':
    for name in sys.argv[1:]:
        {'response': response, 'bench': lambda: response('closed'), 'chain': chain_all}[name]()
