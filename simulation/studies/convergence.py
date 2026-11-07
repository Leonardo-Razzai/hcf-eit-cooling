"""Checks of the approximations behind the converged case-A master equation.

All output goes to results/convergence/ (tables) and is plotted by
plotting/convergence_plots.py into the same folder. One process, sequential.

  elim   full (three-level) vs eliminated master equation, N = 20: cases A, B,
         s = 0.03..5, Delta = Gamma and 3 Gamma, plus an LP01 admixture and a
         two-photon error                       -> elim_check.dat   (~2 min)
  band   exact eliminated vs banded (|n - m| <= B) at N = 40, case A, s = 0.1, 1, 5;
         decay of the steady-state coherences with |n - m|
                                                 -> band_check.dat, coherences.dat (~4 min)
  nconv  n_ss and W versus Fock truncation N: full ME (N <= 40; 60 at s = 1) and
         banded (B = 6, N = 20..160), case A at s = 0.1, 1, 5 and case B at s = 1
                                                 -> nconv.dat       (~15 min)
  tail   Fock populations P(n) at N = 160 (banded)
                                                 -> tail.dat        (~2 min)
  matrices  the elements discarded by the two cutoffs: |rho_(a,n),(b,m)| of the exact
         eliminated steady state (case A, s = 1, N = 40; band |n - m| <= B), and the
         untruncated operators f(x) (recoil exp(i eta_sp u x) for the largest |u|,
         probe and coupling profiles; function_of_x drops elements < 1e-11 x max)
                                                 -> rho_A_s1.dat, op_recoil.dat,
                                                    op_probe.dat, op_coupling.dat (~1 min)

  weak   weak-probe N-scan, case A, s = 0.01, 0.1, 0.3, band 4, N = 60..300 (decides the
         provisional numerics.toml rule)      -> nconv_weak.dat   (~20-30 min)

  rates  cooling drift d<n>/dt and energy diffusion dVar(n)/dt from Fock states n0 = 1..200,
         case A, s = 0.01, 0.1, 1 (banded, N = 260, B = 4): why the tail depends on s
                                                 -> rates_vs_n.dat  (~10 min)
  tailshape  P(n) for case A, s = 0.1 at N = 300 and 400 (band 4)
                                                 -> tail_s0p1.dat   (~5 min)

Usage: python convergence.py [elim] [band] [nconv] [tail] [matrices] [weak] [rates] [tailshape]
"""
import sys

import numpy as np

from _common import E, CASES, RESULTS, log

OUT = RESULTS / 'convergence'
BAND = 6


def save(name, cols, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT / name, np.array(rows), header=' '.join(cols), comments='')


def elim():
    rows = []
    cases = [(key, DG, s, 0.0, 0.0) for key in 'AB' for DG in (1.0, 3.0)
             for s in (0.03, 0.1, 0.3, 1.0, 2.0, 5.0)]
    cases += [('A', 3.0, 1.0, 0.02, 0.0), ('A', 3.0, 0.2, 0.0, 0.2)]
    for key, DG, s, eps, d2 in cases:
        p = E.params(CASES[key], s=s, N=20, Delta_G=DG, eps=eps, d2=d2)
        f, e = E.solve_me(p), E.solve_me(p, eliminate=True)
        rows.append([0 if key == 'A' else 1, DG, s, eps, d2, f['n'], f['W'], e['n'], e['W']])
        log('elim', key, 'Delta %g G s %g eps %g d2 %g' % (DG, s, eps, d2),
            'dn %+.1e dW %+.1e' % (e['n'] / f['n'] - 1, e['W'] / f['W'] - 1))
    save('elim_check.dat', ['case', 'DeltaG', 's', 'eps', 'd2', 'n_full', 'W_full',
                            'n_elim', 'W_elim'], rows)


def band():
    N = 40
    rows, coh = [], []
    for s in (0.1, 1.0, 5.0):
        p = E.params(CASES['A'], s=s, N=N)
        ex = E.solve_me(p, eliminate=True)
        rho = ex['rho']
        k = np.arange(N)
        cmax = [max(np.abs(np.diagonal(rho[i * N:(i + 1) * N, j * N:(j + 1) * N], kk)).max()
                    for i in range(2) for j in range(2) for kk in (kk_, -kk_)) for kk_ in k]
        coh += [[s, kk_, c] for kk_, c in zip(k, cmax)]
        for B in (1, 2, 3, 4, 6, 8):
            r = E.solve_banded(p, B=B)
            rows.append([s, B, ex['n'], ex['W'], r['n'], r['W']])
            log('band', 's %g B %d' % (s, B), 'dn %+.1e dW %+.1e' % (r['n'] / ex['n'] - 1,
                                                                    r['W'] / ex['W'] - 1))
    save('band_check.dat', ['s', 'B', 'n_exact', 'W_exact', 'n_band', 'W_band'], rows)
    save('coherences.dat', ['s', 'k', 'max_abs_rho'], coh)


def nconv():
    rows = []
    for key, s in (('A', 0.1), ('A', 1.0), ('A', 5.0), ('B', 1.0)):
        fc = CASES[key]
        cid = 0 if key == 'A' else 1
        for N in (20, 30, 40) + ((60,) if (key, s) == ('A', 1.0) else ()):
            f = E.solve_me(E.params(fc, s=s, N=N))
            rows.append([cid, s, 0, N, f['n'], f['W'] * fc['om']])
            log('nconv full', key, s, N, '%.5f %.1f' % (f['n'], f['W'] * fc['om']))
        for N in (20, 30, 40, 60, 80, 100, 120, 160):
            r = E.solve_banded(E.params(fc, s=s, N=N), B=BAND)
            rows.append([cid, s, 1, N, r['n'], r['W'] * fc['om']])
            log('nconv band', key, s, N, '%.5f %.1f' % (r['n'], r['W'] * fc['om']))
    save('nconv.dat', ['case', 's', 'method', 'N', 'n', 'W'], rows)


def tail():
    N = 160
    cols, data = ['n'], [np.arange(N)]
    for key, s in (('A', 0.1), ('A', 1.0), ('A', 5.0), ('B', 1.0)):
        r = E.solve_banded(E.params(CASES[key], s=s, N=N), B=BAND)
        cols.append('P_%s_s%s' % (key, ('%g' % s).replace('.', 'p')))
        data.append(np.clip(r['P'], 0, None))
        log('tail', key, s, 'P(n>=20) %.2e' % r['P'][20:].sum())
    save('tail.dat', cols, np.column_stack(data))


def matrices():
    from eitsim.operators import oscillator
    from eitsim.fibre import amp_ratio, U01, U1
    from scipy.special import jv
    N = 40
    fc = CASES['A']
    p = E.params(fc, s=1.0, N=N)
    rho = E.solve_me(p, eliminate=True)['rho']
    save('rho_A_s1.dat', ['col%d' % j for j in range(2 * N)], np.abs(rho))
    b, w, V = oscillator(N)
    full = lambda vals: (V * vals) @ V.conj().T        # f(x) without the 1e-11 cutoff
    u = np.abs(E.NODES).max()
    xa = p['x0a'] * w
    for name, vals in (('op_recoil.dat', np.exp(1j * p['eta_sp'] * u * w)),
                       ('op_probe.dat', amp_ratio(1) * jv(1, U1[0] * xa)),
                       ('op_coupling.dat', jv(0, U01 * xa))):
        save(name, ['col%d' % j for j in range(N)], np.abs(full(vals)))
    log('matrices', 'saved rho (N=%d) and operators; |u|max = %.3f, eta_sp = %.3f'
        % (N, u, p['eta_sp']))


def weak():
    import time
    fc = CASES['A']
    rows = []
    for s in (0.01, 0.1, 0.3):
        for N in (60, 100, 150, 200, 250, 300):
            t = time.time()
            r = E.solve_banded(E.params(fc, s=s, N=N), B=4)
            rows.append([s, N, 4, r['n'], r['W'] * fc['om'], r['P'][-20:].sum()])
            log('weak', 's %g N %d' % (s, N), 'n %.5f W %.1f tail(last 20) %.1e  %.0fs'
                % (r['n'], r['W'] * fc['om'], rows[-1][-1], time.time() - t))
            save('nconv_weak.dat', ['s', 'N', 'band', 'n', 'W', 'P_last20'], rows)
    return rows


RATE_N0 = (1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 120, 160, 200)


def rates():
    import time
    fc = CASES['A']
    rows = []
    for s in (0.01, 0.1, 1.0):
        t = time.time()
        out = E.banded_fock_drift(E.params(fc, s=s, N=260), RATE_N0, B=4)
        for n0, (dr, df) in zip(RATE_N0, out):
            rows.append([s, n0, dr * fc['om'], df * fc['om']])
        log('rates', 's %g' % s, 'drift(n0) /s:', np.round(out[:, 0] * fc['om']), '%.0fs' % (time.time() - t))
        save('rates_vs_n.dat', ['s', 'n0', 'drift', 'diffusion'], rows)


def tailshape():
    fc = CASES['A']
    cols, data = ['n'], [np.arange(400)]
    for N in (300, 400):
        r = E.solve_banded(E.params(fc, s=0.1, N=N), B=4)
        P = np.full(400, np.nan)
        P[:N] = np.clip(r['P'], 0, None)
        cols.append('P_N%d' % N)
        data.append(P)
        log('tailshape', 'N %d n %.5f W %.1f' % (N, r['n'], r['W'] * fc['om']))
        save('tail_s0p1.dat', cols, np.column_stack(data))


STUDIES = {'elim': elim, 'band': band, 'nconv': nconv, 'tail': tail, 'matrices': matrices,
           'weak': weak, 'rates': rates, 'tailshape': tailshape}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
