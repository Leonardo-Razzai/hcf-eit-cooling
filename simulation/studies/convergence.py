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

  evolution  banded_evolve vs the full BDF2 evolve, case A, s = 0.2, thermal <n> = 2, N = 20
                                                 -> evolution_check.dat  (~5 min)
  lines_elim  eliminated vs full ME for the D1/D2 configurations (all excited levels),
         Delta = 1, 3, 12 Gamma; three-level bp = 1/6, 5/6
                                                 -> lines_elim.dat  (~15 min)
  lines_nconv  N-convergence versus bp (5/6, 1/6), s = 0.1, 1, 5
                                                 -> lines_nconv.dat (~10 min)
  lines_nconv_delta  N-convergence versus Delta/Gamma, bp = 1/2, 5/6, s = 1
                                                 -> lines_nconv_delta.dat
  lines_nconv_d2  N-convergence of D2 F' = 2, all excited levels, Delta = 8, 12 Gamma
                                                 -> lines_nconv_d2.dat
  lines_d2scan  W, n_ss versus two-photon detuning around the stored optimum (N = 80)
                                                 -> lines_d2scan.dat (~15 min)
  rate0  W against the ground-state rate W0: case A versus s, N, bp; D1/D2 points where W
         runs away or has a pole                  -> rate_p0.dat (~30 min)

Usage: python convergence.py [elim] [band] [nconv] [tail] [matrices] [weak] [rates] [tailshape]
       [evolution] [lines_elim] [lines_nconv] [lines_nconv_delta] [lines_nconv_d2] [lines_d2scan] [rate0]
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


def evolution():
    """banded_evolve against the full BDF2 evolve (case A, s = 0.2, thermal <n> = 2,
    N = 20, same time grid): band of the study and untruncated band, step sizes."""
    fc = CASES['A']
    s, N = 0.2, 20
    p = E.params(fc, s=s, N=N)
    tmax = 375.0
    th = E.thermal(2.0, N)
    L, d, b = E.liouvillian(p)
    rho0 = np.zeros((d, d), complex)
    rho0[np.arange(N), np.arange(N)] = th             # g_p is the first block
    t, full = E.evolve(L, rho0, tmax, E.number_op(b), h=0.05, n_out=60)
    cols, data = ['t', 'full'], [t, full]
    for B, h in ((N, 0.25), (4, 0.25), (4, 0.125), (4, 0.0625)):
        tb, nb = E.banded_evolve(p, th, tmax, B=B, h=h, n_out=60)
        assert np.allclose(tb, t)
        cols.append('band%d_h%g' % (B, h))
        data.append(nb)
        dev = abs(nb / full - 1)
        log('evolution band', B, 'h', h, 'dev from full: max %.2e at t = %.1f, median %.2e'
            % (dev.max(), t[dev.argmax()], np.median(dev)))
    save('evolution_check.dat', cols, np.column_stack(data))


def lines_elim():
    """Eliminated against full ME for the D1/D2 configurations, N = 20, s = 1:
    multi-level at the optimised two-photon detuning of the stored stage-2 scan
    (Delta = 1, 3, 12 Gamma), and three-level with bp = 1/6, 5/6."""
    from lines import configs
    cool = [c for c in configs() if c['bp'] >= 0.5]
    dl = np.genfromtxt(RESULTS / 'lines_delta_A.dat', names=True)
    rows = []
    for ic, cfg in enumerate(cool):
        fc = E.paper_cases(cfg['line'])['A']
        lv = E.excited_manifold(cfg['line'], cfg['Fe'], cfg['gp'])
        for DG in (1.0, 3.0, 12.0):
            m = (dl['config'] == ic) & (dl['DeltaG'] == DG)
            p = E.params(fc, s=1.0, N=20, Delta_G=DG, d2=float(dl['d2opt'][m][0]), levels=lv)
            a, b_ = E.solve_me(p), E.solve_me(p, eliminate=True)
            rows.append([int(cfg['line'][1]), cfg['Fe'], cfg['gp'], cfg['bp'], DG, len(lv),
                         a['W'], b_['W'], a['n'], b_['n']])
            log('lines_elim', cfg['line'], "F'=%d gp=%d D=%g K=%d" % (cfg['Fe'], cfg['gp'], DG, len(lv)),
                'W %.5f elim %+.2f%%  n %.4f elim %+.2f%%'
                % (a['W'], 100 * (b_['W'] / a['W'] - 1), a['n'], 100 * (b_['n'] / a['n'] - 1)))
    fc = CASES['A']
    for bp in (1 / 6, 5 / 6):
        p = E.params(fc, s=1.0, N=20, bp=bp)
        a, b_ = E.solve_me(p), E.solve_me(p, eliminate=True)
        rows.append([1, 0, 0, bp, 3.0, 1, a['W'], b_['W'], a['n'], b_['n']])
        log('lines_elim 3-level bp %.3f' % bp, 'W %+.2f%% n %+.2f%%'
            % (100 * (b_['W'] / a['W'] - 1), 100 * (b_['n'] / a['n'] - 1)))
    save('lines_elim.dat', ['line', 'Fe', 'gp', 'bp', 'DeltaG', 'K', 'Wfull', 'Welim', 'nfull', 'nelim'], rows)


def lines_nconv():
    """N-convergence of the steady state versus the branching ratio (three-level,
    case A), band 4: does the rule tuned at bp = 1/2 hold for bp = 5/6 and 1/6?"""
    fc = CASES['A']
    rows = []
    for bp in (5 / 6, 1 / 6):
        for s in (0.1, 1.0, 5.0):
            for N in (60, 120, 250):
                r = E.solve_banded(E.params(fc, s=s, N=N, bp=bp), B=4)
                rows.append([bp, s, N, r['n'], r['W']])
                log('lines_nconv bp %.3f s %g N %d' % (bp, s, N), 'n %.5f W %.5f' % (r['n'], r['W']))
    save('lines_nconv.dat', ['bp', 's', 'N', 'n', 'W'], rows)


def lines_nconv_delta():
    """N-convergence at s = 1 away from Delta = 3 Gamma (three-level, case A, bp = 1/2, 5/6):
    is the rule tuned at 3 Gamma valid for the Delta scan of lines.py? Band 6."""
    fc = CASES['A']
    rows = []
    for bp in (1 / 2, 5 / 6):
        for DG in (1.0, 2.0, 5.0, 8.0, 12.0):
            for N in (60, 120, 250) + ((400,) if (bp == 1 / 2 and DG > 4) else ()):
                r = E.solve_banded(E.params(fc, s=1.0, N=N, bp=bp, Delta_G=DG), B=6 if N < 400 else 4)
                rows.append([bp, DG, N, r['n'], r['W']])
                log('lines_nconv_delta bp %.3f D %g N %d' % (bp, DG, N), 'n %.5f W %.5f' % (r['n'], r['W']))
    save('lines_nconv_delta.dat', ['bp', 'DeltaG', 'N', 'n', 'W'], rows)


def lines_nconv_d2():
    """N-convergence of D2 F' = 2 with all excited levels, large Delta, s = 1, at the
    stored optimum of lines_delta_A; D1 F' = 2 as reference. Band 4."""
    dl = np.genfromtxt(RESULTS / 'lines_delta_A.dat', names=True)
    rows = []
    for line, Fe, gp, DG, Ns in (('D1', 2, 2, 12.0, (80, 160, 250, 400)),
                                 ('D2', 2, 1, 12.0, (80, 160, 250, 400)),
                                 ('D2', 2, 2, 8.0, (80, 160, 250, 400)),
                                 ('D2', 2, 2, 12.0, (80, 160, 250, 400, 550))):
        fc = E.paper_cases(line)['A']
        lv = E.excited_manifold(line, Fe, gp)
        m = (dl['line'] == int(line[1])) & (dl['Fe'] == Fe) & (dl['gp'] == gp) & (dl['DeltaG'] == DG)
        d2 = float(dl['d2opt'][m][0])
        popt = E.params(fc, s=1.0, Delta_G=DG, d2=d2, levels=lv)
        beff = float(np.dot(E.excited_weights(popt, 1.0 + d2), [l[3] for l in popt['levels']]))
        for N in Ns:
            r = E.steady(fc, 1.0, cut=(N, 4), Delta_G=DG, d2=d2, levels=lv)
            W = r['W'] * fc['om']
            rows.append([int(line[1]), Fe, gp, DG, d2, beff, N, r['n'], W])
            log('lines_nconv_d2', line, "F'=%d gp=%d D %g beff %.3f N %d" % (Fe, gp, DG, beff, N),
                'n %.5f W %.1f' % (r['n'], W))
    save('lines_nconv_d2.dat', ['line', 'Fe', 'gp', 'DeltaG', 'd2', 'beff', 'N', 'n', 'W'], rows)


def lines_d2scan():
    """W and n_ss versus the two-photon detuning around the stored optimum (s = 1, all
    excited levels, sign +1, N = 80, band 4: shapes only): is the search of lines.py on
    the peak, and is the selection n_ss < 0.2 limiting?"""
    dl = np.genfromtxt(RESULTS / 'lines_delta_A.dat', names=True)
    rows = []
    for line, Fe, gp, DG in (('D2', 2, 1, 3.0), ('D2', 2, 2, 3.0), ('D2', 2, 1, 8.0), ('D2', 2, 2, 8.0),
                             ('D2', 2, 1, 12.0), ('D2', 2, 2, 12.0), ('D1', 1, 2, 12.0),
                             ('D1', 2, 1, 12.0), ('D1', 1, 2, 3.0)):
        fc = E.paper_cases(line)['A']
        lv = E.excited_manifold(line, Fe, gp)
        m = (dl['line'] == int(line[1])) & (dl['Fe'] == Fe) & (dl['gp'] == gp) & (dl['DeltaG'] == DG)
        d0 = float(dl['d2opt'][m][0])
        w = 0.3 * min(1.0, 3.0 / DG)
        for d2 in d0 + np.linspace(-w, w, 17):
            r = E.steady(fc, 1.0, cut=(80, 4), Delta_G=DG, d2=d2, levels=lv)
            rows.append([int(line[1]), Fe, gp, DG, d0, d2, r['W'] * fc['om'], r['n']])
        log('lines_d2scan', line, Fe, gp, DG)
    save('lines_d2scan.dat', ['line', 'Fe', 'gp', 'DeltaG', 'd2stored', 'd2', 'W', 'n'], rows)


def rate0():
    """W against the ground-state rate W0 (solve_banded): case A at the paper settings
    versus s and N, bp = 5/6, and the D1/D2 points where W runs away or has a pole."""
    rows = []

    def row(lab, fc, r, s, DG, bp, d2, N, B):
        W, W0 = r['W'] * fc['om'], r['W0'] * fc['om']
        rows.append([lab, s, DG, bp, d2, N, B, r['n'], r['P0'], W, W0])
        log('rate0 set %d s %g D %g bp %.3f d2 %.3f N %d' % (lab, s, DG, bp, d2, N),
            'n %.4f P0 %.4f W %.0f W0 %.0f' % (r['n'], r['P0'], W, W0))
    fc = CASES['A']
    for s, Ns in ((0.01, (120, 400)), (0.1, (120, 400)), (0.3, (120, 250)), (1.0, (20, 40, 120)),
                  (5.0, (40, 120))):
        for N in Ns:
            B = 4 if N > 160 else 6
            row(0, fc, E.solve_banded(E.params(fc, s=s, N=N), B=B), s, 3.0, 0.5, 0.0, N, B)
    row(1, fc, E.solve_banded(E.params(fc, s=1.0, N=120, bp=5 / 6), B=6), 1.0, 3.0, 5 / 6, 0.0, 120, 6)
    row(1, fc, E.solve_banded(E.params(fc, s=0.1, N=120, bp=5 / 6), B=6), 0.1, 3.0, 5 / 6, 0.0, 120, 6)
    row(1, fc, E.solve_banded(E.params(fc, s=5.0, N=120, bp=5 / 6), B=6), 5.0, 3.0, 5 / 6, 0.0, 120, 6)
    dl = np.genfromtxt(RESULTS / 'lines_delta_A.dat', names=True)
    for line, Fe, gp, DG, Ns in (('D2', 2, 2, 12.0, (80, 160, 250, 400)), ('D2', 2, 1, 3.0, (120,)),
                                 ('D2', 2, 2, 3.0, (120,)), ('D1', 2, 2, 3.0, (120,)),
                                 ('D1', 1, 2, 3.0, (120,))):
        fl = E.paper_cases(line)['A']
        lv = E.excited_manifold(line, Fe, gp)
        m = (dl['line'] == int(line[1])) & (dl['Fe'] == Fe) & (dl['gp'] == gp) & (dl['DeltaG'] == DG)
        d2 = float(dl['d2opt'][m][0])
        for N in Ns:
            B = 4 if N > 160 else 6
            r = E.solve_banded(E.params(fl, s=1.0, N=N, Delta_G=DG, d2=d2, levels=lv), B=B)
            row(2 if line == 'D1' else 3, fl, r, 1.0, DG, lv[0]['bp'], d2, N, B)
    # across the pole of W versus d2: D2 F'=2, gp = 2, 3 Gamma (lines_d2scan.dat)
    fl = E.paper_cases('D2')['A']
    lv = E.excited_manifold('D2', 2, 2)
    for d2 in (-0.05, 0.0, 0.03, 0.06, 0.09, 0.13, 0.2, 0.3):
        r = E.solve_banded(E.params(fl, s=1.0, N=80, d2=d2, levels=lv), B=4)
        row(4, fl, r, 1.0, 3.0, 0.5, d2, 80, 4)
    save('rate_p0.dat', ['set', 's', 'DeltaG', 'bp', 'd2', 'N', 'band', 'n', 'P0', 'W', 'W0'], rows)


STUDIES = {'elim': elim, 'band': band, 'nconv': nconv, 'tail': tail, 'matrices': matrices,
           'weak': weak, 'rates': rates, 'tailshape': tailshape, 'evolution': evolution,
           'lines_elim': lines_elim, 'lines_nconv': lines_nconv,
           'lines_nconv_delta': lines_nconv_delta, 'lines_nconv_d2': lines_nconv_d2,
           'lines_d2scan': lines_d2scan, 'rate0': rate0}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
