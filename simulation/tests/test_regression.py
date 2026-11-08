"""Fast regression tests (~10 s): the package must reproduce stored results.

Run from the repository root:  python -m pytest simulation/tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

SIM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIM))
import eitsim as E  # noqa: E402

RESULTS = SIM / 'results'


def load(name):
    return np.genfromtxt(RESULTS / name, names=True)


def test_mode_constants():
    assert abs(E.gfac(1) - 3.49) < 0.01
    assert abs(E.amp_ratio(1) - 1.82) < 0.01
    assert abs(E.ALPHA - 0.397) < 0.002


def test_fibre_cases():
    c = E.paper_cases()
    assert abs(c['A']['nur'] - 7148.8) < 1
    assert abs(c['B']['nur'] - 30637.6) < 1
    assert abs(c['A']['x0'] - 90.19e-9) < 0.1e-9


def test_master_equation_matches_saturation_results():
    d = load('sat_A.dat')                                  # converged solver (numerics.toml)
    fc = E.paper_cases()['A']
    i = int(np.argmin(np.abs(d['s'] - 1.0)))
    r = E.steady(fc, d['s'][i])
    assert r['N'] == d['N'][i] and r['band'] == d['band'][i]
    assert np.isclose(r['n'], d['nME'][i], rtol=1e-8)
    assert np.isclose(r['W'] * fc['om'], d['WME'][i], rtol=1e-8)


def test_simple_model_matches_saturation_results():
    d = load('sat_B.dat')
    fc = E.paper_cases()['B']
    for i in (0, 6, 12):
        sm = E.simple(E.params(fc, s=d['s'][i], N=20), 2)
        assert np.isclose(sm['W'] * fc['om'], d['WSM2'][i], rtol=1e-9)
        assert np.isclose(sm['n'], d['nSM2'][i], rtol=1e-9)


def test_grm_steady_state_matches_check():
    d = load('grmcheck.dat')
    row = d[(d['case'] == 0) & (d['s'] == 0.01)][0]
    p = E.params(E.paper_cases()['A'], s=0.01, N=20)
    ss = E.grm_steady_state(E.grm_matrix(p, N=60, csat=1.0))
    assert np.isclose(np.arange(60) @ ss, row['nGRM'], rtol=1e-6)



def test_temperature_scan_initial_states():
    c = E.paper_cases()
    assert abs(c['A']['U0'] / E.KB * 1e6 - 532.5) < 0.5          # 530 uK lattice
    for key in 'AB':
        d = load(f'temperature_{key}.dat')
        assert len(d) == 20
        assert np.isclose(d['TU0'][0], 1e-3) and np.isclose(d['TU0'][-1], 0.1)
        n0 = 1 / np.expm1(E.HBAR * c[key]['om'] / (d['TU0'] * c[key]['U0']))
        assert np.allclose(d['n0'], n0, rtol=1e-9)
        assert np.all(d['p0end'] > 0.9)
        tr = load(f'temptrace_{key}.dat')
        assert [c for c in tr.dtype.names if c.startswith('n_')] == \
            ['n_0p001', 'n_0p004', 'n_0p01', 'n_0p04', 'n_0p1']
        # t = 0 mean of the truncated thermal state vs the analytic Bose value
        assert np.allclose([tr['n_0p001'][0], tr['n_0p1'][0]], d['n0'][[0, -1]], rtol=1e-3)


def _temperature_module():
    sys.path.insert(0, str(SIM / 'studies'))
    import temperature
    return temperature


def test_first_crossing():
    fc = _temperature_module().first_crossing
    t = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([4.0, 2.0, 1.0, 0.5])
    assert fc(t, y, 1.5) == 1.5                      # linear interpolation
    assert fc(t, y, 5.0) == 0.0                      # already below at t = 0
    assert np.isnan(fc(t, y, 0.1))                   # never reached


def test_temperature_scan_point_case_B():
    T = _temperature_module()
    fc = E.paper_cases()['B']
    d = load('temperature_B.dat')
    j = 10
    N = T.NFOCK['B']
    Q = E.grm_matrix(E.params(fc, s=1.0, N=N), N=N, csat=E.CSAT)
    p0 = E.thermal(T.bose(fc, d['TU0'][j]), N)
    assert np.isclose(-(np.arange(N) @ Q @ p0) * fc['om'], d['R0'][j], rtol=1e-8)
    assert np.isclose(np.arange(N) @ E.grm_steady_state(Q), d['nss'][j], rtol=1e-6)


def test_grm_drift_beyond_calibration():
    d = load('drift_hi.dat')
    c = E.paper_cases()
    for row in d:
        fc = c['A' if row['case'] == 0 else 'B']
        n0 = int(row['n0'])
        NG = 90 if n0 == 20 else 120
        g = (np.arange(NG) @ E.grm_matrix(E.params(fc, s=1.0, N=NG), N=NG, csat=E.CSAT))[n0]
        assert np.isclose(g * fc['om'], row['GRM'], rtol=1e-8)
        assert row['GRM'] / row['ME'] < 1           # GRM is conservative at high n


# ---------------------------------------------------------------- atomic data
STECK_S = {('D1', 2, 2): 1 / 2, ('D1', 2, 1): 1 / 2, ('D1', 1, 2): 5 / 6, ('D1', 1, 1): 1 / 6,
           ('D2', 2, 3): 7 / 10, ('D2', 2, 2): 1 / 4, ('D2', 2, 1): 1 / 20,
           ('D2', 1, 2): 5 / 12, ('D2', 1, 1): 5 / 12, ('D2', 1, 0): 1 / 6}


def test_hyperfine_strengths_match_steck_table8():
    for (line, F, Fe), v in STECK_S.items():
        assert np.isclose(E.strength(line, F, Fe), v, atol=1e-12)
    for line in E.LINES:
        for F in E.GROUND_F:
            assert np.isclose(sum(E.strength(line, F, f) for f in E.excited_levels(line)), 1)
        for Fe in E.excited_levels(line):
            assert np.isclose(sum(E.branching(line, Fe, F) for F in E.GROUND_F), 1)
    assert np.isclose(E.branching('D1', 1, 1), 1 / 6) and np.isclose(E.branching('D1', 1, 2), 5 / 6)


def test_zeeman_couplings_sum_rules():
    for line in E.LINES:
        Je = float(E.LINES[line]['Je'])
        for Fe in E.excited_levels(line):          # Steck Eq. (40): same total decay for every m_F'
            for me in range(-Fe, Fe + 1):
                tot = sum(E.zeeman_coupling(line, F, me + q, Fe, me, q)**2
                          for F in E.GROUND_F for q in (-1, 0, 1) if abs(me + q) <= F)
                assert np.isclose(tot, 2 / (2 * Je + 1))
        for F in E.GROUND_F:                        # Steck Eq. (41): sum over m_F', q gives S_FF'
            for Fe in E.excited_levels(line):
                tot = sum(E.zeeman_coupling(line, F, 0, Fe, -q, q)**2 for q in (-1, 0, 1) if abs(q) <= Fe)
                assert np.isclose(tot, E.strength(line, F, Fe))
    # D2 cycling transition |2,2> -> |3,3>: d / <J||er||J'> = 2.98931 / 4.22752 (Steck Table 7)
    assert np.isclose(abs(E.zeeman_coupling('D2', 2, 2, 3, 3, -1)), 2.98931 / 4.22752, rtol=1e-5)


def test_lambda_configurations():
    for line, Fe in (('D2', 0), ('D2', 3)):
        for gp in (1, 2):
            with pytest.raises(ValueError):
                E.lambda_config(line, Fe, gp)
    c = E.lambda_config('D1', 1, 2)
    assert np.isclose(c['bp'], 5 / 6) and c['gc'] == 1 and np.isclose(c['Sp'], 1 / 2)


def test_line_dependent_fibre_case():
    a, d1, d2 = 15e-6, E.fibre_case(15e-6), E.fibre_case(15e-6, line='D2')
    assert d1 == E.fibre_case(a, line='D1')                    # D1 default unchanged
    assert np.isclose(d2['eta_sp'] / d1['eta_sp'], E.LINES['D1']['lam'] / E.LINES['D2']['lam'])
    p = E.params(d1, s=1.0)                                    # paper: 8 nW probe in case A
    assert np.isclose(E.lp01_power(d1, p['Wp'] * d1['om']) * 1e9, 8.03, atol=0.01)


def test_lines_study():
    sat = load('sat_B.dat')
    d = load('lines_B.dat')
    base = d[(d['line'] == 1) & (d['Fe'] == 2)]                 # paper configuration
    i = int(np.argmin(abs(sat['s'] - 1)))
    assert np.allclose(base['W1'], sat['WME'][i], rtol=1e-8)
    assert np.allclose(base['n1'], sat['nME'][i], rtol=1e-8)
    row = d[(d['line'] == 2) & (d['Fe'] == 1) & (d['gp'] == 1)][0]
    fc = E.paper_cases('D2')['B']
    r = E.solve_me(E.params(fc, s=1.0, N=20, bp=row['bp']))
    assert np.isclose(r['W'] * fc['om'], row['W1'], rtol=1e-6)
    assert np.isclose(r['n'], row['n1'], rtol=1e-6)


# ---------------------------------------------------------------- multi-level (stage 2)
def test_excited_manifold_amplitudes():
    for line, Fe, gp in (('D1', 2, 1), ('D1', 1, 2), ('D2', 1, 1), ('D2', 2, 2)):
        cfg = E.lambda_config(line, Fe, gp)
        lv = E.excited_manifold(line, Fe, gp)
        assert lv[0]['F'] == Fe and lv[0]['eps'] == 0 and np.isclose(lv[0]['bp'], cfg['bp'])
        for l in lv:
            assert np.isclose(l['cp']**2 * cfg['Sp'], E.strength(line, gp, l['F']))
            assert np.isclose(l['cc']**2 * cfg['Sc'], E.strength(line, cfg['gc'], l['F']))
    p = E.params(E.paper_cases('D2')['A'], levels=E.excited_manifold('D2', 1, 1), sign=-1)
    assert p['levels'][0][1] == 1.0 and p['levels'][2][1] < 0       # sign acts on neighbours only


def test_multilevel_absorption_single_level_is_fano():
    p = E.params(E.paper_cases()['A'], s=1.0)
    d = np.linspace(-3, 3, 41)
    assert np.allclose(E.absorption(p, d), E.fano(p, d), rtol=1e-12)
    lv = [dict(F=2, eps=0.0, cp=1.0, cc=1.0, bp=0.5)]                # explicit one-level manifold
    p1 = E.params(E.paper_cases()['A'], s=1.0, levels=lv)
    p1['levels'].append((1e9, 0.0, 0.0, 0.5))                       # decoupled level: no effect
    assert np.allclose(E.absorption(p1, d), E.fano(p, d), rtol=1e-10)


def test_multilevel_me_far_neighbours_reproduce_three_level():
    fc = E.paper_cases()['B']
    three = E.solve_me(E.params(fc, s=1.0, N=12))
    lv = [dict(F=2, eps=0.0, cp=1.0, cc=1.0, bp=0.5), dict(F=1, eps=1e13, cp=0.5, cc=1.0, bp=0.2)]
    multi = E.solve_me(E.params(fc, s=1.0, N=12, levels=lv))
    assert np.isclose(multi['n'], three['n'], rtol=1e-5)
    assert np.isclose(multi['W'], three['W'], rtol=1e-5)


def test_lines_multilevel_study():
    d = load('lines_ml_B.dat')
    row = d[(d['line'] == 1) & (d['Fe'] == 2) & (d['gp'] == 1) & (d['sign'] == 1)][0]
    fc = E.paper_cases('D1')['B']
    p = E.params(fc, s=1.0, N=20, d2=row['d2opt'], levels=E.excited_manifold('D1', 2, 1))
    r = E.solve_me(p)
    assert np.isclose(r['W'] * fc['om'], row['Wopt'], rtol=1e-6)
    assert np.isclose(r['n'], row['nopt'], rtol=1e-6)
    dl = load('lines_delta_A.dat')
    assert set(np.unique(dl['DeltaG'])) == {1, 2, 3, 5, 8, 12}
    assert np.all(dl['okopt'] == (dl['nopt'] < 0.2))


# ---------------------------------------------------------------- eliminated / banded ME
def test_eliminated_me_matches_full():
    """Adiabatic elimination of |e> (Delta = 3 Gamma): < 1 % in n_ss and W."""
    for key in 'AB':
        p = E.params(E.paper_cases()[key], s=1.0, N=12)
        f, e = E.solve_me(p), E.solve_me(p, eliminate=True)
        assert abs(e['n'] / f['n'] - 1) < 1e-2 and abs(e['W'] / f['W'] - 1) < 1e-2


def test_banded_solver():
    p = E.params(E.paper_cases()['A'], s=1.0, N=16)
    e = E.solve_me(p, eliminate=True)
    full_band = E.solve_banded(p, B=16)                 # band holds everything: exact
    assert np.isclose(full_band['n'], e['n'], rtol=1e-10)
    assert np.isclose(full_band['W'], e['W'], rtol=1e-10)
    b6 = E.solve_banded(p, B=6)
    assert abs(b6['n'] / e['n'] - 1) < 1e-5 and abs(b6['W'] / e['W'] - 1) < 1e-5
    assert np.isclose(b6['P'].sum(), 1.0, atol=1e-10)


def test_convergence_results():
    d = np.genfromtxt(RESULTS / 'convergence' / 'nconv.dat', names=True)
    r = d[(d['case'] == 0) & (d['s'] == 1.0) & (d['method'] == 1)]
    assert np.isclose(r['n'][r['N'] == 160][0], 0.15447, atol=2e-5)
    assert np.isclose(r['W'][r['N'] == 160][0], 1614.1, rtol=1e-3)
    fc = E.paper_cases()['A']
    x = E.solve_banded(E.params(fc, s=1.0, N=40), B=6)
    row = r[r['N'] == 40][0]
    assert np.isclose(x['n'], row['n'], rtol=1e-8) and np.isclose(x['W'] * fc['om'], row['W'], rtol=1e-8)


# ---------------------------------------------------------------- numerics.toml
def test_numerics_config_rules():
    c = E.paper_cases()
    assert E.me_cutoffs(c['A']['eta_sp'], 1.0)['N'] == 120
    weak = E.me_cutoffs(c['A']['eta_sp'], 0.1)
    assert weak['band'] == 4 and weak['N'] == 400 and weak['status'] == 'approximate'
    assert E.me_cutoffs(c['A']['eta_sp'], 0.5)['N'] == 250
    assert E.me_cutoffs(c['B']['eta_sp'], 0.1)['N'] == 40
    with pytest.raises(ValueError):
        E.me_cutoffs(1.2, 1.0)                           # no rule covers it: must fail
    assert E.CSAT == E.CFG['grm']['csat'] and len(E.NODES) == E.CFG['emission']['nodes']


def test_steady_uses_config():
    fc = E.paper_cases()['B']
    r = E.steady(fc, 1.0)
    p = E.params(fc, s=1.0, N=r['N'])
    ref = E.solve_banded(p, B=r['band'])
    assert r['N'] == 40 and r['band'] == 6
    assert np.isclose(r['n'], ref['n'], rtol=1e-12) and np.isclose(r['W'], ref['W'], rtol=1e-12)

if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)


def test_banded_evolution_matches_full():
    """Second order in h and agreement with the full BDF2 evolution (N = 10)."""
    N = 10
    p = E.params(E.paper_cases()['B'], s=0.2, N=N)
    th = E.thermal(1.0, N)
    L, d, b = E.liouvillian(p)
    rho0 = np.zeros((d, d), complex)
    rho0[np.arange(N), np.arange(N)] = th
    t, full = E.evolve(L, rho0, 5.0, E.number_op(b), h=0.05, n_out=5)
    _, n1 = E.banded_evolve(p, th, 5.0, B=N, h=0.0625, n_out=5)
    _, n2 = E.banded_evolve(p, th, 5.0, B=N, h=0.03125, n_out=5)
    _, n3 = E.banded_evolve(p, th, 5.0, B=N, h=0.015625, n_out=5)
    assert np.max(abs(n1 / full - 1)) < 1e-3
    assert abs(n1[-1] - n2[-1]) > 3 * abs(n2[-1] - n3[-1])     # second order: ratio -> 4


def test_cutoff_rules_depend_on_branching():
    eta = 0.71
    assert E.me_cutoffs(eta, 1.0)['N'] == 120
    assert E.me_cutoffs(eta, 0.1, 5 / 6)['N'] == 120
    assert E.me_cutoffs(eta, 1.0, 1 / 6)['status'] == 'unconverged'
    assert E.me_cutoffs(0.34, 1.0, 1 / 6)['status'] == 'validated'
    assert E.me_cutoffs(eta, 1.0, 0.5, 8.0)['N'] == 400
    assert E.me_cutoffs(eta, 1.0, 5 / 6, 8.0)['N'] == 120


def test_rates_match_direct_time_integration():
    """W and W0 of solve_banded against integrals of the exact dynamics (small N)."""
    from scipy.linalg import expm
    N = 8
    p = E.params(E.paper_cases()['B'], s=0.3, N=N)
    r = E.solve_banded(p, B=N)
    M, ri, ci, n_of = E.banded_generator(p, B=N)
    diag = np.nonzero(ri == ci)[0]
    nv = np.zeros(ri.size)
    nv[diag] = n_of[ri[diag]]
    pv = np.zeros(ri.size)
    pv[diag[n_of[ri[diag]] == 0]] = 1.0
    y = np.zeros(ri.size, complex)
    y[int(np.nonzero((ri == 1) & (ci == 1))[0][0])] = 1.0
    dt = 0.5
    U = expm(M * dt)
    tn, tp = 0.0, 0.0
    for k in range(40000):
        a = float(np.real(nv @ y)) - r['n']
        b = r['P0'] - float(np.real(pv @ y))
        tn += a * dt * (0.5 if k == 0 else 1.0)
        tp += b * dt * (0.5 if k == 0 else 1.0)
        y = U @ y
    assert np.isclose(1 / (tn / (1 - r['n'])), r['W'], rtol=1e-3)
    assert np.isclose(1 / (tp / r['P0']), r['W0'], rtol=1e-3)

