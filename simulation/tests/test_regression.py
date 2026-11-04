"""Fast regression tests (~10 s): the package must reproduce stored results.

Run from the repository root:  python -m pytest simulation/tests -q
"""
import sys
from pathlib import Path

import numpy as np

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
    d = load('sat_A.dat')
    fc = E.paper_cases()['A']
    i = int(np.argmin(np.abs(d['s'] - 0.1)))
    r = E.solve_me(E.params(fc, s=d['s'][i], N=20))
    assert np.isclose(r['n'], d['nME'][i], rtol=1e-6)
    assert np.isclose(r['W'] * fc['om'], d['WME'][i], rtol=1e-6)


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

if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
