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


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
