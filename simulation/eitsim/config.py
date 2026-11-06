"""Numerical settings from simulation/numerics.toml (truncations, cutoffs, quadrature).

CFG is the parsed file; me_cutoffs selects the Fock truncation and the band of
the steady-state solver for a given recoil Lamb-Dicke parameter and saturation.
"""
import tomllib
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / 'numerics.toml'

with open(PATH, 'rb') as f:
    CFG = tomllib.load(f)


def me_cutoffs(eta_sp, s):
    """First steady_state.rule with eta_sp <= eta_sp_max and s_min <= s < s_max.

    Returns dict(N, band, status, name). Raises ValueError if no rule applies:
    such parameters need a convergence check and a new rule in numerics.toml.
    """
    for r in CFG['steady_state']['rule']:
        if eta_sp <= r['eta_sp_max'] and r['s_min'] <= s < r['s_max']:
            return dict(N=r['N'], band=r['band'], status=r['status'], name=r['name'])
    raise ValueError('no steady_state rule in %s for eta_sp = %.3f, s = %g'
                     % (PATH.name, eta_sp, s))
