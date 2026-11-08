"""Numerical settings from simulation/numerics.toml (truncations, cutoffs, quadrature).

CFG is the parsed file; me_cutoffs selects the Fock truncation and the band of
the steady-state solver for a given recoil Lamb-Dicke parameter and saturation.
"""
import tomllib
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / 'numerics.toml'

with open(PATH, 'rb') as f:
    CFG = tomllib.load(f)


def me_cutoffs(eta_sp, s, bp=0.5, delta_g=3.0):
    """First steady_state.rule with eta_sp <= eta_sp_max, s_min <= s < s_max and,
    for the optional keys (default: any value), bp_min <= bp <= bp_max and
    delta_min <= delta_g (Delta / Gamma).

    Returns dict(N, band, status, name). Raises ValueError if no rule applies:
    such parameters need a convergence check and a new rule in numerics.toml.
    """
    for r in CFG['steady_state']['rule']:
        if (eta_sp <= r['eta_sp_max'] and r['s_min'] <= s < r['s_max']
                and r.get('bp_min', 0.0) <= bp <= r.get('bp_max', 1.0)
                and r.get('delta_min', 0.0) <= delta_g):
            return dict(N=r['N'], band=r['band'], status=r['status'], name=r['name'])
    raise ValueError('no steady_state rule in %s for eta_sp = %.3f, s = %g, bp = %.3f, Delta = %g Gamma'
                     % (PATH.name, eta_sp, s, bp, delta_g))
