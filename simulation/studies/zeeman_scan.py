"""Zeeman structure of the D1 Lambda configurations, one radial axis, case A trap.

Dark states: a fixed polarisation on F -> F' leaves dark combinations if 2F + 1 > rank of the
coupling, whatever the motion, so atoms pumped into them stop cooling. They must be remixed.

  dark  count the dark combinations of the probe and coupling legs     -> printed
  scan  configurations x polarisation pairs x remixing x delta_e, s = 1, 3 Gamma, N = 12, B = 2;
        rows with NaN have an unmixed dark state and no unique steady state
                                             -> results/zeeman_scan.dat (~1 h, 4 workers)

  refine  the best hyperfine Lambda (D1 F'=1, lin perp lin) and the Zeeman Lambda inside F = 2
          (Qiu et al. 2026; sigma- probe, sigma+ coupling, F = 1 repumper) in delta_e, s, remixing,
          coupling strength, trap                 -> results/zeeman_refine.dat (~30 min, 4 workers)

Usage: python zeeman_scan.py [dark] [scan] [refine]
"""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import E, save, log
from eitsim import zeeman as Z

S2 = 2**-0.5
POL = {'pi': {0: 1.0}, 'lin': {1: S2, -1: S2}, 'sp': {1: 1.0}, 'sm': {-1: 1.0}}
COMP = {'pi': 'lin', 'lin': 'pi', 'sp': 'sm', 'sm': 'sp'}
PAIRS = (('pi', 'pi'), ('pi', 'lin'), ('lin', 'pi'), ('lin', 'lin'), ('sp', 'sp'), ('sp', 'sm'))
CONFIGS = ((1, 2), (2, 1), (2, 2))
REMIX = (('Lx', 0.03), ('Lx', 0.1), ('Lx', 0.3), ('Lz', 0.1), ('Lz', 0.3), ('rep', 0.03), ('rep', 0.1), ('rep', 0.3))
D2S = np.linspace(-0.45, 0.45, 7)
NAME = {'pi': 0, 'lin': 1, 'sp': 2, 'sm': 3}
KIND = {'Lx': 0, 'Lz': 1, 'rep': 2}


def dark():
    for Fe, gp in CONFIGS:
        gc = 3 - gp
        ex = [(Fe, m) for m in range(-Fe, Fe + 1)]
        for name, pol in POL.items():
            for leg, F in (('probe', gp), ('coupling', gc)):
                g = Z.field_matrix('D1', F, ex, pol)
                rk = np.linalg.matrix_rank(g, tol=1e-8)
                log('dark', "D1 F'=%d, probe on F=%d, %s on F=%d, %s:" % (Fe, gp, leg, F, name), '%d dark of %d' % (2 * F + 1 - rk, 2 * F + 1))


def repumpers(st, pp, pc, p, gamma):
    """Incoherent repumpers on both legs with the complementary polarisation, at the same detuning,
    each giving a scattering rate gamma per sublevel on average."""
    out = []
    Wr = np.sqrt(4 * gamma * (p['D']**2 + p['G']**2 / 4) / p['G'])
    na, nE = len(st['ground']), len(st['excited'])
    for F, comp, i0 in ((st['gp'], COMP[pp], 0), (3 - st['gp'], COMP[pc], st['nP'])):
        g = Z.field_matrix('D1', F, st['excited'], POL[comp]) / np.sqrt(Z._norm('D1', F, st['Fe']))
        gm = np.zeros((na, nE), complex)
        gm[i0:i0 + g.shape[0]] = g
        out.append(dict(gm=gm, W=Wr, detuning=p['D']))
    return out


def task(a):
    Fe, gp, pp, pc, kind, val, d2, N, B, s, DG = a
    fc = E.paper_cases()['A']
    st = Z.structure_lambda('D1', Fe, gp, POL[pp], POL[pc])
    p = Z.setup(fc, s, DG, d2, N)
    bv, rp = (0.0, 0.0), ()
    if kind == 'Lx':
        bv = (0.0, val)
    elif kind == 'Lz':
        bv = (val, 0.0)
    else:
        rp = repumpers(st, pp, pc, p, val)
    r = Z.solve(st, p, B=B, b_vec=bv, repump=rp)
    return [Fe, gp, NAME[pp], NAME[pc], KIND[kind], val, d2, N, B, s, DG, r['n'], r['P0'], r['W'], r['W0'], r['pint'][:st['nP']].sum()]


def scan():
    tasks = [(Fe, gp, pp, pc, k, v, d2, 12, 2, 1.0, 3.0) for Fe, gp in CONFIGS for pp, pc in PAIRS for k, v in REMIX for d2 in D2S]
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(task, tasks, chunksize=4))
    save('zeeman_scan.dat', ['Fe', 'gp', 'pp', 'pc', 'kind', 'val', 'd2', 'N', 'B', 's', 'DeltaG', 'n', 'P0', 'W', 'W0', 'pop_probe'], rows)
    log('scan done', len(rows))


def _ref_task(a):
    scheme, nw, pp, pc, kind, val, wc, d2, s = a
    fc = E.reduced_case(nw)
    N, B = 12, 2
    p = Z.setup(fc, s, 3.0, d2, N, wc_scale=wc)
    bv, rp = (0.0, 0.0), ()
    if scheme == 0:
        st = Z.structure_lambda('D1', 1, 2, POL[pp], POL[pc])
        if kind == 'Lx':
            bv = (0.0, val)
        elif kind == 'Lz':
            bv = (val, 0.0)
        else:
            rp = repumpers(st, pp, pc, p, val)
    else:
        st = Z.structure_zeeman_lambda()
        na, nE = len(st['ground']), len(st['excited'])
        g = Z.field_matrix('D1', 1, st['excited'], {1: 1.0}) / np.sqrt(Z._norm('D1', 1, 2))
        gm = np.zeros((na, nE), complex)
        gm[5:] = g
        Wr = np.sqrt(4 * val * (p['D']**2 + p['G']**2 / 4) / p['G'])
        rp = [dict(gm=gm, W=Wr, detuning=p['D'])]
    r = Z.solve(st, p, B=B, b_vec=bv, repump=rp)
    return [scheme, nw, NAME[pp], NAME[pc], KIND[kind], val, wc, d2, s, N, B, r['n'], r['P0'], r['W'], r['W0']]


def refine():
    A = E.paper_cases()['A']['eta_sp']**-2
    tasks = []
    for pp, pc in (('lin', 'pi'), ('pi', 'lin')):
        for kind, val in (('rep', 0.01), ('rep', 0.03), ('rep', 0.1), ('Lx', 0.1), ('Lx', 0.3)):
            for s in (0.3, 1.0, 3.0):
                for d2 in np.arange(0.0, 1.21, 0.15):
                    tasks.append((0, A, pp, pc, kind, val, 1.0, d2, s))
    for nw in (A, 5.0):
        for val in (0.1, 0.3, 1.0):
            for wc in (0.7, 1.0, 1.4, 2.0):
                for s in (0.3, 1.0, 3.0):
                    for d2 in np.arange(-1.0, 1.01, 0.2):
                        tasks.append((1, nw, 'sm', 'sp', 'rep', val, wc, d2, s))
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(_ref_task, tasks, chunksize=4))
    save('zeeman_refine.dat', ['scheme', 'nu_wrec', 'pp', 'pc', 'kind', 'val', 'wc', 'd2', 's', 'N', 'B', 'n', 'P0', 'W', 'W0'], rows)
    log('refine done', len(rows))


if __name__ == '__main__':
    for name in sys.argv[1:]:
        {'dark': dark, 'scan': scan, 'refine': refine}[name]()
