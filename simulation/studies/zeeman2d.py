"""Both radial axes with the Zeeman structure of 87Rb (eitsim.zeeman2d): LP11 switching against TE01.
Meant for the cluster (simulation/cluster/README.md); only `check`, `smoke` and `cost` are cheap.

  check, check_field  reduction tests (zero field; Zeeman Lambda at mu_B B/hbar = 40 nu_r): with N_y = 1 and a one-axis probe the 2D code must reproduce the validated
           one-axis Zeeman solver (zeeman.solve) for the best hyperfine scheme (D1 F'=1, lin/pi,
           repumper); direct dense solver, N_x = 6, seconds
  smoke    every scheme at N = 2 per axis: builds, a few RK4 steps, trace conservation
  cost     dimension, number of jump operators and flops per step for each scheme at N = 4..10
  run      one simulation (one array task): --scheme --case --s --d1 --d2 --gamma --zeeman --N --tau
           --tmax ... -> results/zeeman2d/runs/<tag>.json (skipped if it exists)
  list     write a scan list (one `run` per line) for scripts in simulation/cluster/
  collect  results/zeeman2d/runs/*.json (or --out DIR) -> results/zeeman2d.dat

Protocol of `run`: RK4 from |1,1> with the internal state uniform over F=2 (or the target sublevel,
--init target), until tmax; P00ss is the mean over the last window (the last full switching period, or
50/nu_r), W00 = P00ss / int (P00ss - P00) dt, as in studies/two_axes.py; drift = change of P00 over
the last window (should be << 1e-3).

Usage: python zeeman2d.py check | smoke | cost | run ... | list [--grid main|test] | collect
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

from _common import E, RESULTS, save, log
from eitsim import zeeman as Z
from eitsim import zeeman2d as Z2

RUNS = RESULTS / 'zeeman2d' / 'runs'
CLUSTER = Path(__file__).resolve().parents[1] / 'cluster'
SCHEMES = ('lp11_hf', 'lp11_zl', 'te01_bz', 'te01_bx')


def setup(case, s, DG=3.0, N=6):
    fc = E.paper_cases()[case]
    return Z.setup(fc, s, DG, 0.0, N)


# ------------------------------------------------------------------ check against the 1D solver
def check():
    from zeeman_scan import POL, repumpers
    fc = E.paper_cases()['A']
    N, s, DG, d2, gamma = 6, 1.0, 3.0, 0.45, 0.1
    p = Z.setup(fc, s, DG, d2, N)
    st = Z.structure_lambda('D1', 1, 2, POL['lin'], POL['pi'])
    r1 = Z.solve(st, p, B=N, repump=repumpers(st, 'lin', 'pi', p, gamma))
    lin, pi = (0.0, -1j, 0.0), (0.0, 0.0, 1.0)
    fields = [Z2.field(2, [(lin, 'x', 1.0)], p['Wp'], d2), Z2.field(1, [(pi, 'lp01', 1.0)], p['Wc'], 0.0),
              Z2.field(2, [(pi, 'lp01', 1.0)], Z2.repump_amp(p, gamma), 0.0, incoherent=True),
              Z2.field(1, [(lin, 'lp01', 1.0)], Z2.repump_amp(p, gamma), 0.0, incoherent=True)]
    m = Z2.Model(fields, p, N, 1, bdir='z', zeeman=0.0, fes=(1,), nodes=(8, 16))
    r2 = m.steady_direct()
    log('check 1D solver: P0 %.5f n %.5f W0 %.5f' % (r1['P0'], r1['n'], r1['W0']))
    log('check 2D code, N_y = 1: P0 %.5f n %.5f W0 %.5f (frame mismatch %.1e)' % (r2['P00'], r2['nx'], r2['W00'], m.frame_mismatch))
    ok = abs(r2['P00'] - r1['P0']) < 3e-3 and abs(r2['W00'] / r1['W0'] - 1) < 0.03
    log('check', 'PASS' if ok else 'FAIL')
    return ok


def check_field():
    """Finite field: Zeeman Lambda inside F = 2 (sigma- probe, sigma+ coupling, F = 1 repumper sigma+),
    one axis, against zeeman.structure_zeeman_lambda (which writes the linear Zeeman shift into the
    frame of the two beams)."""
    fc = E.paper_cases()['A']
    N, s, DG, d2, gamma, Zf = 6, 1.0, 3.0, 0.2, 0.3, 40.0
    p = Z.setup(fc, s, DG, d2, N)
    st = Z.structure_zeeman_lambda()
    na, nE = len(st['ground']), len(st['excited'])
    g = Z.field_matrix('D1', 1, st['excited'], {1: 1.0}) / np.sqrt(Z._norm('D1', 1, 2))
    gm = np.zeros((na, nE), complex)
    gm[5:] = g
    r1 = Z.solve(st, p, B=N, repump=[dict(gm=gm, W=Z2.repump_amp(p, gamma), detuning=p['D'])])
    fields = [Z2.field(2, [(Z2.SIGM_Z, 'x', 1.0)], p['Wp'], 0.0), Z2.field(2, [(Z2.SIGP_Z, 'lp01', 1.0)], p['Wc'], Zf - d2),
              Z2.field(1, [(Z2.SIGP_Z, 'lp01', 1.0)], Z2.repump_amp(p, gamma), 0.0, incoherent=True)]
    m = Z2.Model(fields, p, N, 1, bdir='z', zeeman=Zf, fes=(2,), nodes=(8, 16))
    r2 = m.steady_direct()
    log('check_field 1D solver: P0 %.5f n %.5f W0 %.5f' % (r1['P0'], r1['n'], r1['W0']))
    log('check_field 2D code, N_y = 1, Z = %g: P0 %.5f n %.5f W0 %.5f (frame mismatch %.1e, %d classes)' % (
        Zf, r2['P00'], r2['nx'], r2['W00'], m.frame_mismatch, len(m.classes)))
    ok = abs(r2['P00'] - r1['P0']) < 5e-3 and abs(r2['W00'] / r1['W0'] - 1) < 0.05
    log('check_field', 'PASS' if ok else 'FAIL')
    return ok


def smoke():
    for name in SCHEMES:
        p = setup('A', 1.0, N=2)
        m = Z2.build(name, p, 2, 2, zeeman=40.0, d1=0.3, d2=0.3, gamma=0.1)
        r0 = m.initial(1, 1)
        t, P, r = m.evolve(r0, 2.0, dt=0.1, tau=1.0 if m.stages > 1 else None)
        o = m.observables(r)
        log('smoke %-8s d %d classes %d jumps %d mismatch %.1e: trace %.8f P00 %.4f' % (
            name, m.d, len(m.classes), m.cost()['n_jumps'], m.frame_mismatch, o['trace'], o['P00']))


def cost():
    for name in SCHEMES:
        for N in (4, 6, 8, 10):
            p = setup('A', 1.0, N=N)
            m = Z2.build(name, p, N, N, zeeman=40.0, d1=0.3, d2=0.3, gamma=0.1)
            c = m.cost()
            dt = min(0.1, 1.2 / (2 * N))
            steps = 1500 / dt
            log('cost %-8s N %2d: d %4d, jumps %3d, %.1e flop/step, %.0f MB operators; 1500/nu_r = %.0f steps = %.1e flop' % (
                name, N, c['d'], c['n_jumps'], c['flops_per_step'], c['op_memory_MB'], steps, steps * c['flops_per_step']))


# ------------------------------------------------------------------ one run
def tag_of(a):
    return '%s_%s_s%g_d%g_%g_g%g_Z%g_N%d_tau%g_wc%g_DG%g' % (a.scheme, a.case, a.s, a.d1, a.d2, a.gamma, a.zeeman, a.N, a.tau, a.wc, a.DG)


def run(argv):
    ap = argparse.ArgumentParser(prog='zeeman2d.py run')
    ap.add_argument('--scheme', required=True, choices=SCHEMES)
    ap.add_argument('--case', default='A', choices=('A', 'B'))
    ap.add_argument('--s', type=float, default=1.0)
    ap.add_argument('--DG', type=float, default=3.0, help='Delta / Gamma')
    ap.add_argument('--d1', type=float, default=0.0, help='two-photon detuning, channel 1 (nu_r)')
    ap.add_argument('--d2', type=float, default=0.0, help='two-photon detuning, channel 2 (nu_r)')
    ap.add_argument('--gamma', type=float, default=0.0, help='repumper rate (nu_r)')
    ap.add_argument('--zeeman', type=float, default=40.0, help='mu_B B / hbar in units of nu_r')
    ap.add_argument('--wc', type=float, default=1.0, help='coupling Rabi scale')
    ap.add_argument('--N', type=int, default=6, help='Fock states per axis')
    ap.add_argument('--tau', type=float, default=32.0, help='stage length for LP11 switching (1/nu_r)')
    ap.add_argument('--tmax', type=float, default=1500.0)
    ap.add_argument('--dt', type=float, default=0.0)
    ap.add_argument('--init', default='F2', choices=('F2', 'target'))
    ap.add_argument('--tol', type=float, default=2.0)
    ap.add_argument('--out', default=str(RUNS))
    a = ap.parse_args(argv)
    out = Path(a.out) / (tag_of(a) + '.json')
    if out.exists():
        log('exists, skipped:', out)
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    p = setup(a.case, a.s, a.DG, a.N)
    m = Z2.build(a.scheme, p, a.N, a.N, tol=a.tol, zeeman=a.zeeman, d1=a.d1, d2=a.d2, gamma=a.gamma, wc=a.wc)
    log('run', tag_of(a), 'd', m.d, 'classes', len(m.classes), 'frame mismatch %.1e' % m.frame_mismatch,
        'class spread %.2f' % m.class_spread)
    internal = None
    if a.init == 'target':
        internal = np.zeros(m.na)
        internal[m.target] = 1.0
    r0 = m.initial(1, 1, internal)
    tau = a.tau if m.stages > 1 else None
    tmax = a.tmax if tau is None else 2 * tau * np.ceil(a.tmax / (2 * tau))
    dt = a.dt or min(0.1, 1.2 / (2 * a.N))
    t, P, r = m.evolve(r0, tmax, dt=dt, tau=tau, log=log)
    win = int(round((2 * tau if tau else 50.0) / dt))
    pss = float(P[-win:].mean())
    W00 = float(pss / (np.sum(pss - P) * dt))
    drift = float(P[-1] - P[-win - 1])
    o = m.observables(r)
    res = dict(vars(a), tag=tag_of(a), P00=pss, W00=W00, drift=drift, nx=o['nx'], ny=o['ny'], P0x=o['P0x'], P0y=o['P0y'],
               trace=o['trace'], pint=list(o['pint']), frame_mismatch=m.frame_mismatch, class_spread=m.class_spread,
               om=float(p['om']), W00_per_s=W00 * float(p['om']), d=m.d, runtime_s=time.time() - t0,
               P00_t=list(P[::max(1, len(P) // 400)]), t=list(t[::max(1, len(t) // 400)]))
    out.write_text(json.dumps(res, indent=1))
    log('done: P00 %.4f W00 %.4f nu_r (%.0f s^-1), drift %.1e, nx %.3f ny %.3f, %.0f s' % (
        pss, W00, W00 * p['om'], drift, o['nx'], o['ny'], time.time() - t0))


# ------------------------------------------------------------------ scan lists
GRIDS = {
    # quick test of the cluster set-up: two short runs
    'test': [dict(scheme='te01_bz', case='A', s=1.0, d1=0.3, d2=0.3, gamma=0.0, N=4, tmax=300.0),
             dict(scheme='lp11_hf', case='A', s=1.0, d1=0.45, d2=0.0, gamma=0.1, N=4, tmax=300.0)],
}


def main_grid():
    """Case A only (the fibre of interest)."""
    rows = []
    for case in ('A',):
        for s in (0.5, 1.0):
            # TE01, both field directions: two-photon detunings per channel, optional repumper
            for sch in ('te01_bz', 'te01_bx'):
                for d1, d2 in itertools.product((0.0, 0.4, 0.8), repeat=2):
                    for g in (0.0, 0.1):
                        rows.append(dict(scheme=sch, case=case, s=s, d1=d1, d2=d2, gamma=g))
            # LP11 switched: hyperfine Lambda (best 1D: repumper 0.03-0.1, d 0.45-0.75) and Zeeman Lambda
            for d1 in (0.3, 0.6, 0.9):
                for g in (0.03, 0.1):
                    for tau in (16.0, 32.0):
                        rows.append(dict(scheme='lp11_hf', case=case, s=s, d1=d1, gamma=g, tau=tau))
            for d1 in (-0.4, 0.0, 0.4):
                for g in (0.3, 1.0):
                    rows.append(dict(scheme='lp11_zl', case=case, s=s, d1=d1, gamma=g, tau=32.0))
    return rows


def write_list(argv):
    ap = argparse.ArgumentParser(prog='zeeman2d.py list')
    ap.add_argument('--grid', default='main', choices=('main', 'test', 'convergence'))
    ap.add_argument('--N', type=int, default=6)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    if a.grid == 'main':
        rows = main_grid()
        for r in rows:
            r.setdefault('N', a.N)
    elif a.grid == 'convergence':
        rows = [dict(scheme=s, case=c, s=1.0, d1=0.4, d2=0.4, gamma=0.1 if s.startswith('lp11') else 0.0, N=N)
                for s in ('te01_bz', 'lp11_hf') for c in ('A',) for N in (6, 8, 10)]
    else:
        rows = GRIDS['test']
    out = Path(a.out) if a.out else CLUSTER / 'scans' / ('zeeman2d_%s.txt' % a.grid)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as f:
        f.write('# one line per run: arguments of "zeeman2d.py run" (written by zeeman2d.py list --grid %s)\n' % a.grid)
        for r in rows:
            f.write('zeeman2d.py run ' + ' '.join('--%s %s' % (k, v) for k, v in r.items()) + '\n')
    log('wrote %d runs to %s' % (len(rows), out))


def collect(argv=()):
    runs = Path(argv[argv.index('--out') + 1]) if '--out' in argv else RUNS
    keys = ['scheme', 'case', 's', 'd1', 'd2', 'gamma', 'zeeman', 'N', 'tau', 'wc', 'DG', 'P00', 'W00', 'W00_per_s', 'drift',
            'nx', 'ny', 'frame_mismatch', 'runtime_s']
    code = {s: i for i, s in enumerate(SCHEMES)}
    rows = []
    for fjs in sorted(runs.glob('*.json')):
        r = json.loads(fjs.read_text())
        rows.append([code[r['scheme']], 0 if r['case'] == 'A' else 1] + [r[k] for k in keys[2:]])
    keys = ['scheme', 'case'] + keys[2:]
    save('zeeman2d.dat', keys, rows)
    log('collected %d runs (scheme codes: %s)' % (len(rows), ', '.join('%d=%s' % (v, k) for k, v in code.items())))


if __name__ == '__main__':
    cmd, rest = sys.argv[1], sys.argv[2:]
    if cmd == 'run':
        run(rest)
    elif cmd == 'list':
        write_list(rest)
    elif cmd == 'collect':
        collect(rest)
    else:
        {'check': check, 'check_field': check_field, 'smoke': smoke, 'cost': cost}[cmd]()
