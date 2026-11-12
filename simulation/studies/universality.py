"""Strengthening the power-law universality P(n) ~ n^(-1/h) (cluster; see simulation/cluster/README.md).

  fixedh  steady state at FIXED h with Delta and s varied: the exponent should stay at ~1/h, while the
          window moves, n* ~ 1/(c_sat s) (start of the saturation plateau) and n_c = 2 (nu_r/Gamma_n)^2 / s
          (onset of the fast channel). Reduced case, b_p = 1/2, nu_r/omega_rec = alpha/(b_p h).
          --h --DG --s [--N 800 --B 4] -> results/universality/fixedh/<tag>.json (P(n), fit, window)
  drift   local exponent over 2-3 decades without large-N steady states: drift v(n) and energy diffusion
          D2(n) = dVar(n)/dt from short evolutions of Fock states n0 (eitsim banded_fock_drift); the
          zero-flux Fokker-Planck solution gives the local exponent
              gamma(n) = d ln P / d ln n = 2 n v / D2 - d ln D2 / d ln n,
          to be compared with -1/h on the plateau. --h --DG --s [--N 1000 --B 4 --n0max 640]
          -> results/universality/drift/<tag>.json
  list    scan lists for the cluster -> simulation/cluster/scans/universality_{fixedh,drift}.txt
  collect -> results/universality_fixedh.dat, results/universality_drift.dat

Memory (dense banded LU, two ground levels): S ~ 36 N for B = 4, ~ 68 N for B = 8; S^2 x 16 B, twice for drift.
N = 800, B = 4: ~13 GB; N = 1000, B = 4: ~42 GB peak; N = 600, B = 8: ~53 GB peak.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from _common import E, RESULTS, save, log

OUT = RESULTS / 'universality'
CLUSTER = Path(__file__).resolve().parents[1] / 'cluster'
BP = 0.5
CSAT = 2.25
N0_GRID = (5, 7, 10, 14, 20, 28, 40, 56, 80, 113, 160, 226, 320, 452, 640)


def nw_of(h, bp=BP):
    return E.ALPHA / (bp * h)


def tag(kind, a):
    return '%s_h%g_DG%g_s%g_N%d_B%d' % (kind, a.h, a.DG, a.s, a.N, a.B)


def window(p, s, N):
    nstar = 1.0 / (CSAT * s)
    nc = 2.0 / (p['Gn']**2 * s)
    return nstar, nc, max(5.0, 3 * nstar), min(N / 4.0, nc / 3.0)


def fixedh(argv):
    ap = argparse.ArgumentParser(prog='universality.py fixedh')
    ap.add_argument('--h', type=float, required=True)
    ap.add_argument('--DG', type=float, required=True)
    ap.add_argument('--s', type=float, required=True)
    ap.add_argument('--N', type=int, default=800)
    ap.add_argument('--B', type=int, default=4)
    a = ap.parse_args(argv)
    out = OUT / 'fixedh' / (tag('fixedh', a) + '.json')
    if out.exists():
        log('exists, skipped:', out)
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    nw = nw_of(a.h)
    fc = E.reduced_case(nw)
    p = E.params(fc, s=a.s, Delta_G=a.DG, bp=BP)
    r = E.steady(fc, a.s, cut=(a.N, a.B), Delta_G=a.DG, bp=BP)
    P = np.asarray(r['P'])
    k = np.arange(a.N)
    nstar, nc, lo, hi = window(p, a.s, a.N)
    sel = (k >= lo) & (k <= hi) & (P > 0)
    slope = np.polyfit(np.log(k[sel]), np.log(P[sel]), 1)[0] if sel.sum() >= 5 and hi / lo >= 2 else np.nan
    res = dict(h=a.h, DG=a.DG, s=a.s, N=a.N, B=a.B, nu_wrec=nw, bp=BP, P0=r['P0'], n=r['n'], W0=r['W0'],
               nstar=nstar, nc=nc, fit_lo=lo, fit_hi=hi, slope=float(slope), decades=float(np.log10(hi / lo)) if hi > lo else 0.0,
               Gn=float(p['Gn']), P=P.tolist(), runtime_s=time.time() - t0)
    out.write_text(json.dumps(res))
    log('fixedh h %g DG %g s %g: P0 %.4f n %.3f slope %.3f (-1/h %.3f) over n = %.1f-%.1f (%.2f decades), n* %.2f n_c %.0f, %.0f s' % (
        a.h, a.DG, a.s, r['P0'], r['n'], slope, -1 / a.h, lo, hi, res['decades'], nstar, nc, time.time() - t0))


def drift(argv):
    ap = argparse.ArgumentParser(prog='universality.py drift')
    ap.add_argument('--h', type=float, required=True)
    ap.add_argument('--DG', type=float, default=8.0)
    ap.add_argument('--s', type=float, default=1.0)
    ap.add_argument('--N', type=int, default=1000)
    ap.add_argument('--B', type=int, default=4)
    ap.add_argument('--n0max', type=int, default=640)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--t1', type=int, default=100)
    ap.add_argument('--t2', type=int, default=400)
    a = ap.parse_args(argv)
    out = OUT / 'drift' / (tag('drift', a) + '.json')
    if out.exists():
        log('exists, skipped:', out)
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    nw = nw_of(a.h)
    fc = E.reduced_case(nw)
    p = E.params(fc, s=a.s, Delta_G=a.DG, bp=BP, N=a.N)
    n0s = [n for n in N0_GRID if n <= a.n0max]
    vd = E.banded_fock_drift(p, n0s, B=a.B, h=a.dt, t1=a.t1, t2=a.t2)
    v, D2 = vd[:, 0], vd[:, 1]
    ln = np.log(np.asarray(n0s, float))
    dlnD = np.gradient(np.log(np.abs(D2)), ln)
    gamma = 2 * np.asarray(n0s) * v / D2 - dlnD
    nstar, nc, _, _ = window(p, a.s, a.N)
    res = dict(h=a.h, DG=a.DG, s=a.s, N=a.N, B=a.B, nu_wrec=nw, bp=BP, n0=n0s, v=v.tolist(), D2=D2.tolist(),
               gamma=gamma.tolist(), nstar=nstar, nc=nc, Gn=float(p['Gn']), dt=a.dt, t1=a.t1, t2=a.t2, runtime_s=time.time() - t0)
    out.write_text(json.dumps(res))
    for n0, vv, dd, g in zip(n0s, v, D2, gamma):
        log('drift h %g DG %g n0 %4d: v %+.4e D2 %.4e local exponent %.3f (-1/h %.3f)' % (a.h, a.DG, n0, vv, dd, g, -1 / a.h))
    log('done in %.0f s' % (time.time() - t0))


def write_lists(argv):
    (CLUSTER / 'scans').mkdir(parents=True, exist_ok=True)
    f1 = CLUSTER / 'scans' / 'universality_fixedh.txt'
    with open(f1, 'w') as f:
        f.write('# fixed h, Delta and s varied (universality.py list); N = 800, B = 4\n')
        for h in (0.4, 0.7):
            for DG in (3.0, 5.0, 8.0, 12.0):
                for s in (0.1, 0.3, 1.0, 3.0):
                    f.write('universality.py fixedh --h %g --DG %g --s %g --N 800 --B 4\n' % (h, DG, s))
    f2 = CLUSTER / 'scans' / 'universality_drift.txt'
    with open(f2, 'w') as f:
        f.write('# local exponent from drift and diffusion (universality.py list)\n')
        for h in (0.3, 0.5, 0.7, 1.0):
            f.write('universality.py drift --h %g --DG 8 --s 1 --N 1000 --B 4 --n0max 640\n' % h)
        f.write('universality.py drift --h 0.4 --DG 3 --s 1 --N 1000 --B 4 --n0max 640\n')
        f.write('universality.py drift --h 0.5 --DG 8 --s 1 --N 600 --B 4 --n0max 320\n')     # band check: same, B = 4
        f.write('universality.py drift --h 0.5 --DG 8 --s 1 --N 600 --B 8 --n0max 320\n')     # band check: B = 8
    log('wrote', f1, 'and', f2)


def collect(argv):
    rows = []
    for fj in sorted((OUT / 'fixedh').glob('*.json')):
        r = json.loads(fj.read_text())
        rows.append([r[k] for k in ('h', 'DG', 's', 'N', 'B', 'P0', 'n', 'W0', 'nstar', 'nc', 'fit_lo', 'fit_hi', 'slope', 'decades')])
    if rows:
        save('universality_fixedh.dat', ['h', 'DG', 's', 'N', 'B', 'P0', 'n', 'W0', 'nstar', 'nc', 'fit_lo', 'fit_hi', 'slope', 'decades'], rows)
    rows = []
    for fj in sorted((OUT / 'drift').glob('*.json')):
        r = json.loads(fj.read_text())
        for n0, v, d, g in zip(r['n0'], r['v'], r['D2'], r['gamma']):
            rows.append([r['h'], r['DG'], r['s'], r['N'], r['B'], n0, v, d, g])
    if rows:
        save('universality_drift.dat', ['h', 'DG', 's', 'N', 'B', 'n0', 'v', 'D2', 'gamma'], rows)
    log('collected')


if __name__ == '__main__':
    {'fixedh': fixedh, 'drift': drift, 'list': write_lists, 'collect': collect}[sys.argv[1]](sys.argv[2:])
