"""Cooling dynamics (paper Fig. 5).

  relaxation  case A, s = 0.2, from <n> = 2: master equation (BDF2, ~6 min),
              GRM (calibrated and c_sat = 1) and simple model -> results/dyn_A.dat
  capture     cooling from 10 uK with the calibrated GRM, four protocols
              (~2 min)                                       -> results/capture_{A,B}.dat

Usage: python dynamics.py [relaxation] [capture]   (default: both)
"""
import sys

import numpy as np
from scipy.linalg import expm

from _common import E, CASES, save, log


def relaxation():
    fc = CASES['A']
    N = 30
    p = E.params(fc, s=0.2, N=N)
    L, d, b = E.liouvillian(p)
    th = E.thermal(2.0, N)
    rho0 = np.zeros((d, d), complex)
    rho0[np.arange(N), np.arange(N)] = th            # |g_p> x thermal(2)
    sm = E.simple(p, 2)
    t, nME = E.evolve(L, rho0, 6 / sm['W'], E.number_op(b), h=0.02, n_out=60)
    out = {}
    for lab, c in (('nGRM', E.CSAT), ('nGRM1', 1.0)):
        Q = E.grm_matrix(p, N=N, csat=c)
        out[lab] = [np.arange(N) @ (expm(Q * x) @ th) for x in t]
    nS = sm['n'] + (th @ np.arange(N) - sm['n']) * np.exp(-sm['W'] * t)
    save('dyn_A.dat', ['tms', 'nME', 'nGRM', 'nGRM1', 'nSM'],
         np.column_stack([t / fc['om'] * 1e3, nME, out['nGRM'], out['nGRM1'], nS]))
    log('relaxation: final <n> ME', nME[-1], 'GRM', out['nGRM'][-1])


# protocols: label -> list of stages (Delta/Gamma, s, m); the second stage
# starts at t_switch
PROTOCOLS = {'s1': [(3.0, 1.0, 1)], 's03': [(3.0, 0.3, 1)],
             'staged': [(1.0, 1.0, 1), (3.0, 0.3, 1)], 'LP12': [(3.0, 1.0, 2)]}


def capture(T=10e-6):
    for key, fc in CASES.items():
        nbar = E.KB * T / (E.HBAR * fc['om']) - 0.5
        N = 1600 if key == 'A' else 400        # up to the trap depth (temperature.py)
        tmax, tsw = (150e-3, 40e-3) if key == 'A' else (6e-3, 1.5e-3)
        tw = np.linspace(0, tmax, 301)
        cols, data = ['tms'], [tw * 1e3]
        for lab, stages in PROTOCOLS.items():
            Us = []
            for Dg, s, m in stages:
                p = E.params(fc, s=s, m=m, N=N, Delta_G=Dg)
                Us.append(expm(E.grm_matrix(p, N=N, csat=E.CSAT) * (tw[1] * fc['om'])))
            v = E.thermal(nbar, N)
            nn, p0 = [], []
            for t in tw:
                nn.append(np.arange(N) @ v)
                p0.append(v[0])
                v = (Us[0] if (len(Us) == 1 or t < tsw) else Us[1]) @ v
            cols += ['n_' + lab, 'p0_' + lab]
            data += [nn, p0]
            log('capture', key, lab, 'n0 %.2f' % nbar, 'n(tmax) %.3f' % nn[-1],
                'P(n=0) %.3f' % p0[-1])
        save(f'capture_{key}.dat', cols, np.column_stack(data))


STUDIES = {'relaxation': relaxation, 'capture': capture}

if __name__ == '__main__':
    for name in sys.argv[1:] or STUDIES:
        STUDIES[name]()
