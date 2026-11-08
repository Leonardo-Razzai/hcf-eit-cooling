"""Cooling dynamics (paper Fig. 5).

  relaxation  case A, s = 0.2, from <n> = 2: converged banded master equation (BDF2, ~10 min),
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
    """Case A, s = 0.2, from a thermal state with <n> = 2: converged master
    equation (banded, eliminated; N and band from numerics.toml), GRM with
    calibrated and unit c_sat, simple model."""
    fc = CASES['A']
    s = 0.2
    cut = E.me_cutoffs(fc['eta_sp'], s)
    N = cut['N']
    p = E.params(fc, s=s, N=N)
    th = E.thermal(2.0, N)
    sm = E.simple(p, 2)
    t, nME = E.banded_evolve(p, th, 6 / sm['W'], B=cut['band'], h=0.25, n_out=60)
    NG = 200
    thg = E.thermal(2.0, NG)
    out = {}
    for lab, c in (('nGRM', E.CSAT), ('nGRM1', 1.0)):
        Q = E.grm_matrix(E.params(fc, s=s, N=NG), N=NG, csat=c)
        out[lab] = [np.arange(NG) @ (expm(Q * x) @ thg) for x in t]
    nS = sm['n'] + (th @ np.arange(N) - sm['n']) * np.exp(-sm['W'] * t)
    save('dyn_A.dat', ['tms', 'nME', 'nGRM', 'nGRM1', 'nSM'],
         np.column_stack([t / fc['om'] * 1e3, nME, out['nGRM'], out['nGRM1'], nS]))
    log('relaxation (N %d, band %d, %s): final <n> ME' % (N, cut['band'], cut['status']),
        nME[-1], 'GRM', out['nGRM'][-1])


# protocols: label -> list of stages (Delta/Gamma, s, m); the second stage
# starts at t_switch
PROTOCOLS = {'s1': [(3.0, 1.0, 1)], 's03': [(3.0, 0.3, 1)],
             'staged': [(1.0, 1.0, 1), (3.0, 0.3, 1)], 'LP12': [(3.0, 1.0, 2)]}


def capture(T=10e-6):
    for key, fc in CASES.items():
        nbar = E.KB * T / (E.HBAR * fc['om']) - 0.5
        N = E.CFG['grm']['N_case_' + key]      # up to the trap depth (numerics.toml)
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
