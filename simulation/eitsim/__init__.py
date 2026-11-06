"""Radial EIT cooling of lattice-trapped 87Rb in a hollow-core fibre.

Probe in LP11 (node on the fibre axis), coupling beam in LP01. One radial axis
is treated quantum mechanically; the axial motion is frozen. Internal units:
angular radial trap frequency nu_r = 1.

Modules
  atoms            87Rb D1/D2 data, hyperfine strengths, branching, Zeeman couplings
  fibre            constants, LP mode profiles, lattice trap parameters
  eit              Lambda-system parameters, emission pattern, Fano profile
  operators        oscillator operators, field profiles as operators of x
  master_equation  Liouvillian, steady state, cooling rate, time evolution
  rate_models      simple sideband model and generalised rate model (GRM)
"""
from .atoms import (LINES, GROUND_F, excited_levels, strength, branching, zeeman_coupling,
                    lambda_config, neighbour_shift, reduced_hyperfine,
                    excited_manifold)
from .fibre import (HBAR, KB, M_RB, LAM_D1, GAMMA, U01, U1, amp_ratio, gfac,
                    fibre_case, lp01_power, paper_cases)
from .eit import (GP, GC, EX, NODES, WEIGHTS, ALPHA, params, fano, absorption,
                  excited_weights)
from .master_equation import (liouvillian, liouvillian_eliminated, effective_ops, solve_banded, number_op, expect, factorize,
                              steady_state, relax_rate, solve_me, evolve,
                              fock_drift)
from .rate_models import (CSAT, simple, unified, grm_matrix, grm_steady_state,
                          thermal)
