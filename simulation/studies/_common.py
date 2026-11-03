"""Shared paths and helpers for the study scripts."""
import sys
from pathlib import Path

import numpy as np

SIM = Path(__file__).resolve().parents[1]          # simulation/
RESULTS = SIM / 'results'
sys.path.insert(0, str(SIM))

import eitsim as E  # noqa: E402

CASES = E.paper_cases()
NME = 20            # Fock truncation for master-equation steady states


def save(name, cols, rows):
    RESULTS.mkdir(parents=True, exist_ok=True)
    np.savetxt(RESULTS / name, np.array(rows), header=' '.join(cols), comments='')


def load(name):
    return np.genfromtxt(RESULTS / name, names=True)


def log(*a):
    print(*a, flush=True)
