# Radial EIT cooling of lattice-trapped atoms in a hollow-core fibre

This repository holds the manuscript and the simulation code for a theory paper. The paper
proposes EIT cooling of the *radial* motion of ⁸⁷Rb atoms trapped in a fibre lattice. The coupling
beam travels in the LP₀₁ mode and the probe in the LP₁₁ mode. The probe's node sits on the trap
axis, so the probe couples to the motion only through its gradient and never drives the carrier.

- `paper/manuscript/`: the paper (`main.tex`, Nature-like style in `style/`, figures in `figures/`)
- `paper/background_report/`: a longer technical report with the theory background
- `simulation/eitsim/`: Python package (fibre modes, EIT parameters, operators, master equation,
  rate models)
- `simulation/studies/`: scripts that produce `simulation/results/*.dat`
- `simulation/plotting/`: scripts that turn the results into the paper figures
- `simulation/tests/`: fast regression tests

## Quick start

```bash
make venv          # create .venv and install requirements
make test          # regression tests (~10 s)
make figures       # regenerate the paper figures from the stored results
make paper         # build paper/manuscript/main.pdf (needs TeX Live + latexmk)
make results       # rerun every simulation (~15 min)
```

See `CLAUDE.md` for the project overview, conventions and details.
