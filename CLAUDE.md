# CLAUDE.md — radial EIT cooling in a hollow-core fibre

Guidance for Claude Code working in this repository. Read this first.

## What the project is

A theory/simulation paper (target: Nature Communications / Nature Physics style) on
**radial EIT cooling of lattice-trapped ⁸⁷Rb atoms inside a hollow-core fibre (HCF)**.

- Atoms sit in a 1064 nm standing-wave lattice formed by the fibre's LP₀₁ mode. Axial motion
  can already be ground-state cooled (Raman sideband cooling); radial motion cannot, because
  every guided beam propagates along the axis.
- Idea: EIT cooling with the **coupling beam in LP₀₁** and the **probe in LP₁₁**. The LP₁₁ field has a
  node on the fibre axis, exactly at the trap centre, so the carrier vanishes and the probe couples
  to radial motion only through its gradient, with η⊥ = G₁₁ x₀/a and G₁₁ = 3.49. This is the fibre analogue of an ion
  at a standing-wave node (Cirac et al., PRA 46, 2668, 1992).
- The paper gives the intuitive dressed-state picture and a simple rate model (Morigi et al. 2000).
  It then compares them with master-equation simulations for two fibres:
  - **case A**: antiresonant fibre, core radius a = 15 µm, ν_r = 7.15 kHz;
  - **case B**: bandgap fibre, a = 3.5 µm, ν_r = 30.6 kHz.

  Both use ν_z = 300 kHz, D1 line, Δ = 3Γ. The paper states where the simple model fails, and
  discusses implementations and applications (long-lived spin-wave memories, waveguide QED, …).

### Key results (keep the text consistent with these)

Case A is from the converged solver.

| Quantity | Case A | Case B |
|---|---|---|
| W at s = 1 (s⁻¹) | 1.6×10³ | 1.2×10⁴ |
| Ground-state rate W₀ at s = 1 (s⁻¹) | 2.2×10³ | 1.4×10⁴ |
| n_ss at s = 1 | 0.15 | 0.04 |
| Weak-probe floor n_ss | 0.018 | 0.009 |
| Maximum W at s ≈ 5 (s⁻¹; n_ss = 0.41 in A, above the 0.2 flag) | 3.3×10³ | 1.8×10⁴ |
| Probe power at s = 1 | 8 nW | 0.4 nW |
| Cooling from 10 µK | P₀ = 97.5% after 150 ms (staged Δ = Γ → 3Γ); ⟨n⟩ = 1.45, carried by a 0.4% hot tail | ⟨n⟩ = 0.04 in 2.4 ms |
| Time to P₀ ≥ 0.9 vs k_BT_r/U₀ = 10⁻³–10⁻¹ (s = 1) | 2.6–255 ms, ≈ 2.5 s × k_BT_r/U₀ | 0–7.7 ms, ≈ 80 ms × k_BT_r/U₀ |
| Saturated phonon removal (GRM, s = 1) | 1.9 ms⁻¹ | 12 ms⁻¹ |

- **Simple model:**
  - It is valid at weak probe only if it counts κ = 1/b_p ≈ 2 photons per cycle.
  - At s ~ 1 it overestimates W by up to 40% (38% in A at s ≈ 0.6; 43% in B at s = 10) and underestimates n_ss by up to 5×.
    The cause is Raman splitting of the narrow line, since s = 2(Ω_R/Γ_n)², and the recoil tail.
- **Mode impurity:** the rate halves at an LP₀₁ admixture ε ≈ 1.3 η⊥ ≈ 0.03 (linear interpolation between ε = 0.02 and 0.04; 1.2 η⊥ in log-log). The simple model's threshold,
  2.8 η⊥, is about 2× too optimistic. A trap offset x_s is equivalent to ε_eff = η⊥ x_s/x₀.
  **Caveat (purity.dat, fixed s, converged):** the 1.3 η⊥ figure is for the maximum W (⟨n⟩-based, s ≈ 3), which in case A
  is dominated by the tail. The ground-state rate W₀ halves at ε ≈ 3.6 η⊥ (s = 1) and 2.2 η⊥ (s = 3) in both fibres, and P₀ is
  unchanged (slightly better) up to ε ≈ 3 η⊥. LP₀₁ power < 10⁻³ keeps ≥ 90 % of W₀ at s ≤ 1 in case A. Paper not updated.
- **Saturated phonon-removal rate:** 30–45% of (1−h)Γ_n/2. The generalised rate model (GRM) is
  calibrated with c_sat = 2.25 (17% rms).
- **LP₁₂ vs LP₁₁:** identical rate at equal s; a higher mode only lowers the probe power.
- **D1 vs D2, stage 1 (effective Λ, three-level ME):** only Γ, η_sp, b_p and the leg strengths
  differ. b_p = 5/6 (D1 F′=1 with probe on F=2; D2 F′=1 with probe on F=1) gives W(s=1) = 2.9×10³ s⁻¹
  and n_ss = 0.038 in case A (1.8× faster in W, 1.55× in W₀, 4× colder than b_p = ½), 1.7×10⁴ s⁻¹ / 0.020 in case B.
  b_p = 1/6 does not cool in case A (h > 1). D2 F′=1 (probe F=1) has a weak coupling leg (S = 1/20)
  and a neighbour light shift of −1.3× δ₊, so stage 2 must decide it; D1 F′=1 has −0.11×.
  D2 is otherwise ~5% slower (larger Γ, η_sp). Paper probe powers use S = 1; with S_FF′ they
  are 1.2–2× higher.
- **D1 vs D2, stage 2 (all excited hyperfine levels in the ME, Δ = 3Γ, s = 1, two-photon detuning
  re-optimised, both signs of the Raman interference through the other levels):**
  Case A (branch `universal`), δ_e selected on W₀ (W₀, P₀ at Δ = 3Γ, both signs): D1 F′=2 (paper) W₀ = 2.2×10³ s⁻¹
  within 1% of three levels, P₀ = 0.953–0.955; D1 F′=1 with probe on F=2: 3.2–3.6×10³ s⁻¹, P₀ = 0.974–0.979 (1.5–1.6×
  the paper rate); D2 F′=2: 2.1–2.2×10³ s⁻¹ (within 3%), P₀ = 0.941–0.955; D2 F′=1 (probe on F=1, b_p = 5/6): needs
  δ_e ≈ −1.3 ν_r, 2.0–4.2×10³ s⁻¹ by interference sign, stops cooling near Δ = 8Γ: unreliable without the Zeeman
  model. For P₀ = 0.99 (sign +1): D1 F′=1 at 5Γ, W₀ ≈ 1.9×10³ s⁻¹; D1 F′=2 near 6Γ, ≈ 1.0×10³ s⁻¹. W₀ is meaningful
  only for cold states (P₀ ≳ 0.85). W values in `main` come from selection on W; case-B values are unchanged:
  D1 F′=1 W = 1.6–1.8×10⁴ s⁻¹, n_ss = 0.018–0.022 at δ_e ≈ −0.1 ν_r; D2 F′=2 within ±5% of three levels.
- **Initial temperature:** above n̄₀ ≈ 3 the sideband saturates and the cooling time grows linearly
  with k_BT_r/U₀ (U₀/k_B = 530 µK). Recoil makes hot atoms diffuse in energy (≈ η_sp√n per photon),
  so ⟨n⟩ lags behind P₀. Beyond its calibration range (n ≤ 15) the GRM underestimates the
  master-equation drift (GRM/ME = 0.75, 0.42 at n₀ = 20, 30 for A; 0.83, 0.56 for B), so its
  thermal cooling times are conservative.

## Repository layout

```
paper/
  manuscript/            main paper (Nature-like layout)
    main.tex             text; Results / Discussion / Methods, unnumbered sections
    refs.bib             32 references, no DOI fields (Nature reference style)
    style/
      nature-article.sty layout: 2 columns, Times/Helvetica, "Fig. 1 |" captions, headings
      macros.sty         notation (\etap, \gp, \gc, \SI, \um, ...)
    figures/             PDFs generated by simulation/plotting (tracked: inputs of the build)
  background_report/     long technical report (theory background, HOM trends, core-size trade-off,
                         NANFs, Raman comparison, soft-trap results; pgfplots read results/report/)
                         and core_facts.tex (short: equations and facts, links into report.pdf)
  tools/trim_bib.py      keep only cited .bib entries and strip DOIs
simulation/
  eitsim/                Python package (5 topic modules, see below)
  studies/               scripts that run simulations and write results/
  results/               numerical results (.dat, whitespace columns with a header line)
    report/              data for the background report
  plotting/              style.py (Nature figure style) + make_figures.py
  tests/                 fast regression tests against results/
  docs/                  code_guide.tex: implementation and motivation of every module
                         (make docs); update it when the code changes
  requirements.txt
presentation/            18-slide talk for an atomic-physics audience
  make_assets.py         slide images: setup schematic, level scheme, equations
                         (matplotlib mathtext), result plots, memory dephasing plot
  build_deck.js          pptxgenjs deck (theme, layouts, speaker notes); apply_theme.js
                         writes the theme colours into the .pptx
  assets/                generated PNGs (300 dpi)
  radial_eit_cooling_talk.pptx
notes/                   research notes; anomalous-cooling/ (theory note + to-do);
                         quantum-memory/ (standalone memory note + light scripts, for a separate project)
Makefile                 make test | results | figures | paper | report | slides | clean
```

## Setup (WSL / Linux)

```bash
make venv                       # python3 -m venv .venv; pip install -r simulation/requirements.txt
make test                       # ~30 s, must pass before committing simulation changes
```

LaTeX needs TeX Live with `latexmk`, `pgfplots`/`tikz`, `natbib`, `booktabs`, `caption`,
`hyperref`, `upgreek`, and the `nature` package (provides `naturemag.bst`). On Debian/Ubuntu:

```bash
sudo apt install texlive-latex-extra texlive-pictures texlive-publishers texlive-fonts-recommended latexmk
```

`nature-article.sty` uses the URW Times/Helvetica clones (`utm`/`uhv`) if present, else the psnfss
fonts (`ptm`/`phv`). On Windows/MiKTeX a missing font package made pdflatex hang silently.
If a build stalls, check fonts first.

## Common commands

| Task | Command | Time |
|---|---|---|
| Regression tests | `make test` | ~30 s |
| Rerun all studies | `make results` | ~45 min |
| One study | `cd simulation/studies && python robustness.py twophoton` | — |
| Regenerate figures | `make figures` | seconds |
| Build paper | `make paper`, or in `paper/manuscript`: `pdflatex main; bibtex main; pdflatex main` ×2 | — |
| Build background report | `make report` | — |
| Build code guide | `make docs` | — |
| Build the talk | `make slides` (needs Node.js; `cd presentation && npm install` once) | ~30 s |

The talk reuses `simulation/results`. On WSL, render figures and slide images with the Windows
fonts registered (Arial/Calibri from `/mnt/c/Windows/Fonts`); without them matplotlib falls back to
DejaVu Sans and every regenerated image changes. Slide text uses `x_sub` notation, which `build_deck.js`
turns into real subscripts (`rich()`). ν is set in Cambria, because Calibri's ν looks like a v.
After layout changes, export the slides as images and check them: no text may overlap the
figures, and margins must be at least 0.5 in.

Study scripts (run from `simulation/studies/`):
- `saturation.py` → `sat_A.dat`, `sat_B.dat` (Fig. 3).
- `robustness.py [carrier|twophoton|mismatch|purity]` → `carrier_A.dat`, `twophoton_A.dat`,
  `mismatch_{A,B}.dat` (Fig. 4); `purity.dat` (P₀, ⟨n⟩, W, W₀ against ε at s = 0.3, 1, 3, both cases, ~12 min).
- `dynamics.py [relaxation|capture]` → `dyn_A.dat` (converged master equation, N = 400, ~4 min), `capture_{A,B}.dat` (Fig. 5).
- `calibration.py [steady|drift]` → `grmcheck.dat`, `drift_calib.dat`. These fix `rate_models.CSAT`
  and Methods Table 3.
- `temperature.py [scan|drift]` → `temperature_{A,B}.dat` (metrics, 20 log steps) and
  `temptrace_{A,B}.dat` (⟨n⟩(t), P₀(t) at k_BT_r/U₀ = 1e-3, 4e-3, 1e-2, 4e-2, 1e-1; column
  `n_0p004` = 0.004) (Fig. 6, ~2 min);
  `drift_hi.dat` (master-equation check of the GRM at n₀ = 20, 30, ~25 min).
- `lines.py [effective|multilevel|delta]` → case A only, converged solver (case-B files are kept from the
  full ME at N = 20): `lines_A.dat`, `lines_scan_A.dat` (stage 1, ~40 min); `lines_ml_A.dat` (stage 2 at Δ = 3Γ,
  ~20 min, 4 workers, `OPENBLAS_NUM_THREADS=2`); `lines_delta_A.dat` (stage 2 vs Δ, ~70 min, 2 workers,
  `OPENBLAS_NUM_THREADS=4`).
- `convergence.py` → `results/convergence/`: evidence for `numerics.toml` (elimination, band, N, tail, drift,
  time evolution, D1/D2 checks, W against W0); `plotting/convergence_plots.py` plots it.
- `universal.py [thermal|map|mapband]` → `universal_thermal*.dat` (⟨n⟩, P₀ from thermal states, ~1.5 h), `universal_map.dat` (plane (b_p, ν_r/ω_rec), ~1 h,
  4 workers), `universal_mapband.dat`; `hscan` → `universal_hscan*.dat` (h = 0.1–1.2 by two routes, Δ = 8Γ, N = 400, ~25 min, 2 workers: tail exponent collapses on 1/h); plots by `plotting/universal_plots.py` in `results/universal/`. Older subcommands (steady, drift, capture, robust) still use N = 20.
- `species.py [points|sims]` → `species_points.txt` (literature ions/tweezers and Li–Cs, ¹⁷¹Yb, ⁸⁷Sr in the fibre on (ν_r/ω_rec, b_p)), `species_sims.dat` (Yb 556 nm narrow line and K, ME, ~10 min); plot `plotting/species_plots.py` → `results/universal/species.pdf`; note `notes/species_map.md`.
- `vector_modes.py [fields|carrier|twoaxes|shared|twoaxes_n|rb_sublevels]` → `vector_*.dat`: on-axis E_z of linear LP₁₁/TM₀₁ (ε_z = 1.40 η⊥ in A, 2.90 η⊥ in B; zero for TE₀₁, HE₂₁ and LP₁₁ polarised across its lobes); TE/TM/HE₂₁ splitting and loss in the Zeisberger–Schmidt tube model; 1D ME with E_z on the probe transition (W₀ −6 % A, −31 % B at s = 1; P₀ unchanged); 2D two-channel ME (TE₀₁ components on two Λs): both axes at once, P₀₀ 0.93–0.96, W₀₀ = 0.35–0.46 W₀ (1.85× alternation) only with separate g_c per channel; shared g_c cools only x + y. `rb_sublevels` → `vector_rb_sublevels.dat` (⁸⁷Rb D1 realisation: target |2,1⟩ via F′=2, two probe frequencies, 6× coupling imbalance); figures `plotting/te01_plots.py` → `results/vector/` (report Sec. TE₀₁). ¹⁷¹Yb (F=½) cannot do the TE₀₁ scheme. Note `notes/vector-modes/vector_modes_theory.md`.
- `zeeman2d.py [check|check_field|smoke|cost|run|list|collect]` → 2D Zeeman master equation for ⁸⁷Rb (`eitsim/zeeman2d.py`): `lp11_hf`, `lp11_zl` (LP₁₁ switched), `te01_bz`, `te01_bx`. **Cluster only** (N = 6: ~1 h per run, N = 8: ~6 h): scan lists in `simulation/cluster/scans/`, sbatch scripts and instructions in `simulation/cluster/` (README.md); results → `results/zeeman2d/runs/*.json`, `collect` → `results/zeeman2d.dat`. Submit with `simulation/cluster/submit_scan.sh` (array + dependent `analyse.sbatch`: collect, `plotting/zeeman2d_plots.py` → `results/zeeman2d/summary.txt` + figures). From here: `simulation/cluster/remote.sh status|summary|fetch` (needs an open ssh ControlMaster to `ophfe1`; the bastion asks for a password; if `ssh -o BatchMode=yes ophfe1 true` fails, ask the author to log in once). `check`/`check_field` (reduction to the 1D Zeeman solver) and `smoke` are cheap.
- `ramp.py [drift|chain|optimise[=smax]|twostage]` → `ramp_*.dat`: s-ramp step 3 (chain calibrated on ME drift, optimal control); the optimum is two stages (s_max, then s ≈ 1), 2× faster to P₀ = 0.9 from n̄₀ = 30 on the chain (case A, b_p = ½). `validate` (step 4, banded ME, piecewise s via `E.banded_evolve_piecewise`): 5.8, 7.6, 10.2 ms for s_max = 10, 5, 3 against 15.6 ms at constant s = 2 (2.7, 2.1, 1.5×).
- `two_axes.py [validate|steady|thermal|exact|switch|switch_n]` → `two_axes_*.dat`: 2D rate chain on P(n_x, n_y) (P₀ reliable, ⟨n⟩ an upper bound), and the exact 2D ME at N = 6 per axis; `switch` gives P₀₀ and the 2D rate W₀₀ against the stage length τ of the alternation (`two_axes_switch.dat`, ~7 min): τ ≳ 5/Γ_n needed, W₀₀ peaks at ≈ 0.2 W₀ for τΓ_n ≈ 5–15; `switch_n` confirms within 1 % at N = 8, 10.
- `fast.py [sign|delta]` → `fast_sign.dat`, `fast_delta.dat`: Fock drift at Δ = −3Γ (params accepts Δ < 0; narrow peak at δ = −ν_r, probe offset d2 ≈ +0.8) and at Δ = 2, 5Γ. The fast channel heats at Δ < 0 (capture range below n ≈ 8); its onset and the dip (n ≈ 0.5 n_c) move with n_c.
- `report_scaling.py` → `results/report/*.dat`. It is pure Python and independent of `eitsim`;
  it produces the background-report numbers.

`plotting/make_figures.py [mechanism|saturation|robust|dynamics|temperature|lines]` writes
`paper/manuscript/figures/{mechanism,saturation,robust,dynamics,temperature}.pdf` (paper Figs 1, 3, 4, 5, 6)
and `lines.pdf` (D1 vs D2, stage 2, W₀ and P₀; not yet in the paper).
Fig. 2, the level scheme, is TikZ inside `main.tex`.

## Numerical settings: `simulation/numerics.toml`

All numerical choices live in `simulation/numerics.toml` (loaded by `eitsim/config.py` as `E.CFG`),
each with its justification: steady-state cutoff rules, legacy full-ME N, `function_of_x` threshold,
emission quadrature, GRM `csat` and Fock spaces, `n_rate_max`. **Steady states should use
`E.steady(fc, s, **kw)`**: it takes N and the band from the first `[[steady_state.rule]]` matching the
actual η_sp and s (`E.me_cutoffs`) and calls `solve_banded`; parameters no rule covers raise an
error (run a convergence check, then add a rule). Rules have a `status` (`validated` /
`provisional`). Change a value only with convergence evidence, then rerun the affected studies.

## Simulation package `eitsim`

Internal units: **angular radial trap frequency ν_r = 1**. Physical rates = value × `fc['om']`.

| Module | Topic |
|---|---|
| `atoms.py` | ⁸⁷Rb D1/D2 data (Steck rev. 2.3.4; D1 Γ kept at 5.746 MHz), `strength` S_FF′, `branching`, `zeeman_coupling` (6j/CG via sympy), `lambda_config(line, Fe, gp)`, `neighbour_shift` |
| `fibre.py` | SI constants, ⁸⁷Rb D1 data, LP mode constants (`U01`, `U1`), `amp_ratio`, `gfac`, `fibre_case(a, line=)` → trap parameters + line data, `lp01_power`, `paper_cases(line)` |
| `eit.py` | Λ-system parameters `params(fc, s, Delta_G, m, eps, xs, d2, nu_rel, bp, N)` → dimensionless dict; emission-pattern quadrature `NODES/WEIGHTS/ALPHA`; Fano profile `fano` |
| `operators.py` | oscillator operators, functions of x̂ via the eigenbasis of the truncated b+b†, probe/coupling profiles and recoil operators |
| `master_equation.py` | `liouvillian`, `factorize`, `steady_state`, `relax_rate`, `solve_me`, BDF2 `evolve`, `fock_drift` |
| `zeeman.py` | all ground sublevels, eliminated and banded: `structure_lambda`, `structure_zeeman_lambda`, `solve` |
| `rate_models.py` | simple model `simple`, `unified`, generalised rate model `grm_matrix`, `grm_steady_state`, `thermal`, `CSAT = 2.25` |

Conventions and gotchas:
- **Vectorisation.** ρ is vectorised column-stacked, so vec(AρB) = (Bᵀ⊗A)vec(ρ).
- **Steady state.** The first row of L is replaced by the trace condition. Factorise with
  `permc_spec='MMD_AT_PLUS_A'`, which is about 4× faster than the default.
- **Rates.** `solve_banded` returns W (integrated ⟨n⟩ relaxation from |g_p,1⟩) and W0 (the same for P₀). W gets a
  negative contribution while the recoil tail fills: it depends on N, has a pole and changes sign near n_ss ≈ 0.2.
  W0 is stable (`results/convergence/rate_p0.dat`), and above W by ≈ 1.2× even at weak probe. Studies still report W.
- **Cooling rate.** `relax_rate` returns the inverse *integrated* relaxation time from |g_p,1⟩.
  Do **not** use the slowest Liouvillian eigenvalue. It picks up slow modes at the truncation edge
  and drifts with N, because saturated high-n states cool linearly in time.
- **Truncation.** N = 20 Fock states is **not** enough for case A (η_sp ≈ 0.71) at s ≲ 1: a smooth,
  roughly exponential recoil-diffusion tail (P(n) falls ≈ ×0.89 per quantum; 0.3 % of atoms above
  n = 10 carry a third of ⟨n⟩) is cut off. Case A, s = 1: n_ss = 0.133 / 0.147 / 0.152 / 0.154 and
  W = 1319 / 1446 / 1546 / 1605 s⁻¹ for N = 20 / 30 / 40 / 60 (n_ss converged ≈ 0.155, W still
 rising). Case A at s = 5 and case B are converged at N = 20. Always check N-convergence when
  η_sp ≳ 0.5. The tail also depends on b_p and Δ, so the `numerics.toml` rules do too; D2 F′=2 at Δ ≥ 8Γ is
  not converged (effective b_p falls, h ≈ 0.57 for g_p = 1 at 12Γ). Use N ≥ 30 for dynamics. Thermal GRM dynamics need N up to the trap depth, U₀/ħν_r = 1552 (A) / 362 (B):
  the code uses N = 1600 / 400. Recoil diffusion carries a few atoms to high n, and a smaller N
  acts as an artificial sink (N = 300 gave ⟨n⟩ = 0.24 instead of 1.45 for case A from 10 µK).
- **Recoil.** The master equation scatters κ = 1/b_p photons per cycle automatically, since decays
  into g_c are re-pumped. Any rate model must include this.
- **Detuning convention.** `params` keeps lasers at the design values when `nu_rel ≠ 1`
  (trap-frequency mismatch). Results are rescaled to the actual ν_r.
- **Multi-level ME.** `params(..., levels=E.excited_manifold(line, Fe, gp), sign=±1)`; K = 1
  reproduces the three-level model. LU fill dominates (≈ 20 s per solve for D2, K = 4, N = 20).
- **Determinism.** The emission quadrature uses a fixed RNG seed, so results are reproducible to
  ~1e-8. The tests rely on this.

## Paper conventions (requested by the author — keep them)

- **Layout.** Nature-family style:
  - The introduction has no heading.
  - **Results** (subsections), **Discussion** (challenges, two implementations, applications,
    conclusion) and **Methods** (derivations, numerics, GRM, experimental estimates).
  - Then Data/Code availability, References, Acknowledgements, Author contributions,
    Competing interests.
- **Main text vs Methods.** The main text keeps only what is needed. Details go to Methods.
  The standard EIT-cooling derivations are *not* repeated; they cite Morigi et al. 2000.
- **Abstract.** Short (~150 words) and to the point, with the main numbers.
- **Citations.**
  - Superscript numbers only (`\cite{key}`, no `~` before it).
  - Use `ref.~\citenum{key}` when a reference is a noun.
  - Never write full citations or DOIs in the text.
- **References.**
  - At most ~30 (hard limit 50).
  - Each statement gets one main reference, at most two.
  - Models and implementations must be supported by a citation.
  - **Verify every new reference on Crossref** (api.crossref.org) before adding it.
  - Run `tools/trim_bib.py` to drop unused entries and DOIs.
- **Figures.**
  - Nature style (`plotting/style.py`): Arial 7 pt, 183 mm wide, bold lowercase panel letters,
    Okabe–Ito colours.
  - **No text may overlap data.** Put legends and labels outside the axes, or in regions kept
    empty by fixed axis limits.
  - Render and inspect every figure after changes.
- **Placeholders.** `\note{...}` (red) marks things the author must fill in: authors, affiliations,
  repository DOI, acknowledgements, contributions, competing interests.
- **Workflow for numbers.** Change physics → rerun the affected study → `make figures` → update
  numbers in `main.tex` (and in this file's key results) → `make paper` and check the PDF.
  Never hand-edit files in `results/`.

## Open items / next steps

- **Zeeman structure (one axis, `notes/zeeman_results.md`, `eitsim/zeeman.py`, `studies/zeeman_scan.py`).** Dark states at every n force remixing; best W₀ about 3× below the closed Λ. The effective 2D chain (`two_axes_zeeman.py`) passes the closed-Λ benchmark but fails the P(n) test for Zeeman schemes: do not quote its Zeeman P₀₀. **The full 2D Zeeman ME is written (`eitsim/zeeman2d.py`, validated against the 1D solver) but not run: the author runs it on the cluster (`simulation/cluster/README.md`, scan lists `zeeman2d_{test,main,convergence}.txt`, case A only: 108 + 6 runs). Do not run the scans locally.**
- **Prior work:** carrier-free EIT at a node is Xing & Mehta, PR Applied 24, 014034 (2025), and the heavy-tail law is Dechant et al. (2015, 2016). Do not present either as new; see `notes/paper_draft_ideas.md`.
- **Both radial axes (theory and chain done; `notes/anomalous-cooling/two_axes_theory.md`).** One probe cools one motional mode; recoil heats both axes, so a cold state needs h < 1/2. h = 0.188 ζ/b_p with ζ = k_L a ω_rec/ν_z: the limits on a scale with ν_z. Numbers in the key-results table are for ONE axis.
- **Case-A convergence fix (done).** The case-A
  master-equation results at N = 20 (key results table, Figs 3, 4, 5a, GRM calibration) are not
  converged (see Truncation). Brute force is too slow (one N = 90 solve > 1 h, 7 GB). Plan:
  **adiabatic elimination of the excited state** (effective-operator formalism, Reiter & Sørensen,
  PRA 85, 032111, 2012) as an option of the ME solver. Valid because |Δ + iΓ/2| ≫ Ω_c, ν_r, Γ_n,
  Ω⊥√n (Ω_c/Δ ≈ 0.04, ν_r/Δ ≈ 4×10⁻⁴); it keeps the EIT interference exactly (effective Hamiltonian
  and jump operators are built from the bright combination Ω̄_p f_p⟨g_p| + Ω_c⟨g_c|, so the dark
  state stays dark), reproduces the Fano profile for |δ| ≪ Δ (narrow peak, zero at δ = 0, blue-sideband
  suppression (Γ/4Δ)²) and removes the stiffness (fastest rate Γ_n instead of Γ ≈ 800 ν_r) and a third
  of the dimension ((2N)² instead of (3N)²). Steps: implement → validate against the full ME at
  N = 20, 40 (A and B, several s; < 1 % in n_ss and W) → converge case A in N (60, 90, 120) → rerun
  affected case-A studies and figures (results only; paper text changes need the author's approval)
  → merge into `main`. Then redo the `universal` studies with converged N.
  **Status:** elimination validated (< 0.5 % vs full ME, N = 20 and 40) but gives no speed-up for
  direct solves (recoil operator nearly dense, nnz ∝ N⁴). Added `solve_banded` (eliminated ME,
  coherences |n − m| ≤ B = 6, dense LU): exact for B ≥ N, 2×10⁻⁸ error at B = 6, N = 160 in 34 s.
  Converged case A, s = 1: n_ss = 0.1545, W = 1614 s⁻¹ (paper N = 20: 0.133, 1319). Checks and plots:
  `studies/convergence.py` → `results/convergence/`, `plotting/convergence_plots.py`.
  Then: time evolution validated, Fig. 5a data, D1/D2 rerun, figures, code guide done. W0 added (Rates, above).
  Merged into `main`, with the paper numbers updated. **Pending:** rerun `universal` with converged cut-offs;
  talk (slide 12 and result plots) still has the old case-A numbers.

- **Possible separate paper: anomalous (scale-free) Fock statistics** — hypothesis from the case-A
  convergence problem: recoil diffusion with h = κα·η_sp² may give a power-law tail P(n) ∝ n^(−1/h);
  s·n = 1 is a crossover (tail onset n* ≈ 1/s), a genuine transition only in h as N → ∞ (thresholds
  h = 1/2 for ⟨n⟩, h = 1 for a steady state), rounded by physical cut-offs (sideband resolution,
  Bessel non-linearity, anharmonicity, trap depth). Full plan, evidence and to-do list:
  **`notes/anomalous-cooling/`** (theory note `anomalous_cooling.tex`, to-do `anomalous_cooling_todo.md`). Nothing there is established yet.
- **Optimal s ramp (done, steps 1–4).** Two stages, s_max then s = 1, are 1.5–2.7× faster to P₀ = 0.9 (master equation). Theory and results in `notes/anomalous-cooling/s_ramp_theory.md`. Open: a recipe across the map.
- Simulations use an ideal closed Λ system and one radial axis. D1 vs D2: stage 1 (effective Λ)
  and stage 2 (excited hyperfine manifold) done; next the full Zeeman structure (dark states of
  F=2→F′=1, interference signs, tensor light shift of 5P₃/₂, leakage increases κ) and both
  radial axes, e.g. alternating LP₁₁ lobes.
- For submission, move the text into the official Springer Nature template (`sn-jnl.cls`, option
  `sn-nature`).
- The background report is updated (purity, ramp, switching, fast channel, memory OD); the paper is not.
- Fill the `\note{}` placeholders.

## Git

- **Branch per stage (default workflow, requested by the author).** Do each piece of work on a
  short-lived branch from `main` (e.g. `lines-stage1`), in small, focused commits with imperative
  messages. Merge into `main` with `git merge --no-ff` only when: all tests pass, stored results
  that should not change are reproduced exactly, new results are committed, and the code guide and
  this file are updated. `main` must always have code, results, paper and tests in agreement.
- **Delete the stage branch after merging (requested by the author).** Right after a `--no-ff` merge into
  `main`, run `git branch -d <branch>` (the safe form: it refuses unmerged branches), and delete its remote copy
  if it was pushed. Nothing is lost: the merge commit keeps the branch tip as its second parent, so a branch can
  be recreated with `git branch <name> <merge-commit>^2` (find it with `git log --merges --oneline main`). Keep
  only `main`, `universal` and branches with unmerged work.
- **Paper and talk need explicit approval.** After producing new results, always ask the author
  whether they should go into the paper and/or the presentation; never add them unasked.
- **Long-lived branch `universal`** (from `main`): reduced, fibre-specific results (universal design
  map, generalised robustness, Rb case-A imperfections) and the eventual paper/talk rework. Commit
  there in focused steps; merge into `main` only when the author decides.
- **One simulation at a time (requested by the author).** Worker processes inside one study are allowed
  (`lines.py`) if they fit in memory. Otherwise never run simulations in parallel: do not
  start a study, check or diagnostic while another is still running (check with
  `ps -eo pid,etime,rss,args | grep python`), and wait for it to finish first. Concurrent jobs
  compete for CPU and memory (11 GB) and make every run slower.
- Do not commit LaTeX build products (see `.gitignore`).
- `figures/*.pdf` and `results/*.dat` *are* tracked: they are inputs or outputs the paper depends on.
- Run `make test` before committing changes to `simulation/`.
