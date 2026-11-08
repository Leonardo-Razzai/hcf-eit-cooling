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

Case A is from the converged solver (branch `fix-truncation`); `paper/manuscript/main.tex` still has the
N = 20 values until the author approves the changes.

| Quantity | Case A | Case B |
|---|---|---|
| W at s = 1 (s⁻¹) | 1.6×10³ | 1.2×10⁴ |
| Ground-state rate W₀ at s = 1 (s⁻¹) | 2.2×10³ | not computed |
| n_ss at s = 1 | 0.15 | 0.04 |
| Weak-probe floor n_ss | 0.018 | 0.009 |
| Maximum W at s ≈ 5 (s⁻¹; n_ss = 0.41 in A, above the 0.2 flag) | 3.3×10³ | 1.8×10⁴ |
| Probe power at s = 1 | 8 nW | 0.4 nW |
| Cooling from 10 µK | P₀ = 97.5% after 150 ms (staged Δ = Γ → 3Γ); ⟨n⟩ = 1.45, carried by a 0.4% hot tail | ⟨n⟩ = 0.04 in 2.4 ms |
| Time to P₀ ≥ 0.9 vs k_BT_r/U₀ = 10⁻³–10⁻¹ (s = 1) | 2.6–255 ms, ≈ 2.5 s × k_BT_r/U₀ | 0–7.7 ms, ≈ 80 ms × k_BT_r/U₀ |
| Saturated phonon removal (GRM, s = 1) | 1.9 ms⁻¹ | 12 ms⁻¹ |

- **Simple model:**
  - It is valid at weak probe only if it counts κ = 1/b_p ≈ 2 photons per cycle.
  - At s ~ 1 it overestimates W by up to 38% (A; 43% in B at s = 10) and underestimates n_ss by up to 5×.
    The cause is Raman splitting of the narrow line, since s = 2(Ω_R/Γ_n)², and the recoil tail.
- **Mode impurity:** the rate halves at an LP₀₁ admixture ε ≈ 1.3 η⊥. The simple model's threshold,
  2.8 η⊥, is about 2× too optimistic. A trap offset x_s is equivalent to ε_eff = η⊥ x_s/x₀.
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
  D1 is robust: F′=2 (paper) changes by < 2.5% (retune ≤ 0.05 ν_r); F′=1 with probe on F=2 gives
  W = 2.8–3.1×10³ s⁻¹, n_ss = 0.034–0.043 (A) and 1.6–1.8×10⁴ s⁻¹, 0.018–0.022 (B) after retuning
  δ_e ≈ −0.1 ν_r, i.e. 1.7–1.9× the paper rate in A. D2 F′=2 (A, Δ = 3Γ): W = 1.3–1.8×10³ s⁻¹ (retune
  −0.18 to +0.29 ν_r), limited by the selection n_ss < 0.2, since W has a pole there; W₀ at sign +1 is
  2.2×10³ s⁻¹, within 1% of D1 F′=2. At Δ ≥ 8Γ its W is not converged in N (see Truncation).
  D2 F′=1 (probe on F=1, b_p = 5/6) needs δ_e ≈ −1.3 ν_r, depends on the interference sign by a factor 2
  (A: 1.7–3.3×10³ s⁻¹) and stops cooling for Δ ≥ 8Γ: unreliable without the Zeeman model. Rates with
  n_ss > 0.2 are ill-defined (relaxation from |g_p,1⟩) and flagged.
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
  background_report/     earlier long technical report (theory background, HOM trends,
                         core-size trade-off, NANFs, Raman sideband comparison); pgfplots
                         read data from simulation/results/report/
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
notes/                   research notes; anomalous-cooling/ (theory note + to-do)
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
- `robustness.py [carrier|twophoton|mismatch]` → `carrier_A.dat`, `twophoton_A.dat`,
  `mismatch_{A,B}.dat` (Fig. 4).
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
- `report_scaling.py` → `results/report/*.dat`. It is pure Python and independent of `eitsim`;
  it produces the background-report numbers.

`plotting/make_figures.py [mechanism|saturation|robust|dynamics|temperature|lines]` writes
`paper/manuscript/figures/{mechanism,saturation,robust,dynamics,temperature}.pdf` (paper Figs 1, 3, 4, 5, 6)
and `lines.pdf` (D1 vs D2, stage 2; not yet in the paper).
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

- **Case-A convergence fix in progress (branch `fix-truncation`, from `main`).** The case-A
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
  **Pending:** author's approval of the merge and of the paper numbers; decide W or W0 in the paper; rerun
  `universal` with converged cut-offs; talk (slide 12 and result plots) still has the old case-A numbers.

- **Possible separate paper: anomalous (scale-free) Fock statistics** — hypothesis from the case-A
  convergence problem: recoil diffusion with h = κα·η_sp² may give a power-law tail P(n) ∝ n^(−1/h);
  s·n = 1 is a crossover (tail onset n* ≈ 1/s), a genuine transition only in h as N → ∞ (thresholds
  h = 1/2 for ⟨n⟩, h = 1 for a steady state), rounded by physical cut-offs (sideband resolution,
  Bessel non-linearity, anharmonicity, trap depth). Full plan, evidence and to-do list:
  **`notes/anomalous-cooling/`** (theory note `anomalous_cooling.tex`, to-do `anomalous_cooling_todo.md`). Nothing there is established yet.
- Simulations use an ideal closed Λ system and one radial axis. D1 vs D2: stage 1 (effective Λ)
  and stage 2 (excited hyperfine manifold) done; next the full Zeeman structure (dark states of
  F=2→F′=1, interference signs, tensor light shift of 5P₃/₂, leakage increases κ) and both
  radial axes, e.g. alternating LP₁₁ lobes.
- For submission, move the text into the official Springer Nature template (`sn-jnl.cls`, option
  `sn-nature`).
- The background report still quotes the old purity estimate ε ≲ 2.8 η⊥. The master equation gives
  rate halving at ≈ 1.4 η⊥. Update it if the report is reused.
- Fill the `\note{}` placeholders.

## Git

- **Branch per stage (default workflow, requested by the author).** Do each piece of work on a
  short-lived branch from `main` (e.g. `lines-stage1`), in small, focused commits with imperative
  messages. Merge into `main` with `git merge --no-ff` only when: all tests pass, stored results
  that should not change are reproduced exactly, new results are committed, and the code guide and
  this file are updated. `main` must always have code, results, paper and tests in agreement.
- **One simulation at a time (requested by the author).** Worker processes inside one study are allowed
  (`lines.py`) if they fit in memory. Otherwise never run simulations in parallel: do not
  start a study, check or diagnostic while another is still running (check with
  `ps -eo pid,etime,rss,args | grep python`), and wait for it to finish first. Concurrent jobs
  compete for CPU and memory (11 GB) and make every run slower.
- Do not commit LaTeX build products (see `.gitignore`).
- `figures/*.pdf` and `results/*.dat` *are* tracked: they are inputs or outputs the paper depends on.
- Run `make test` before committing changes to `simulation/`.
