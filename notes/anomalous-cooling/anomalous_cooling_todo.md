# Anomalous (scale-free) Fock statistics in radial EIT cooling — research to-do

Possible separate paper. **Theory note: `anomalous_cooling.tex`** (figures from `make_plots.py`,
run from this folder). Status: **hypothesis**, from the case-A convergence problem
(branch `fix-truncation`, November 2026). Nothing here is established yet.

## Observation that started it

- Case A (a = 15 µm, η_sp = 0.71) needs very large Fock truncations N; case B (η_sp = 0.34) is
  converged at N = 20. The steady state has a long recoil-diffusion tail in P(n).
- Weak-probe N-scan (banded eliminated ME, B = 4, `results/convergence/nconv_weak.dat`):
  - s = 0.3: converged by N ≈ 250 (< 0.1 %).
  - s = 0.01: +0.1 % per 50 states at N = 300 (n_ss ≈ 0.0180).
  - s = 0.1: still +0.65 % per 50 states at N = 300 (n_ss ≈ 0.0324 → extrapolated 0.0329);
    increments shrink ×0.67 per 50 states between N ≈ 225 and 275.
  - W converges much faster than n_ss (W is measured from n = 1).
- Convergence is worst at intermediate s (~0.1), better at both smaller and larger s.

## Minimal model (to be checked)

One two-level system per Fock state: |g_p,n⟩ ↔ |+,n−1⟩ with Raman coupling Ω_R√n, narrow-line
width Γ_n ≈ (Γ/Δ)ν_r, effective saturation **s_n = s·n**.

- Removal rate R(n) ≈ R_sat · s n / (1 + c s n), R_sat ∝ Γ_n/(2κ), κ = 1/b_p, c ≈ 2.25.
- Drift v(n) = −(1 − h) R(n), with **h = κ α η_sp²** (recoil quanta per removed phonon).
- Recoil energy diffusion D(n) ≈ h n R(n).
- Fokker–Planck: P(n) ∝ exp(∫v/D dn)/D → if v/D ≈ −(1−h)/(h n): **P(n) ∝ n^(−1/h)**,
  exponent independent of s at this order. Case A h ≈ 0.40 → P ~ n^−2.5; case B h ≈ 0.09 → n^−11.
- Crude: recoil jumps are large and discrete; the exponent may be modified.

## Hypothesis: crossover in s, transition in h

- **s·n = 1 is a crossover** per Fock state (linear → saturated removal); s sets where the tail starts
  (n* ≈ 1/s) and its weight, which explains the worst convergence at intermediate s.
- **A genuine transition only in h, in the limit N → ∞** (the unbounded ladder is the thermodynamic
  limit, like the photon number at a laser threshold).

Features that would characterise it:
1. Scale-free (power-law) tail P(n) ∝ n^(−1/h) instead of exponential.
2. Successive divergence of moments: h < 1/3 all finite; 1/3 ≤ h < 1/2 Var(n) diverges;
   1/2 ≤ h < 1 normalisable, P₀ finite, ⟨n⟩ = ∞ ("anomalous cooling"); h ≥ 1 no steady state,
   P₀ → 0. Two thresholds: **h = 1/2 (mean energy), h = 1 (existence of a steady state)** —
   stricter than the naive h < 1. Case A (h ≈ 0.40) would be in the "finite ⟨n⟩, divergent Var" window.
3. Finite-size scaling: missing ⟨n⟩ ∝ N^(2 − 1/h); N acts as system size; convergence exponent
   measures h. First rough estimate from s = 0.1: h_eff ≈ 0.33 (nominal 0.40).
4. Universality: exponent depends only on h; tails for different (ν_r/ω_rec, b_p, s) at equal h
   should collapse.
5. Critical slowing down: non-exponential (power-law) relaxation of hot states; Liouvillian gap
   closing with N near threshold; W (from n = 1) stays finite.
6. Rounding by physical cut-offs (act as finite size), each bringing a new parameter:

   | scale | condition | case A, s = 0.1 | parameter |
   |---|---|---|---|
   | saturation | s n = 1 | n ≈ 10 | s |
   | loss of sideband resolution (blue sideband resonant; high-intensity EIT breakdown) | s n ≈ 2(ν_r/Γ_n)² ≈ 2(Δ/Γ)² | n ≈ 180 | Δ/Γ |
   | Bessel non-linearity (capture range) | u₁₁·2x₀√n/a ≈ 1 | n ≈ 470 | x₀/a (core size) |
   | trap anharmonicity shift ~ Γ_n | (3/8)(n ħν_r/U₀)ν_r ≈ Γ_n | n ≈ 1380 | U₀/ħν_r (lattice depth) |
   | trap depth (loss) | n ħν_r = U₀ | n ≈ 1550 | U₀ |

Consequences: W insensitive to the tail; n_ss = core (Γ/4Δ)²(1+h)/(1−h)(1+s) + tail term set by h
(exponent) and s, Δ/Γ, core size, lattice depth (cut-offs); n_ss very sensitive near h → 1/2.
P₀ is the robust figure of merit.

Physical analogy (verify on Crossref before citing): Lévy statistics in subrecoil laser cooling
(VSCPT, Raman cooling; Bardou, Bouchaud, Aspect, Cohen-Tannoudji) — dark state plus random recoil →
broad, non-ergodic statistics; random multiplicative ("Kesten") processes with power-law
stationary tails. Related high-intensity EIT breakdown: Bartolotta et al., arXiv:2411.18818 (2024);
Morillas-Rozas et al., arXiv:2607.22243 (2026, preprint).

## First diagnostics (7 Nov 2026, `results/convergence/tail_mechanism.png`)

Banded eliminated ME, case A (h = κα η_sp² = 0.40), drift and diffusion from Fock states
(`rates_vs_n.dat`), tail at s = 0.1 (`tail_s0p1.dat`):

- **Scale-free tail confirmed at s = 0.1:** P(n) is a clean power law over n = 30–200,
  d ln P / d ln n = −2.31 … −2.39 (log-linear slope varies 4×, so not exponential). Prediction
  −1/h = −2.48. Steepens to −2.9 above n ≈ 200. N = 300 and 400 give identical tails;
  n_ss = 0.03243 → 0.03265 (+0.7 %), W = 335 → 338 s⁻¹.
- **Saturation crossover confirmed:** removal rate grows ∝ n up to n ≈ 1/s and saturates at
  ≈ 2000 s⁻¹ (s = 0.01: linear to n ≈ 50–80, plateau from n ≈ 120; s = 0.1: plateau from n ≈ 10–20;
  s = 1: plateau from n = 1). Diffusion per quantum saturates at ≈ 2800 s⁻¹, i.e. D ∝ n R(n) on
  the plateau, as assumed.
- **Unexpected: no collapse at loss of sideband resolution — a fast high-n cooling channel.**
  Beyond s·n ≈ 20 (≈ 2(Δ/Γ)² = 18, the predicted resolution-loss scale) the removal rate grows again,
  roughly ∝ n, far above the narrow-line bound (s = 1: 4.6×10³ s⁻¹ at n = 30, 4.4×10⁴ s⁻¹ at n = 200;
  s = 0.1: rising from n ≈ 120). This channel cuts off the power-law tail. Hierarchy:
  linear (n < 1/s) → saturated plateau with power-law tail (1/s < n < ~20/s) → fast cooling
  (n > ~20/s). The plateau spans n ≈ 10–200 at s = 0.1 (inside the computed range → slowest
  convergence), only n ≈ 1–20 at s = 1 (tail cut early → fast convergence), and n ≈ 100–2000 at
  s = 0.01 (tail weight tiny, but still drifting +0.1 %/50 states). Mechanism to be understood
  (hot atoms leave the narrow EIT resonance; Doppler-like cooling on the broad line?).
- **Fokker–Planck is only qualitative:** the measured local ratio −v n/D is 0.5–1.2 on the plateau,
  not (1−h)/h = 1.48; the large discrete recoil jumps matter. The tail exponent nevertheless is
  close to −1/h.

## The fast high-n channel: Sisyphus-like (gray-molasses) cooling? (7 Nov 2026)

Figure: `results/convergence/regimes.png` (from `rates_vs_n.dat`). Plotted against s·n₀ the drift
curves **collapse**: s = 0.01 and 0.1 coincide in the resolved and plateau regions; s = 0.1 and 1
coincide at the onset of the fast channel, at s·n ≈ 2(ν_r/Γ_n)² = 18.0 (Γ_n/ν_r = 0.333 for
Δ = 3Γ). Three regimes of one trapped atom:

| regime | condition | removal rate | limited by |
|---|---|---|---|
| resolved EIT (quantum, sideband) | s·n < 1 | R ∝ n | narrow line, unsaturated |
| saturated plateau (power-law tail lives here) | 1 < s·n < 2(ν_r/Γ_n)² | R ≈ const ≈ 2×10³ s⁻¹ | narrow line Γ_n |
| fast, Sisyphus-like | s·n > 2(ν_r/Γ_n)², i.e. Ω_R√n > ν_r | R ∝ n, far above Γ_n bound | optical pumping / light shift |

**Onset, three equivalent forms (derived in the note, Eq. onset).** With tan θ(x) =
(Γ_n/ν_r)√(s/8)·x/x₀ along the LP₁₁ profile and the turning point x_t = 2x₀√n:
Ω_R√n = ν_r ⇔ tan θ(x_t) = 1 (dark state rotated by π/4 over the orbit) ⇔ bright-state light-shift
hill δ_B(x_t) − δ_B(0) = ħν_r ⇔ s·n = 2(ν_r/Γ_n)². In general the hill is ħν_r tan²θ(x_t) ∝ n, so
regime III removes ∝ n quanta per cycle (s = 1, n = 200: ~11 ħν_r vs measured ~22× plateau rate).
The rate model with the calibrated c = 2.25 reproduces regimes I and II of the measured drift.

**Intuitive picture (hypothesis).** When the Raman coupling at level n exceeds the trap frequency
(Ω_R√n > ν_r) the internal state evolves faster than the motion: the atom follows the local dark
state |D(x)⟩ ∝ Ω_c|g_p⟩ − Ω̄_p f_p(x)|g_c⟩, whose mixing angle θ(x) ≈ Ω̄_p G x/(aΩ_c) varies along
the LP₁₁ gradient. Motion couples it non-adiabatically (∝ velocity) to the bright state, which for
blue detuning (Δ > 0) is light-shifted up: the atom climbs a potential hill before optical pumping
returns it to the dark state. That is dark-state Sisyphus / gray-molasses cooling (friction ∝ v², so
energy removal ∝ E ∝ n). It also explains why the GRM (c_sat calibrated for n ≤ 15) underestimates
the ME drift at n₀ = 20–30 in case A (GRM/ME = 0.42 at n₀ = 30, s = 1).

**Literature (verified on Crossref/arXiv; full texts not yet read).** Known separately:
- gray molasses: Grynberg & Courtois, EPL 27, 41 (1994); Weidemüller et al., EPL 27, 109 (1994);
- unified theory of cooling of a trapped atom (resolved/unresolved sideband, gray molasses,
  Λ-enhanced gray molasses): Phatak et al., PRA 110, 043116 (2024), arXiv:2406.19153;
- EIT sub-Doppler friction for free atoms: He, Tengdin, Anderson & Rey, PRA 95, 053403 (2017);
- high-energy breakdown / capture range of EIT cooling for ions: Bartolotta et al., arXiv:2411.18818;
- Λ-enhanced gray molasses in tweezers: Janse van Rensburg et al., PRA 113, 033123 (2026).
Not found: the crossover resolved EIT → plateau → Sisyphus within one trapped atom vs n and s, in a
gradient-only (node) geometry. To check in full texts before claiming novelty.

**Fibre-specific aspects.** The mechanism is general Λ physics; the regime is set by the fibre:
(i) the Sisyphus "lattice" is the transverse mode profile (scale a ≈ 15 µm), not λ/2 — a smooth hill
over the whole radial excursion, large capture range, no corrugation; (ii) the soft trap puts the
crossover inside the thermal population (n̄ ≈ 30 at 10 µK): hot atoms are captured in the
Sisyphus regime and finished in the EIT regime (in tight tweezers/ion traps the crossover is at much
larger n); (iii) the node keeps the low-n regime carrier-free. Consequence: hot-atom capture may be
much faster than the GRM-based capture times of the paper.

## To do

- [x] Analyse the running diagnostics (`studies/convergence.py rates tailshape`) — see above:
      drift v(n), diffusion D(n) vs n₀ at s = 0.01, 0.1, 1 (`rates_vs_n.dat`); P(n) at s = 0.1,
      N = 300/400 (`tail_s0p1.dat`); figure `results/convergence/tail_mechanism.png`.
      Check: R(n) saturation at n ≈ 1/s; D ∝ n R; local ratio −v n/D ≈ (1−h)/h, s-independent;
      log–log tail slope vs −1/h.
- [ ] Write the minimal analytic model properly (R(n), v, D, tail exponent, crossover scales,
      tail contribution to n_ss); improve the Fokker–Planck treatment for large discrete jumps
      (e.g. master equation for P(n) with the exact recoil kernel).
- [ ] Scan h via ν_r/ω_rec and b_p with the banded solver (reduced cases, branch `universal`):
      tail exponents, finite-size-scaling exponents, collapse in h, thresholds h = 1/2 and h = 1.
- [ ] Fast channel — test the Sisyphus picture:
  - [ ] map R(n, s) on a 2D grid (n ≈ 1–300, s = 0.01–10; banded drift, ~6 min per s); collapse vs
        s·n and vs Ω_R√n/ν_r;
  - [ ] vary Δ/Γ: onset should move as 2(ν_r/Γ_n)²/s ≈ 2(Δ/Γ)²/s;
  - [ ] flip the sign of Δ: gray molasses needs blue detuning, so red detuning should turn the
        high-n regime into heating while the resolved regime survives;
  - [ ] scaling of the high-n rate with Ω_c, Δ, Γ (Sisyphus friction, not Γ_n); bright-state
        population and dark-state following at high n;
  - [ ] semiclassical friction/diffusion of an atom moving through the LP₁₁ gradient (linear
        response of the optical Bloch equations) vs ME drift at high n;
  - [ ] consequence: capture from thermal states with banded time evolution (n̄ ≈ 30, N ≈ 400)
        instead of the GRM — true capture time from 10 µK including this channel.
- [ ] Optional: Liouvillian gap vs N near threshold; time-domain relaxation from hot states
      (power-law vs exponential).
- [ ] Role of the cut-offs: include anharmonicity / trap depth / Bessel non-linearity and see how
      they round the transition; dependence on Δ/Γ and core size.
- [ ] Literature check (Crossref): Lévy statistics in subrecoil cooling, Kesten processes,
      dissipative phase transitions in open quantum systems, high-intensity EIT cooling breakdown;
      full texts of Phatak 2024, Bartolotta 2024, He 2017 for the EIT/gray-molasses crossover.
- [ ] Decide: section of the fibre paper, or a separate paper on anomalous Fock statistics in
      laser cooling with soft traps.
