# Anomalous (scale-free) Fock statistics in radial EIT cooling — research to-do

Possible separate paper. Status: **hypothesis**, from the case-A convergence problem
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
   | Bessel non-linearity (capture range) | u₁₁x₀√(2n)/a ≈ 1 | n ≈ 700–900 | x₀/a (core size) |
   | trap anharmonicity shift ~ Γ_n | (3/8)(n ħν_r/U₀)ν_r ≈ Γ_n | n ≈ 600–700 | U₀/ħν_r (lattice depth) |
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
- [ ] Understand the fast high-n cooling channel beyond s·n ≈ 2(Δ/Γ)² (vary Δ/Γ: does the onset
      follow 2(Δ/Γ)²/s?); it truncates the tail and sets where convergence is reached.
- [ ] Optional: Liouvillian gap vs N near threshold; time-domain relaxation from hot states
      (power-law vs exponential).
- [ ] Role of the cut-offs: include anharmonicity / trap depth / Bessel non-linearity and see how
      they round the transition; dependence on Δ/Γ and core size.
- [ ] Literature check (Crossref): Lévy statistics in subrecoil cooling, Kesten processes,
      dissipative phase transitions in open quantum systems, high-intensity EIT cooling breakdown.
- [ ] Decide: section of the fibre paper, or a separate paper on anomalous Fock statistics in
      laser cooling with soft traps.
