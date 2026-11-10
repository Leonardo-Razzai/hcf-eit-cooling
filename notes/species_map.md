# Other atoms, ions and tweezers on the cooling map

Branch `map-species`. Script `simulation/studies/species.py` (`points`, `sims`), figure
`simulation/plotting/species_plots.py` → `simulation/results/universal/species.{pdf,png}`.
Data: `results/species_points.txt` (one row per point, with source), `results/species_sims.dat`.

## Coordinates

The cold-atom problem depends on ν_r/ω_rec and b_p (and Δ/Γ, s); the tail on
h = α ω_rec/(b_p ν_r), with α = 0.397.

- ω_rec is taken on the **cooling (probe) line**.
- For EIT schemes, b_p is the branching of the excited level into the probe ground level:
  - D1 hyperfine Λ: from 6j symbols for each nuclear spin;
  - J = ½ → ½ Zeeman Λ in ions: 1/3;
  - ¹⁷¹Yb ¹S₀ F = ½ → F′ = ½: 2/3 for a σ probe and π coupling, 1/3 the other way round.
- Schemes without a dark resonance (gray molasses, Raman sideband, narrow-line sideband) are placed at 1/κ, where κ is the estimated number of photons per removed phonon. These points are approximate (marked "approx" in the data file).
- Fibre points use the same geometry as the paper: ν_z = 300 kHz, a = 15 µm (A) and 3.5 µm (B).
  - Alkalis use the 1064 nm lattice.
  - Yb uses the 759 nm magic lattice and Sr the 813 nm one, so ν_r is 0.71× and 0.76× the 1064 nm value.

## Points

| Point | ν/ω_rec | b_p (1/κ) | h | Measured (literature) |
|---|---|---|---|---|
| Ca⁺ EIT, Roos 2000 (1.62 MHz) | 51 | 1/3 | 0.023 | ⟨n⟩ ≈ 0.1 |
| Ca⁺ strings, Lechner 2016 (2.08 MHz) | 66 | 1/3 | 0.018 | ⟨n⟩ = 0.01–0.02 |
| Ca⁺ standing-wave EIT, Xing 2025 (3 MHz) | 95 | 1/3 | 0.013 | ⟨n⟩ ≈ 0.05 in 150 µs |
| ¹⁷¹Yb⁺ tripod, Feng 2020 (4.45 MHz) | 520 | ≈1/3 | 0.002 | ⟨n⟩ = 0.04–0.08 |
| Be⁺ sympathetic, Lin 2013 (2.97 MHz) | 13 | ≈1/3 | 0.09 | ⟨n⟩ = 0.04–0.08 |
| Mg⁺ tripod, Wu 2024 (1.5 MHz) | 15 | ≈1/3 | 0.08 | ⟨n⟩ = 0.26 |
| ⁸⁷Rb EIT tweezer, Qiu 2025/26, radial 45 kHz | 12 | 1/6 | 0.19 | ⟨n_r⟩ ≈ 3.7 while imaging |
| same, axial 4 kHz | 1.1 | 1/6 | 2.2 | not in Lamb–Dicke regime |
| ⁶Li ΛGM tweezer, Blodgett 2023 (700 kHz) | 9.5 | ≈½ | 0.08 | ≈ 11 % ground state while imaging |
| ³⁹K GM tweezer, Ang'ong'a 2022 (180 kHz) | 21 | ≈½ | 0.04 | imaging, no n̄ |
| Cs nanofibre ΛGM and EIT, Pache 2026 (46 kHz) | 22 | ≈½ | 0.04 | loading/OD, no n̄ |
| ⁸⁸Sr 689 nm sideband, Cooper 2018 (211 kHz) | 44 | ≈1 | 0.009 | ground fraction 0.69–0.96 |
| ¹⁷¹Yb Raman sideband, Jenkins 2022 (139 kHz) | 37 | ≈½ | 0.02 | n̄ from sideband ratio |

Fibre points use the best D1 Λ for b_p (worst Λ in brackets):

| Species | Case A: ν/ω_rec | Case A: h | Case B: ν/ω_rec | Case B: h | b_p (best, worst) |
|---|---|---|---|---|---|
| ⁶Li D1 | 0.10 | 4.6 | 0.42 | 1.07 | 8/9 (1/9) |
| ²³Na D1 | 0.29 | 1.66 | 1.23 | 0.39 | 5/6 (1/6) |
| ³⁹K D1 | 0.83 | 0.58 | 3.55 | 0.13 | 5/6 (1/6) |
| ⁸⁷Rb D1 | 1.97 | 0.24 | 8.43 | 0.056 | 5/6 (1/6) |
| ¹³³Cs D1 | 3.81 | 0.14 | 16.3 | 0.032 | 3/4 (1/4) |
| ¹⁷¹Yb 556 nm (³P₁) | 1.35 | 0.44 | 5.79 | 0.10 | 2/3 (1/3) |
| ¹⁷¹Yb 399 nm (¹P₁) | 0.70 | 0.86 | 2.98 | 0.20 | 2/3 |
| ⁸⁷Sr 689 nm (³P₁) | 1.13 | 0.70 | 4.85 | 0.16 | ≈½ |
| ⁸⁷Sr 461 nm (¹P₁) | 0.51 | 1.57 | 2.17 | 0.37 | ≈½ |

## What the comparison says

1. **Ions and tweezers sit at h ≲ 0.1, two to three orders of magnitude stiffer than the fibre in ν/ω_rec.**
   - Every published ion EIT experiment has h = 0.002–0.09. By the minimal model the tail exponent there is 1/h ≳ 11, so the recoil tail is invisible. This is why the tail was never reported, even though the ion and fibre dressed-state physics is identical.
   - The fibre trap is the first EIT setting with h ~ 0.1–1.
2. **The one EIT tweezer experiment comes closest.**
   - Qiu et al. (Rb, Zeeman Λ with b_p = 1/6) is at h = 0.19 radially.
   - When both radial axes are cooled, the effective value is 2h = 0.38, above the variance threshold 1/3.
   - Their axial direction (h = 2.2) is outside both the Lamb–Dicke and the steady-state range.
   - This is the setting in which the heavy tail should be measurable first: a P(n) or ⟨n⟩-versus-P₀ discrepancy on the radial modes.
3. **Choice of species in the fibre (case A geometry):**
   - Cs (h = 0.14) is better than Rb (0.24 with b_p = 5/6).
   - K is marginal (0.58).
   - Na and Li cannot be cooled to a steady state by this mechanism (h > 1).
   - In case B everything except Li is below h = 1/2.
   - Heavier atoms on longer wavelengths gain as m λ²; b_p matters as much as the species.
4. **Yb and Sr.**
   - ¹⁷¹Yb has a clean two-level ground state (F = ½, no hyperfine leakage) and an exact Zeeman Λ.
     - On the 556 nm line it sits at h = 0.44 (case A) and 0.10 (case B), the same as Rb with b_p = ½ in case A.
     - On the 399 nm line it is at h = 0.86 in A, too hot.
   - ⁸⁷Sr on 461 nm is too hot in A (h = 1.6). Its 689 nm line has Γ = 7.5 kHz ≈ ν_r (5.5 kHz in A, 23 kHz in B): that is resolved-sideband cooling, not EIT, and is outside the eliminated model.
     - Point (idea): at the LP₁₁ node the carrier is suppressed for *any* scheme, so narrow-line sideband cooling on 689 nm with a probe in LP₁₁ would be carrier-free too. This is not simulated.

## Master-equation checks (`species.py sims`)

All checks use Δ = 3Γ, eliminated banded solver, N = 250, B = 4, and case-A geometry.

| Case | h | s | P₀ | ⟨n⟩ | W₀ (s⁻¹) |
|---|---|---|---|---|---|
| ¹⁷¹Yb 556, b_p = 2/3, Γ/ν_r = 36 | 0.44 | 0.3 / 1 | 0.982 / 0.962 | 0.062 / 0.143 | 790 / 1720 |
| same, Γ/ν_r = 800 (map default) | 0.44 | 0.3 / 1 | 0.981 / 0.961 | 0.063 / 0.145 | 790 / 1740 |
| ¹⁷¹Yb 556, b_p = 1/3, Γ/ν_r = 36 | 0.88 | 0.3 / 1 | 0.940 / 0.872 | 1.21* / 1.21* | 460 / 1010 |
| ³⁹K D1, b_p = 5/6 | 0.58 | 0.3 / 1 | 0.980 / 0.958 | 0.14 / 0.28 | 1010 / 2310 |
| ³⁹K D1, b_p = ½ | 0.96 | 0.3 / 1 | 0.947 / 0.891 | 1.40* / 1.42* | 640 / 1420 |

W₀ is converted with ν_r = 5.10 kHz for Yb (759 nm lattice) and 7.15 kHz for K. Entries marked * are not converged in N: for h > ½ the mean is carried by the cut-off, and N = 250 is a lower bound. P₀ is reliable.

- **The narrow Yb line does not matter.** At Γ/ν_r = 36 the results agree with the broad-line default within 1–2 %. The reason is that Δ = 3Γ ≈ 110 ν_r still satisfies the elimination conditions. The stored map therefore applies directly to ¹⁷¹Yb on 556 nm.
  - Case-A Yb with the σ probe is close to Rb case A, with b_p = ½ (W₀ = 1.7×10³ against 2.2×10³ s⁻¹, P₀ = 0.96 against 0.955).
  - The two-level ground state removes the Zeeman dark-state problem that costs Rb a factor ≈ 3 (`notes/zeeman_results.md`). Yb may therefore be the better species in a fibre once the Zeeman structure is included.
- **K in case A** with b_p = 5/6 cools (P₀ = 0.96 at s = 1), with a heavier tail (⟨n⟩ = 0.28).
  - Caveat: the ³⁹K 4P₁/₂ hyperfine splitting is 55 MHz ≈ 9Γ, so the neighbour level enters at Δ = 3Γ (stage-2 effect). This is not included here.
- **P₀ stays high up to h ≈ 0.9** (0.87–0.94), consistent with the universality scan: P₀ is set by the core, ⟨n⟩ by the tail.

## Caveats

- Literature b_p and κ for gray-molasses and Raman schemes are estimates; the map is built for a closed EIT Λ.
- Ion points are placed on the map for scale only. Ions have η ≪ 1 and Γ ≫ ν, so they sit far in the "easy" corner, and the model applies.
- Fibre ν_r assumes the same ν_z = 300 kHz for every species. At fixed lattice power, U₀ and ν_z depend on polarisability and mass, so a real design would rescale ν_z. h scales as 1/ν_r ∝ 1/ν_z.
- Sr and Yb have no D1 analogue with hyperfine Λ in the bosonic isotopes. Only the fermions (¹⁷¹Yb, ⁸⁷Sr) have a ground Zeeman Λ. ⁸⁷Sr (I = 9/2) has an unresolved ¹P₁ hyperfine structure, so its b_p is only indicative.

## Open

- Yb 556 nm in the full Zeeman model (it is exact: only two ground sublevels). Check the polarisation requirements of an LP₁₁ σ probe (ties to the vector-mode branch).
- Sr 689 nm resolved sideband at the LP₁₁ node: needs the full (non-eliminated) master equation, since Γ ~ ν_r.
- Cs D1 F′ = 3/4 stage-2 check (neighbour level at 1.17 GHz ≈ 250 Γ: harmless).
