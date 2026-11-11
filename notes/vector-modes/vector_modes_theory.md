# Real fibre modes: polarisation at the node, beating along z, one vector mode for both axes

Branch `fibre-vector-modes`. Study `simulation/studies/vector_modes.py`
(`fields`, `carrier`, `twoaxes`); results `vector_fields.dat`, `vector_carrier.dat`, `vector_twoaxes.dat`.

The paper and all simulations so far treat LP₁₁ as a scalar field: transverse, linearly polarised,
∝ J₁(u₁r/a) cos φ. In a real fibre, LP₁₁ is a quasi-degenerate group of four vector modes. Every guided
field also has a longitudinal component E_z. This note covers three questions:
1. Does the node survive?
2. What does the near-degeneracy do along the atom column?
3. Can one vector mode cool both radial axes without switching?

## 1. The LP₁₁ group near the axis

The weakly guiding group is {TE₀₁, TM₀₁, HE₂₁ᵉ, HE₂₁ᵒ}. To first order in r, the transverse fields near the axis are:

| Mode | E_t ∝ | ∇·E_t | E_z on the axis |
|---|---|---|---|
| TE₀₁ | (−y, x), azimuthal | 0 | 0 (TE: E_z ≡ 0) |
| TM₀₁ | (x, y), radial | 2 | ≠ 0 (E_z ∝ J₀) |
| HE₂₁ᵉ | (x, −y) | 0 | 0 (E_z ∝ r² cos 2φ) |
| HE₂₁ᵒ | (y, x) | 0 | 0 |

In the air core, Maxwell's equations give E_z = (i/β) ∇·E_t exactly. The two kinds of linear LP₁₁ therefore differ:
- **Polarisation along the lobe axis**, x x̂ = ½(TM₀₁ + HE₂₁ᵉ): ∇·E_t ≠ 0, so **E_z ≠ 0 on the axis**.
- **Polarisation across the lobe axis**, y x̂ = ½(HE₂₁ᵒ − TE₀₁): ∇·E_t = 0, so **E_z = 0 on the axis**. The node is exact in all three components.

**E_z as a carrier.** With E_x = A J₁(u₁r/a) cos φ, normalised to the LP₀₁ axis field at equal power (A = amp_ratio = 1.82), the on-axis field is

  ε_z = A u₁/(2ka),

in the same units as the LP₀₁ admixture ε of the purity study. Relative to the probe gradient η⊥ = G₁₁x₀/a,

  ε_z/η⊥ = A u₁/(2 G₁₁ η_sp) ≈ 1.0/η_sp.

The ratio does not depend on a except through η_sp.

| | ε_z | η⊥ | ε_z/η⊥ |
|---|---|---|---|
| Case A (a = 15 µm) | 0.029 | 0.021 | 1.40 |
| Case B (a = 3.5 µm) | 0.126 | 0.044 | 2.90 |

**Comparison with the purity results.** These are at the scale where a same-polarisation LP₀₁ admixture matters:
- the W₀ of case A halves at 3.6 η⊥ (s = 1);
- the maximum W halves at 1.3 η⊥.

**Where E_z acts.** E_z is polarised along the fibre axis, so its effect depends on the quantisation axis:
- **Quantisation along the probe polarisation** (π probe in the transverse plane), or **along the fibre axis** (σ probe): E_z drives *other* transitions from g_p, not the Λ probe transition. It is then a weak off-resonant scatterer, at most ε_z² × the probe scattering rate, and is not EIT-suppressed.
- **Quantisation transverse but not along the probe polarisation** (e.g. B ∥ y with an x̂ probe): the probe is a σ⁺/σ⁻ mixture with a ẑ part, so E_z enters the probe transition itself. This is the worst case simulated in `carrier`.
  - E_z is in quadrature (factor i) with the transverse field, so the profile becomes A x + iε_z. The node is *filled* (minimum intensity ε_z²), not shifted. A real ε instead shifts the node to x = −ε/A (the purity study and the trap-offset equivalence).

**Design rules.**
- Polarise the LP₁₁ probe across its lobe axis, or use TE₀₁ or HE₂₁.
- Put the quantisation axis along the probe polarisation, or along the fibre.
- The coupling beam (HE₁₁) has E_z ∝ (i/β)∂_x J₀ ∝ x: zero on the axis, a weak gradient on the coupling field only. It does not create a probe carrier.

## 2. Splitting and beating along z

**Tube model** (Zeisberger & Schmidt, Sci. Rep. 7, 11761 (2017), Eqs. 25–26). For an antiresonant tube with wall phase φ = k d √(n²−1), up to (kR)⁻⁴:

  Re n_eff = 1 − j²/2(kR)² − j²/(kR)³ · cot φ/√(ε−1) · {1 (TE), ε (TM), (ε+1)/2 (HE)}.

Here j = j₁,₁ = 3.832 for all four LP₁₁ members and ε = n_glass².
- **At antiresonance** (cot φ = 0) TE₀₁ and TM₀₁ are degenerate to this order. HE₂₁ − TE₀₁ = j²/(kR)⁴ comes from the σ⁴ term.
- **Away from antiresonance** the group splits linearly in cot φ.
- **Loss at any phase**, relative to HE₁₁: TE₀₁ 1.49, HE₂₁ 4.05, TM₀₁ 6.6 (ε = 2.1). **TE₀₁ is the lowest-loss member**, in agreement with the hollow-core bandgap fibres that guide TE₀₁ alone over 50–100 nm.

Beat lengths L = λ/Δn at 795 nm (`vector_fields.dat`):

| Phase offset from antiresonance | Case A: TM–TE | Case A: HE₂₁–TE | Case B (tube scaling only): TM–TE | Case B: HE₂₁–TE |
|---|---|---|---|---|
| 0 | ∞ (σ⁴) | 10.7 m | ∞ | 3.2 cm |
| 0.05π | 0.54 m | 0.99 m | 6.9 mm | 9.6 mm |
| 0.1π | 0.26 m | 0.50 m | 3.4 mm | 5.5 mm |
| 0.2π | 0.12 m | 0.23 m | 1.5 mm | 2.7 mm |
| 0.4π | 2.8 cm | 5.6 cm | 0.36 mm | 0.69 mm |

- Case B is a bandgap fibre, so the tube model gives only the (kR)⁻³ scaling. The splitting is larger there; that is what makes TE₀₁-only bandgap fibres possible.
- Bends, ellipticity and twist add splitting within the HE₂₁ pair and between the members. Only a mode solver for the real structure, or a measurement, can fix these numbers.

**What the beating does.** A linear LP₁₁ launched at z = 0 is a sum of two vector modes. The relative phase Δβ z rotates the local pattern along the atom column. For example:

  ½(TM₀₁ + HE₂₁ᵉ e^{iΔβz}) goes from x x̂ (z = 0) to y ŷ (Δβ z = π),

with elliptical mixtures between. With a polarisation-selected Λ, an atom at z is cooled along an axis set by z, at a rate ∝ cos²(Δβz/2). The lobe-parallel launch also carries a z-dependent E_z.
- **Case A:** at most a few 10⁻² of a beat over a ~1 cm column near antiresonance; the scalar picture holds.
- **Case B:** the column can span several beat lengths. Every site is cooled along a different axis, and some sites see the probe in the wrong polarisation.
- **TE₀₁ is non-degenerate:** a pure TE₀₁ probe does not beat within the group (only with the HE₂₁ pair, through imperfections), and keeps E_z = 0.

## 3. One vector mode for both axes

**Why it can work.** TE₀₁ has E_x ∝ −y and E_y ∝ x: each transverse polarisation component carries the gradient along the *other* axis. If the two components drive two Λ channels with **distinguishable Raman final states**, one channel has coupling operator ∝ y and the other ∝ x, and both axes are cooled at the same time, without switching.

What "distinguishable" requires: a different excited sublevel is **not** enough. A separate coupling ground state g_c per channel is needed (simulation below). EIT cooling is the Raman process g_p,n → g_c,n−1. If both channels end in the same g_c, their coherent (Hamiltonian) parts add as (f_x + f_y)†-type operators, and only the diagonal mode x + y is cooled. Only the spontaneous jumps stay separate, through the distinct excited levels.

**Two realisations:**
- **B ∥ x:** E_x is π (∝ y → cools y); E_y drives σ⁺ and σ⁻ (∝ x → cools x). The coupling beam (LP₀₁) needs both π and σ components, e.g. linear at 45° in the transverse plane.
- **B ∥ z (fibre axis):** the σ± components of TE₀₁ are ∝ x ∓ iy = r e^{∓iφ}. Each couples to one *circular* motional mode, since x + iy ∝ a_L + a_R†: σ⁻ removes L quanta, σ⁺ removes R quanta. Together they cool both, and L and R span the plane.
  - The coupling beam is guided, hence σ±. To keep the final states apart it should be **circular**, e.g. pure σ⁺:
    - channel 1: g_p(m) →σ⁺→ e(m+1) ←σ⁺← g_c(m);
    - channel 2: g_p(m) →σ⁻→ e(m−1) ←σ⁺← g_c(m−2).
  - A linearly polarised coupling would also connect both e(m±1) to g_c(m): a shared final state.
  - The two Raman resonances differ by the Zeeman splitting of g_c(m) and g_c(m−2), 2g_F μ_B B. They need either a two-frequency probe (both frequencies in TE₀₁) or a field small compared with Γ_n but large enough to split the Zeeman dark states. This is to be checked in the Zeeman model.

**Why one transition is not enough.** Both components on one transition (or a circularly polarised LP₁₁ / OAM beam on one transition) give a single coupling operator ∝ x + iy. That is one motional mode: L cools while R is only heated by recoil. This is the `coh` scheme in the simulation.

**Recoil.** Recoil is unchanged: photons from both channels heat both axes. The two-axis condition h → 2h of the minimal model still holds. What changes is the *rate*: alternation loses ≈ 5× (W₀₀ ≈ 0.2 W₀ at the optimal stage length τΓ_n ≈ 5–15), because each axis is cooled only half the time and every switch costs ~1/Γ_n of re-equilibration of the dark state.

**Minimal model** (`twoaxes`). Eliminated 2D master equation, N = 6 per axis. Ground states g_p plus one shared g_c (`te01`) or two (`te01s`); one excited level per channel. Channel profiles kron(f_p, 1) and kron(1, f_p); coupling kron(f_c, f_c); the same dipole emission pattern as `two_axes.exact_2d`. Each channel is at saturation s. At equal total power, TE₀₁ has s/2 per component compared with a linear LP₁₁ at s, since each component has amplitude 1/√2.

**Results** (`vector_twoaxes.dat`, `vector_shared.dat`). N = 6 per axis, Δ = 3Γ. W₀₀ = P₀₀ss / ∫(P₀₀ss − P₀₀)dt from |1,1⟩.

| Case | s per channel | Alternation, τ = 32 (P₀₀, W₀₀/W₀) | Circular profile, one channel | Two channels, shared g_c (also with g_c light shift retuned) | Two channels, separate g_c |
|---|---|---|---|---|---|
| B (ν/ω_rec = 8.4, b_p = ½) | 0.5 | 0.954, 0.25 | 0.24, no steady cooling | 0.16 / 0.28 | **0.960, 0.46** |
| B | 1 | 0.931, 0.22 | 0.17 | 0.15 / 0.27 | **0.937, 0.41** |
| A (ν/ω_rec = 2.0, b_p = 5/6) | 0.5 | 0.947, 0.22 | 0.23 | 0.13 / 0.30 | **0.953, 0.41** |
| A | 1 | 0.925, 0.19 | 0.21 | 0.13 / 0.29 | **0.931, 0.35** |

- **Separate final states: one vector mode cools both axes at once.** P₀₀ is slightly above alternation (0.93–0.96), and W₀₀ is 1.8–1.9× larger at equal s per channel.
- **Converged in N:** at N = 8, P₀₀ = 0.9367 / 0.9294 and W₀₀ = 0.0295 / 0.0265 ν_r (B / A, s = 1), within 0.2 % and 2 % of N = 6 (`vector_twoaxes_n.dat`).
- **The circular-profile single channel cools one circular mode only.** P₀₀ ≈ 0.2, and ⟨n_x⟩ = ⟨n_y⟩ ≈ 1 at the truncation N = 6.
- **The shared-g_c scheme fails** whether or not the doubled g_c light shift is compensated (coupling Rabi frequency /√2 per channel), so it is not a detuning artefact.
- **Why W₀₀ is 0.4 W₀, not ≈ W₀.** From |1,1⟩, two independent channels at rate W give P₀₀ = (1 − e^{−Wt})² and W₀₀ = W/1.5 = 0.67 W₀. Alternation, with each axis cooled half the time, gives 0.33 W₀. Both measured values are 55–70 % of these ideals; the remainder is the doubled recoil (h → 2h) and, for alternation, the switching transient.
  - This partly answers the open question of why alternation gives W₀₀ ≈ 0.2 W₀: a third is the ceiling set by the definition and the half duty cycle.
- **At equal total probe power** a TE₀₁ probe has s/2 per component. Comparing the s = 0.5 row of TE₀₁ with alternation at s = 1 (the same power in one LP₁₁ lobe): P₀₀ = 0.96 vs 0.93 (B) and 0.95 vs 0.93 (A); W₀₀ = 0.023 vs 0.016 ν_r (B) and 0.020 vs 0.015 ν_r (A). **TE₀₁ is ≈ 1.4× faster and colder at equal power.**

**Longitudinal field, worst case** (`vector_carrier.dat`: E_z on the probe transition, converged 1D master equation):

| Case | s | W₀ without → with E_z (s⁻¹) | P₀ |
|---|---|---|---|
| A | 0.3 / 1 / 3 | 1079 → 1067 / 2225 → 2086 / 3694 → 2760 | unchanged or higher (0.954 → 0.956 at s = 1) |
| B | 0.3 / 1 / 3 | 6561 → 6077 / 13 838 → 9560 / 20 705 → 6673 | unchanged or higher (0.972 → 0.973) |

- The quadrature (iε_z) and real (ε_z) forms agree within a few percent.
- Case A loses 6 % at s = 1 and is harmless.
- Case B loses 31 % at s = 1 and 68 % at s = 3: there the polarisation geometry matters (Section 1 rules), or TE₀₁/HE₂₁ should be used.
- P₀ is not degraded, consistent with the purity study (P₀ unchanged up to ε ≈ 3 η⊥).

## Open

- A mode solver for the real ARF/NANF and bandgap structures: splitting, E_z, and the purity of a launched TE₀₁.
- Steady P₀₀ at small s and the one-axis floor at weak probe in the two-channel scheme.
- The real Zeeman structure, now worked out at the level of couplings in the background report (Sec. "One vector mode for both radial axes", figures `simulation/results/vector/`, data `vector_rb_sublevels.dat`). Remaining: the 2D Zeeman master equation (cluster).
  - **⁸⁷Rb D1, B ∥ z, σ⁺ coupling on F = 1:** the only valid target is g_p = |2,1⟩ via F′ = 2.
    - Channel 1: → |2′,2⟩ ← |1,1⟩, probe d² 1/6, coupling d² 1/2, Raman shift −μ_B B.
    - Channel 2: → |2′,0⟩ ← |1,−1⟩, probe d² 1/4, coupling d² 1/12, Raman shift 0.
    - So the probe needs two frequencies, and the coupling strengths are unbalanced by 6×.
    - Spectators: |2,2⟩ (σ⁻ only, resonant with ω₁), |2,0⟩ (σ⁺ only, resonant with ω₂), |2,−1⟩ (off resonance), |2,−2⟩ (no Λ).
  - **B ∥ x with π coupling:** also valid, with better-balanced legs (coupling 1/4 and 1/3).
  - **Expectation:** several times below the minimal model (remixing ~3×, imbalance ≲ 2.5×). With κ = 2–5, case A is excluded for two axes at 300 kHz (2h = 0.8–2) and case B is allowed (2h = 0.19–0.47).
- **¹⁷¹Yb cannot do the TE₀₁ scheme:** its F = ½ ground state has only two sublevels, so there are no two separate final states. It remains a good species for one-axis cooling (`notes/species_map.md`).
- Launching TE₀₁: s-waveplate / q-plate, or a photonic lantern. TE₀₁ also offers the azimuthal symmetry needed for a rotating-frame interpretation.
