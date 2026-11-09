# Cooling both radial axes: theory (9 Nov 2026)

Analytic only; nothing here has been simulated. **All results so far are for one radial axis** with the other axis ignored, including its heating.

## 1. One probe cools one motional mode

Near the node a probe field is f(x,y) ≈ c_x x + c_y y, so with x ∝ a + a†, y ∝ b + b†:

  f ∝ c_x a + c_y b + c_x a† + c_y b†.

On the red sideband (rotating-wave, Ω_R√n ≪ ν_r) only the annihilation part is resonant, so **a single coherent probe couples to one motional mode,** A = (c_x a + c_y b)/√(|c_x|²+|c_y|²). The orthogonal mode is dark in the resolved regime.

| Probe | c_x, c_y | Mode cooled | Left uncooled |
|---|---|---|---|
| LP₁₁, lobes along x | 1, 0 | a | b |
| LP₁₁ rotated by φ | cos φ, sin φ | x′ | y′ |
| Scalar OAM, l = +1, x + iy | 1, i | a₋ = (a + ib)/√2 | a₊ = (a − ib)/√2 |

**A single scalar OAM beam therefore does not solve it.** Angular momentum shows why: absorbing l = +1 and removing energy can only remove a quantum of the mode with L = −1. Removing an L = +1 quantum would need l = −1.
The creation part of x + iy adds L = +1 quanta, but it is the blue sideband, suppressed by (Γ/4Δ)².

## 2. Recoil couples the axes: the main change

One cycle removes one quantum from the addressed mode, but recoil heats **both** axes (α_x + α_y = 0.8 for the present dipole pattern). So per axis, with equal cycle rates:

  v = −(R/2)(1 − 2h), D = h n R → v/D = −(1 − 2h)/(2h n), so **effective h₂ = 2h** when both axes are saturated (R is shared).

The other limit: x hot and saturated, y cold, where y tracks R(n_y) = R_x h/(1−h). This gives h₂ = h/(1−h). The two estimates bracket the tail exponent 1/h₂. They coincide at h = 1/2.

**Energy balance:** net removal per cycle is 1 − 2h, so a cold steady state needs **h < 1/2, not h < 1.** (Fast channel not included, as in the 1D argument.)

**Core floor.** With rates A∓ ∝ (n, n+1), equal axes:

  n_core = ξ (1 + 2h)/(1 − 2h), ξ = (Γ/4Δ)²,  against ξ(1 + h)/(1 − h) in 1D.

| a (µm) | b_p | h | 2h | h/(1−h) | 2D/1D floor | P₀ both axes, 3Γ | 6Γ |
|---|---|---|---|---|---|---|---|
| 3.5 | 1/2 | 0.09 | 0.19 | 0.10 | 1.2 | 0.98 | 0.995 |
| 3.5 | 5/6 | 0.06 | 0.11 | 0.06 | 1.1 | 0.98 | 0.996 |
| 7 | 1/2 | 0.19 | 0.38 | 0.23 | 1.5 | 0.97 | 0.99 |
| 7 | 5/6 | 0.11 | 0.23 | 0.13 | 1.3 | 0.98 | 0.99 |
| 10 | 1/2 | 0.27 | 0.54 | 0.37 | 1.9 | 0.96 | 0.99 |
| 10 | 5/6 | 0.16 | 0.32 | 0.19 | 1.4 | 0.97 | 0.99 |
| 15 | 1/2 | 0.40 | 0.81 | 0.67 | 3.9 | 0.88 | 0.97 |
| 15 | 5/6 | 0.24 | 0.48 | 0.32 | 1.7 | 0.96 | 0.99 |
| 21 | 1/2 | 0.56 | 1.1 | 1.3 | no cold state | | |
| 21 | 5/6 | 0.34 | 0.68 | 0.51 | 3.2 | 0.93 | 0.98 |

Weak probe, core only, so upper bounds on P₀: no s dependence, no tail. **Case A with b_p = 1/2 (the current paper) moves into the window where the mean is set by the cut-off,** and the floor rises almost 4×.
With 2h: h₂ < 1/3 needs a < 10 µm for b_p = 5/6 and 6 µm for b_p = 1/2.

**Consequences.** The argument for small a and large b_p gets much stronger. It agrees with the expectation of a ≈ 7–10 µm from optical depth.

**Polarisation matters.** α depends on the emission pattern: dipole along z gives α_x = α_y = 0.4, radial sum 0.8; along x gives 0.2 and 0.4, sum 0.6; averaged, 2/3. Guided light is transverse, so the sum may be below 0.8. Must come from the Zeeman model.

## 3. Ways to cool both axes

| Scheme | Field | Works when | Rate per axis | Comment |
|---|---|---|---|---|
| **A. Alternating lobes** | LP₁₁a, then LP₁₁b | any switching | ≈ W/2 | simplest; needs mode switching ≲ 1/W ~ 0.1–1 ms |
| **B. Rotating pattern** | LP₁₁ at angle ω_rot t | W ≪ ω_rot ≪ Γ_n | ≈ W/2 | two frequencies δf = 2ω_rot within the narrow line; continuous |
| **C. Trap anisotropy** | one LP₁₁ at 45° | W ≪ abs(ν_x − ν_y) ≲ Γ_n | ≈ W/2 | beat between axes rotates the dark mode into the bright one; cw, single beam |
| **D. Two OAM** | l = +1 and l = −1, incoherent | no interference between them | ≈ W/2 | same as B if the frequencies differ |
| **E. Vector vortex** | TE₀₁ or TM₀₁ | the two circular polarisations drive independent Raman paths | possibly ≈ W | **most interesting** |

**C, numbers.** W ≈ 0.05–0.08 ν_r, Γ_n = 0.33 ν_r at 3Γ: needs anisotropy of ≈ 10–30 %. At 6Γ, W ≈ 0.03, Γ_n = 0.17: ≈ 5–15 %. A fixed lobe aligned with a trap axis cools only that axis, so the angle matters, including any rotation along z.

**E, why.** TE₀₁ has E ∝ J₁(ur/a)(−sin φ, cos φ). In the circular basis, σ⁺ carries e^{+iφ} ∝ x + iy and σ⁻ carries e^{−iφ} ∝ x − iy. Total angular momentum J = 0. So **σ⁺ cools a₋ and σ⁻ cools a₊: one beam addresses both modes,** provided the two polarisations drive distinct Raman paths with no interference. It is a true fibre eigenmode, with an exact axial zero, unlike LP₁₁, which is a superposition. **This needs the Zeeman model,** which is needed anyway. A magnetic field or atomic structure that mixes the paths would reduce it to one mode.

**Hot stage [hypothesis].** The one-mode rule is for resolved sidebands. For Ω_R√n ≫ ν_r the hill δ_B ∝ 1 + tan²θ(r) of a donut depends on r, so the Sisyphus channel should act on both axes, even with a scalar OAM beam, while the cold stage does not.
A possible protocol: donut or OAM when hot, scheme A, C or E when cold.

## 4. Numerical plan (later)

1. 2D rate chain P(n_x, n_y) with R(n), D from the 1D data and the recoil kernel on both axes; N ≈ 200 per axis, sparse. Tests h₂, thresholds, floors, t_prep.
2. Small exact check: eliminated ME on (n_x, n_y) with N ≈ 8–10 per axis; dimension (2N²)², sparse: compare with the chain, schemes A and C, at stiff h where N is small.
3. Rerun the map with h₂ in mind; which (a, b_p) still give P₀ ≥ 0.95 on both axes.
4. Scheme E with the Zeeman model.
