# Do higher-order modes help a fibre spin-wave memory?

Theory only; no master-equation simulation. All numbers are **estimates** for ⁸⁷Rb, $a = 15$ µm, $x_{0} = 0.09$ µm,
1 atom per lattice site over 1 cm unless stated. Companion to `paper_draft_ideas.md`, Section 6.

**Caveat (see `quantum-memory/memory_theory.md`):** the OD values below are for a closed two-level transition (r = 1) at full filling. A usable Λ has r = ¼ (clock) to ½ (stretched), and f ≈ ½, so divide them by 4–8. The ratios g are unaffected.

## Short answer

A HOM can play three roles, and they behave very differently.

| Metric | HOM as signal or control | HOM as the trap (dark trap) |
|---|---|---|
| Storage time | no direct gain | **potential gain,** ≈ 5–7× beyond a lowered red lattice |
| Retrieval efficiency | gain $\propto$ 1/OD, only if OD is small and the mode conversion is excellent | indirect |
| Bandwidth | gain $\propto$ OD; control power also falls | none |

The largest possible advantage is the **dark trap**, not storing the light in a HOM. It is also the riskiest, and a simply lowered red lattice gets much of the way without any HOM.
Separately, signal and control in *different* modes may give a practical gain in noise.

## 1. Why LP0m, and not LP11, for the signal

Atoms cooled to the axis behave as point scatterers, so what counts is the mode's on-axis intensity per unit power.

**LP₁₁ on the node: no.** Ground-state atoms see ≈ 4×10⁻⁴ of the LP₀₁ peak intensity, so OD ≈ 0.01. The node that makes
cooling carrier-free kills the memory coupling. Each absorption would also change n by one.

**LP₁₁ off-axis: not worth it.** The lobe peak is only 1.13× the LP₀₁ peak, and it would need off-axis traps and lose the node cooling.

**LP₀m: yes.** Normalising the on-axis intensity gives $g = J_1^2(u_{01})/J_1^2(u_{0m})$:

| Mode | g (cold atoms) | Same OD as LP₀₁ with a = | Central lobe FWHM (µm) | g relative to LP₀₁, both at 10 µK | OD (1 cm, two-level) | k difference, period | Loss scale (u/$u_{01}$)² |
|---|---|---|---|---|---|---|---|
| LP₀₁ | 1 | 15 µm | 14.1 | 1 | ≈ 30 | — | 1 |
| LP₀₂ | 2.33 | 9.8 µm | 6.1 | 2.21 | ≈ 70 | 0.91 mm | 5.3 |
| LP₀₃ | 3.66 | 7.8 µm | 3.9 | 3.18 | ≈ 110 | 0.32 mm | 12.9 |
| LP₀₄ | 4.99 | 6.7 µm | 2.9 | 3.85 | ≈ 150 | 0.17 mm | 24 |
| LP₀₅ | 6.32 | 6.0 µm | 2.3 | 4.23 | ≈ 190 | 0.10 mm | 38.5 |

Two conclusions:

- **Radial cooling and HOMs need each other.** Cold atoms get the full g; at 10 µK LP₀₅ loses a third of it, and the hot
  atoms also see an inhomogeneous coupling. The tail matters little here: $\langle I\rangle/I_{0}$ stays 0.99 at $\langle n\rangle = 0.15$, even for LP₀₅.
- **A HOM mimics a smaller core.** LP₀₂ at 15 µm has the OD of LP₀₁ at 9.8 µm. The smaller core also raises $\nu_r$ and lowers h.
  So a HOM signal is useful only when something else fixes the core large: loading, wall distance, loss or cooling power.

## 2. Retrieval efficiency

Optimal storage has an error that scales as $\kappa/d$ (Gorshkov 2007). **The coefficient $\kappa \approx 3$ below is from memory and must be
checked,** along with forward or backward retrieval and the real $\Lambda$ line strength, which lowers d.

The gain must pay for mode conversion, in and out. With fibre-to-HOM conversion efficiency T per pass, a HOM helps only if

$$\frac{\kappa}{d}\left(1 - \frac{1}{g}\right) > 1 - T^2$$

| LP₀₁ OD d | Mode | Error without HOM | Needed T per pass |
|---|---|---|---|
| 30 | LP₀₂ | ≈ 10 % → 4.3 % | ≥ 0.971 |
| 30 | LP₀₃ | ≈ 10 % → 2.7 % | ≥ 0.963 |
| 70 | LP₀₂ | ≈ 4.3 % → 1.8 % | ≥ 0.988 |
| 140 | LP₀₂ | ≈ 2.1 % → 0.9 % | ≥ 0.994 |

**So at OD of about 30 the gain needs better than ≈ 96–97 % conversion per pass; above about 100 it is not worth it.** In a
small core or a long column there is no efficiency case. The P₀⁴ ceiling is independent of this.

Other costs: a LP₀₁ admixture in a control beam makes on-axis stripes along z (period 0.9 mm for LP₀₂), of $\pm 2\sqrt{p/g}$ in
intensity for power fraction p (≈ ±4 % for LP₀₂, $p = 10^{-3}$). In a signal it is a loss. Loss per metre scales at least as u², though 1 cm tolerates a lot;
many low-loss fibres strip HOMs by design, so LP₀m may simply not be guided in the best fibres.

## 3. Bandwidth

- Adiabatic optimal storage needs $Td\gamma \gg 1$ (Gorshkov 2007), so at fixed error the **pulse can be shorter by g**.
  Scale: $d = 30$, $\gamma = \Gamma$/2 gives $d\gamma \approx 5\times10^8$ s$^{-1}$; at $T \approx 10/(d\gamma)$, about 50 MHz. Not a precise limit.
- **The same power gives g times the control Rabi frequency squared.** Power is not the bottleneck at these OD values, but it lowers heating of the fibre and scatter.
- **Temporal multimode capacity** $\propto \sqrt d$ for EIT and Raman (Nunn 2008), so only $\sqrt g$.
- Modal dispersion over 1 cm is negligible.
- No cavity linewidth applies. That advantage belongs to the fibre, not to the HOM.

**New idea, speculative: wavevector multiplexing.** On-axis atoms cannot store several transverse modes independently, as they are
point-like. But a signal in LP₀m with control in LP₀m′ writes a spin wave of wavevector $\Delta\beta = \beta_m - \beta_{m'}$. Over L = 1 cm,
$\Delta\beta L \approx 69$ rad for LP₀₂ with LP₀₁ (6.94 mm⁻¹), so the spin waves are nearly orthogonal and could form
parallel channels, addressed by mode pair. Pairs (m, m′) give more channels than modes. Needs a crosstalk calculation, and the
channels share the atoms. In a lattice the axial Lamb–Dicke parameter stays $\ll$ 1, so these spin waves do not dephase; this would
not work with free atoms.

## 4. Storage time: HOM as the trap

### The problem with a red lattice, which also concerns the present paper

At an antinode the absolute differential shift is the full $\delta_{\rm diff} \approx 800$ Hz at 530 µK. Then:

1. **Any variation of lattice intensity along the column dephases the spin wave.** 1 % → **≈ 8 Hz,** far larger than the
   ≈ 0.01 Hz radial value quoted in the paper, and a limit radial cooling cannot remove. It needs compensation (magic-field or
   light-shift compensation, Dudin 2010), a spin echo, or a low depth.
2. **Scattering ≈ 4 s⁻¹** heats ≈ 0.47 radial quanta s⁻¹ per axis, destroying the state that was cooled.

### The dark trap

Blue-detuned TE₀₁-type donut for radial confinement, plus a blue LP₀₁ standing wave for axial confinement: the atom sits at a
3D intensity zero. Blue traps are known to reduce scattering and extend coherence (Ozeri 1999; bottle beam 2010); we found no
in-fibre version in a short search.

In a harmonic dark trap, the seen intensity is $\langle V\rangle/U_0 = \hbar(2\nu_r + \nu_z)/4U_0$, so ≈ 7×10⁻³ at $\nu_z = 300$ kHz. Equivalent to one
lattice depth $U_0/k_B$ = 530 µK:

| Trap | Scattering (s⁻¹) | Radial heating (quanta s⁻¹ per axis) | Absolute shift (Hz) | Effect of 1 % intensity (Hz) |
|---|---|---|---|---|
| Red 1064 nm, antinode | 4.2 | 0.47 | 800 | 8 |
| Red, lowered ×100 | 0.04 | 0.05 | 8 | 0.08 |
| Blue 760 nm, dark | 0.25 | 0.055 | 50 | 0.25 |
| Blue 532 nm, dark | 0.016 | 0.007 | 3 | 0.016 |

**Read this carefully.** Only lowering the red lattice by ×100 already matches the 760 nm dark trap, and it needs no HOM. Against that
lowered lattice, 532 nm gains only **≈ 7× in heating and ≈ 5× in shift sensitivity** (65× and 500× against the full-depth lattice). The dark trap gains
more when combined with its own lowering: in a dark trap the heating in quanta does not depend on depth, and shift sensitivity falls as $\nu$,
and lowering $\nu_z$ reduces scattering further. Against a lowered lattice, the case for a dark trap is real but not overwhelming.

A 3D ground-state atom would then be limited by other things: magnetic noise, trap noise and, possibly first of all, background gas in the fibre.

**A bonus for cooling.** For the same depth, the TE₀₁ radial frequency is 1.94× higher than for LP₀₁, so **h falls by 1.94** (0.40 → 0.21
at $a = 15$ µm, $b_p = 1/2$), which puts it below the h = ⅓ threshold. The trap light shift vanishes at the centre, so the cooling resonance is unshifted.

**Risks.**
- All TE₀₁, TM₀₁, HE₂₁ have an axial zero, so the node is robust, but a superposition that becomes a linear LP₁₁-like pattern
  loses confinement along one axis. The polarisation beat along the fibre must be long compared with the column.
- A LP₀₁ admixture p adds intensity ≈ p at the axis and displaces the trap centre, periodically along z. Needs p ≲ 10⁻³, comparable
  to the LP₁₁ probe purity already required.
- The fibre must guide 760 or 532 nm with a clean donut, plus 795 nm. Likely a hard fibre problem at 532 nm.
- Loading into a dark guide; transfer after cooling is an option. $\epsilon$ is 9× larger at 760 nm.
- The 760 nm numbers use a two-line model; near D2 the estimate is rough. Check 532 nm for other Rb transitions and for photoionisation.

## 5. Noise: possibly the most practical gain

Single-photon operation in a collinear fibre geometry must separate a few photons from a strong control beam, and suffers
four-wave-mixing noise. With signal and control in orthogonal modes:

- **A mode filter** adds to the frequency and polarisation filtering. Its extinction would need estimating, and scattering can undo it.
- **Phase mismatch.** FWM is known to be suppressible by mismatch (see the discussion in Thomas 2019); here $\Delta\beta L \approx 70$ rad.
  **Whether FWM is really mismatched depends on the mode of the noise field; to be calculated.** We found no paper on this control by mode.

## Verdict

| Idea | Advantage | Confidence | Next step |
|---|---|---|---|
| Dark HOM trap, 532 nm | storage time, ≈ 5–7× beyond a lowered lattice, plus lower h | moderate; feasibility unknown | fibre feasibility; full $\eta$(T) |
| Lowered red lattice (no HOM) | storage time | high | **use as the baseline** |
| LP₀m signal | efficiency, bandwidth | good, but only for low OD and $T \gtrsim 0.97$ | fix $\kappa$; realistic T |
| Signal and control in different modes | noise | unknown | FWM calculation |
| Wavevector multiplexing | capacity | speculative | crosstalk calculation |
| LP₁₁ signal on the node | none | high | drop |

**A fair claim must compare with the best no-HOM alternative,** in particular a smaller core and a lowered lattice.

## What could go in the paper

- A short subsection: LP₀m gain table, and the dark-trap estimate as an outlook.
- The longitudinal inhomogeneity (8 Hz per 1 %) must go in the memory section whatever else is decided.

## References

Verified on Crossref (DOI, title): Gorshkov 2007 PRL 98, 123601; Nunn 2008 PRL 101, 260502; Ozeri 1999 PRA 59, R1750 (blue-detuned trap);
Trapping a single atom in a blue-detuned bottle beam, Opt. Lett. 35, 2164 (2010); Thomas 2019 PRA 100, 033801; Dudin 2010 PRL 105, 260502;
Gouraud 2015 PRL 114, 180503 (nanofibre memory, for context).
Full texts not read for this note. $\kappa \approx 3$ is from memory.
