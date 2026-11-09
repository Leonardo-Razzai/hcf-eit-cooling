# Paper draft: ideas and structure

Working document for the author and collaborators. **Nothing here is final.**
Status tags: **[done]** computed and converged; **[model]** analytic estimate, not simulated;
**[to do]** not computed yet; **[decide]** needs a decision.

## Focus (author's decision, 9 Nov 2026)

**Radial ground-state cooling of atoms inside a hollow-core fibre with an LP₁₁ probe.** Everything else serves that claim. The heavy-tailed
statistics explain how the cooling behaves in this soft trap; the memory is why it matters.

**Claim we could make, if the open items hold:** all-guided radial ground-state cooling of a cm-long atomic column in a hollow-core fibre, with nanowatts of
LP₁₁ light and no transverse optical access.

**Prior work to state clearly.** Targeted search only, no paper found on LP₁₁ or node-based ground-state cooling in a hollow-core fibre. But radial cooling in a fibre is
not new in general: in-fibre gray molasses and delta-kick cooling (arXiv:2112.10088), a hyperfine-selective radial scheme and its Cs study (Waterloo, DAMOP 2024; full
paper not found). Guided-light ground-state cooling exists for nanofibres (Meng et al. arXiv:1712.05749; Østfeldt et al. arXiv:1708.03488), and the node idea for ions in
integrated waveguides (Xing & Mehta 2025; arXiv:2512.10900). So the claim is specifically: **ground state, radial, inside a hollow-core fibre, through a higher-order guided mode.**

**Relation to Xing & Mehta, PR Applied 24, 014034 (2025), and arXiv:2512.10900 (compared 9 Nov 2026).** Already theirs, for ions: carrier-free EIT at a node, a rate independent of the gradient, wider bandwidth, tolerance to modal impurity and position, master-equation simulation with Zeeman levels, and an experiment ($\bar n = 0.05$). So the concept is not ours, and most of the present manuscript's Figs 3 and 4 is a translation of it to another system. They work at $\eta \approx 0.1$, MHz traps, where recoil is a small correction. **Ours:** the recoil-dominated regime ($h \sim 0.1$–$1$, heavy tail, failure of $W$, one law for all regimes), two axes with recoil coupling, alkali dark states and remixing, the fibre (modes, purity along a column, $\zeta$), capture from about 10 µK, memory.

**Weak points of this claim, to fix first:**
1. **Only one radial axis is simulated, and an LP₁₁ lobe pattern addresses one axis.** Without the second axis there is no radial ground state.
   Theory (9 Nov 2026, `notes/anomalous-cooling/two_axes_theory.md` on `universal`): one probe cools one motional mode, so a single scalar OAM beam is not enough. Recoil
   heats both axes, so h → about 2h, and a cold state needs h < 1/2. At ν_z = 300 kHz, case A with b_p = 1/2 becomes marginal. **The limits scale as a/ν_z:** a deeper lattice relaxes them (costs: power ∝ ν_z²a², purity ∝ 1/(aν_z), excited-state light shift). Cool deep, store shallow.
   TE₀₁/TM₀₁ may cool both modes with one beam (to be checked with the Zeeman model).
2. Hyperfine-averaged Λ: Zeeman dark states of the best configuration.
3. Real-fibre modes: TE₀₁/TM₀₁/HE₂₁ structure, polarisation, beating along z.

## 1. Introduction

- Cm-long atom columns in a guided mode: large optical depth without a cavity. Axial motion can be
  ground-state cooled, radial motion cannot be reached by co-propagating beams.
- The idea: coupling in LP₀₁, probe in LP₁₁, atoms on the node. Carrier-free, the fibre analogue of an
  ion at a standing-wave node (Cirac 1992).
- Do not repeat Morigi or Cirac. State at once what is new: the soft-trap regime.

## 2. Models, from simplest to numerics (universal)

Everything in reduced form: rates in units of ν_r; only ν_r/ω_rec, b_p, Δ/Γ and s matter (Γ/ν_r and x₀/a
drop out, < 0.3 %).

| Level | Model | Gets right | Fails where | Status |
|---|---|---|---|---|
| 0 | Sideband rate equation (Morigi) | weak-probe core, with κ = 1/b_p | s ≳ 0.1; n_ss up to 5× low, W up to 40 % high | [done] |
| 1 | Add saturation of the red sideband, recoil per cycle h | rate scale, threshold h = 1 | no tail, no Raman splitting | [done] |
| 2 | Generalised rate model (c_sat = 2.25) | drift at n ≲ 15, b_p = ½ | b_p ≠ ½, near threshold, drift above n ≈ 15 | [done] |
| 3 | Fokker–Planck from v(n), D(n) | tail for n ≳ 5, exponent ≈ 1/h | n ≲ 5, where single jumps ≈ n | [done] with measured v, D |
| 4 | Markov chain with exact recoil kernel and measured R(n) | should reach N of several thousand and thresholds | not built | [to do] |
| 5 | Eliminated, banded master equation | everything in the closed Λ, N up to ≈ 550 | band in regime III checked only up to n̄₀ = 10 | [done] |
| 6 | Full master equation | reference | N ≲ 60 | [done] |

**Main figure:** one quantity (P₀, ⟨n⟩ or the drift) against s or n for all levels, with a shaded region
where each stops being valid. Avoid re-deriving the standard theory.

**A new point for this section: which rate?** W from ⟨n⟩ depends on N and has a pole near n_ss ≈ 0.2. W₀
from P₀ does not. At case A, s = 1: W = 1319 → 1614 s⁻¹ and W₀ = 2188 → 2225 s⁻¹ from N = 20 to 120.
The two measure different things; W₀/W ≈ 1.2 at weak probe has no agreed explanation yet. W₀ is meaningful
only for cold states (P₀ ≳ 0.85). [done, on `universal`] **[to do]:** does 1/W₀ give the time from thermal
states? (Study running.)

### The cooling map

Plane (ν_r/ω_rec, b_p), with h = 1/3, 1/2 and 1 drawn on it and the real configurations marked. Grid chosen:
ν_r/ω_rec = 1, 1.4, 1.97 (case A), 3, 5, 8.6 (case B), 15; b_p = 1/6, 1/3, 1/2, 2/3, 5/6; Δ = 3, 6 Γ;
s = 0.3, 1, 3. **[to do; study written, queued on `universal`.]**

Panels: P₀; W₀/ν_r; ⟨n⟩/(1−P₀) (tail diagnostic); N-sensitivity; t_prep(P₀*) if the thermal study supports it.

Atoms other than Rb enter as points on this map only (Sr later). Only b_p = 1/6, 1/2 and 5/6 exist for Rb;
1/3 and 2/3 are there to show the map is continuous.

## 3. Anomalous cooling (main physics result)

- **Three regimes of removal rate against s·n:** R ∝ n (s n < 1); saturated plateau; a fast, Sisyphus-like
  channel above s n = 2(ν_r/Γ_n)². The curves collapse against s n. [done, numerics] The Sisyphus reading is a
  hypothesis; supporting it: from n = 20 to 200 at s = 1 the drift grows 17×, the diffusion per quantum only 1.8×.
  **[to do]** the sign-of-Δ test, the Δ scan of the onset, the semiclassical friction.
- **Power-law tail,** P(n) ∝ n^(−1/h): measured −2.34 against −2.48 at s = 0.1. [done]
- **Three thresholds, four windows:** h = 1/3 (variance), 1/2 (mean), 1 (ground-state cooling). They are
  crossovers, of width δh ≈ h²/ln(2(Δ/Γ)²). [model] Variance grows as N^0.54 against the predicted N^0.5 in case A. [done]
- **Above h = 1,** the fast channel may give a hot steady state, ⟨n⟩ ≈ (h−1)n_c, not a runaway. [model; scaling
  not confirmed: 3.8 vs 2.69 at s = 1, 0.76 vs 1.44 at s = 5]
- **Why it is fibre-relevant:** the soft trap puts the crossovers inside the thermal population; the hill is
  the smooth mode profile, not λ/2; the node keeps low n carrier-free.
- **Literature check done (9 Nov 2026): the tail law itself is not new.** Sisyphus cooling in a harmonic trap has P(E) ∝ E^(−1/D), the mean
  diverging at D > 1/2 and no stationary state at D > 1 (Dechant, Kessler, Barkai PRL 115, 173006; Dechant et al. PRE 94, 022151): the same
  structure with D ↔ h, from constant energy removal and diffusion ∝ E. **Cite these prominently.**
- **What we can claim as new, after a targeted search:** the heavy tail in *resolved-sideband ground-state cooling*, with h set by the trap and b_p; coexistence
  with P₀ > 0.95 and a tail starting near n ≈ 3 (Sisyphus: above ~10³ E_r); the cut-off, threshold widths and hot state above h = 1; the failure of W; the memory
  consequences; the LP₁₁ fibre geometry.
- **Also cite:** carrier-free EIT at a node for ions, Xing & Mehta PR Applied 24, 014034 (2025), and experiment arXiv:2512.10900; saturated rate in the dark-polaron
  theory arXiv:2607.22243; EIT for atoms beyond Lamb–Dicke, Roghani & Helm PRA 77, 043418. Our "fibre analogue of an ion at a node" is now also an analogue of these.
- **Name:** prefer "heavy-tailed Fock statistics" to "anomalous cooling".
- **Open:** read Roghani & Helm in full, and the Waterloo transverse-cooling result in a fibre (DAMOP 2024).

**Figure:** drift collapse, tail with the fit, window map with data points.

## 4. ⁸⁷Rb: best core radius and D1 vs D2

### What "best core radius" should mean [decide]

Scalings for ⁸⁷Rb D1, ν_z = 300 kHz. [model]

| a (µm) | ν_r (kHz) | h (b_p = ½) | h (b_p = 5/6) | OD, 1 atom/site, 1 cm, two-level | ε (Hz per quantum) |
|---|---|---|---|---|---|
| 3.5 | 30.6 | 0.09 | 0.06 | ≈ 550 | 1.1 |
| 7 | 15.3 | 0.19 | 0.11 | ≈ 140 | 0.55 |
| 10 | 10.7 | 0.27 | 0.16 | ≈ 70 | 0.39 |
| 15 | 7.15 | 0.40 | 0.24 | ≈ 30 | 0.26 |
| 21 | 5.1 | 0.56 | 0.34 | ≈ 15 | 0.18 |
| 31 | 3.5 | 0.83 | 0.50 | ≈ 7 | 0.13 |

OD scales from the background report (16 at a = 15 µm for N = 10⁴) with ≈ 1.9×10⁴ sites per cm. The real Λ
line strength will lower it. ε is the differential light shift per radial quantum per axis.

Candidate definitions:

1. **Largest a with a clean tail:** h < 1/3, so a ≈ 21 µm for b_p = 5/6. Simple, but ignores the application.
2. **Fastest t_prep(P₀*):** favours small a at fixed depth, since rates ∝ ν_r and h falls. No interior optimum
   from cooling alone; a is limited by loading.
3. **Maximum memory efficiency at a target storage time** under fibre constraints. **Recommended:**

   η_mem(a; T) = η_RW(OD(a)) × |φ(T)|⁴ × e^(−2T/τ_other)

   with constraints: loading and atom number, wall distance, LP₁₁ existence and loss at 795 nm, cooling power.

Expectation, not a result: OD pushes towards a ≈ 7–10 µm, the range between our two cases; h does not fix a by itself, because ν_z can compensate (h ∝ a/(b_p ν_z)). Case A was
never optimised. Needs fibre input (loss, purity, loading) that the simulations cannot give.

### D1 vs D2 (case A, s = 1, all excited levels, δ_e selected on W₀) [done, on `universal`]

| Configuration | W₀ at 3Γ (s⁻¹, both signs) | P₀ at 3Γ | Note |
|---|---|---|---|
| D1 F′=2 (present paper) | 2.2×10³ | 0.953–0.955 | robust; within 1 % of three levels |
| D1 F′=1, probe on F=2 | 3.2–3.6×10³ | 0.974–0.979 | b_p = 5/6, h = 0.24 |
| D2 F′=2 | 2.1–2.2×10³ | 0.941–0.955 | b_eff falls at large Δ; tail worsens (h ≈ 0.57 at 12Γ) |
| D2 F′=1, probe on F=1 | 2.0–4.2×10³ | 0.947–0.974 | factor 2 by interference sign, stops near 8Γ |

For P₀ = 0.99 (sign +1): D1 F′=1 at 5Γ, W₀ ≈ 1.9×10³ s⁻¹; D1 F′=2 near 6Γ, ≈ 1.0×10³ s⁻¹.
**D1 F′=1 is the candidate, subject to Section 5.**

Rerun all of this at the chosen a. Rb imperfections (mode purity, two-photon error, ν_r mismatch) are still at
N = 20 and must be redone. **[to do]**

## 5. Zeeman structure and repumping [to do; the main open risk]

All results above use a hyperfine-averaged Λ.

**Why it may matter.** F=2→F′=1 is an F→F−1 transition. For a fixed polarisation it has Zeeman dark states
(for example |2,±2⟩ for π light), independent of motion, which accumulate population and stop the cooling.
b_p = 5/6 then no longer describes the cycle. D1 F′=1 is at least closed in the hyperfine sense: both grounds
are in the Λ. [standard argument; to be worked out for our polarisations]

**Questions:**
- Which polarisations do LP₀₁ and LP₁₁ give at the node? LP₁₁ is a near-degenerate TE₀₁/TM₀₁/HE₂₁ family; its
  polarisation structure in a real fibre must be included.
- Does the gradient-only coupling still have a useful dark state among the Zeeman sublevels?
- Interference signs through F′=1, 2, now computed instead of bracketed.

**Options for destabilising dark states (compare, then choose):**
1. Weak repumper on F=2→F′=2, off-resonant from the cooling Λ. Costs extra scattering, raising κ and h.
2. Transverse magnetic field: Larmor mixing faster than the pumping. Cheap, but it shifts the two-photon resonance.
3. Mixed or switched polarisation.
4. Microwave or RF mixing between sublevels.

Figure of merit: P₀ and t_prep including the extra photons. A repumper raises κ, so it moves the point on the
map. That is a useful link between the sections.

**Numerics.** 8 ground sublevels. The banded system has S ≈ 64 N(2B+1) unknowns, dense:
N = 30, B = 3 → 1.3×10⁴ (≈ 3 GB, feasible); N = 60, B = 4 → 3.5×10⁴ (≈ 19 GB, too large here).
Since b_p = 5/6 has a light tail, try N ≈ 30 first, with a convergence check. Beyond that: Krylov solver or
quantum trajectories. `zeeman_coupling` already exists. Also include both radial axes at least by alternating stages.

Reference for ion tripod schemes: Lechner et al. 2016. Dark-state destabilisation: Berkeland & Boshier 2002 (from memory; verify).

## 6. Application: spin-wave memory in a cm-long column

**The honest position: we have not shown that this outperforms existing memories.** What follows is the
argument and what must be computed.

### State of the art

| Memory | Efficiency | Lifetime | Notes |
|---|---|---|---|
| Yang 2016, 3D lattice, ring cavity | 76(5) % | 0.22(1) s (1/e) | best joint result found; cavity finesse 52 |
| Bao 2012, ring cavity | 73(2) % | 3.2(1) ms | |
| Hsiao 2018, EIT in dense cold gas | 92.0(1.5) % | time-bandwidth 1200 | coherent pulses |
| Wang 2019, cold Rb, EIT | > 85 %, fidelity > 99 % | not extracted | single-photon polarisation qubits |
| PRL 2025, spin-wave compaction | 94.6 % reported | not extracted | taken from an APS synopsis; original not read |
| Dudin 2013, 1D lattice | low (to confirm) | "time scale of a minute" (light storage) | abstract not retrieved |
| Peters 2019/2020, cold atoms in HCF | to extract | short (free fall) | SNR 59 at single-photon level |

Lifetime × efficiency is the gap: the best efficiencies and the best lifetimes are in different experiments.

### What the fibre can offer

1. **No cavity:** high OD from confinement, so bandwidth is not limited by a cavity linewidth (9.2 MHz in Yang 2016),
   which also leaves room for temporal multimode storage. Fibre-coupled output, no free-space mode matching.
2. **Motional state:** in the ground state the radial light-shift dephasing disappears as a *decay*. [model]
   Coherence |φ(T)|² → ≈ P₀² for T ≳ 1/(2ε) and stays there, so the efficiency is capped at ≈ P₀⁴, not reduced further.
   P₀ = 0.99 gives 0.96, P₀ = 0.954 gives 0.83. **So the P₀ gain from D1 F′=1 directly matters.**
3. **Long-wavelength spin wave:** co-propagating Λ beams, no motional dephasing from axial transfer.

### What limits it [model unless stated; all to check]

- **Lattice photon scattering, likely the real limit.** At U₀/k_B = 530 µK, about 4 s⁻¹. Rayleigh events preserve the
  hyperfine coherence (Ozeri 2005) but kick the atom: ≈ 0.45 radial quanta s⁻¹ per axis in case A. So P₀ decays
  within ≈ 1 s, and the coherence with it. **Our ground state would not last unless the trap is lowered.**
- **Adiabatic lowering after cooling.** n is conserved, ε ∝ √(U₀), scattering ∝ U₀, heating in quanta ∝ √(U₀).
  A factor 100 in depth leaves ≈ 55 lattice recoil energies (still deep, but check tunnelling and gravity): heating ≈ 0.045 quanta s⁻¹ per axis, ~20 s scale. **Idea to simulate.**
- Second radial axis (not simulated) and axial cooling (assumed, from Raman sideband).
- **Lattice intensity variation along the column, possibly the first limit.** At an antinode 1 % gives ≈ 8 Hz (δ_diff ≈ 800 Hz),
  far above the ≈ 0.01 Hz radial value in the present paper. Needs compensation, echo or a low depth.
- Magnetic field noise: clock states or magic field. Lattice intensity and pointing noise. Background gas inside
  the fibre (a poor vacuum conductance): **may be the actual limit; check the fibre experiments.**
- Collisions if > 1 atom per site, hence the OD limit used above.
- Inhomogeneity along the column, in particular from fibre loss.

### What to compute

- η_RW(OD) for the Λ used (Gorshkov 2007): the 1/OD scaling and forward vs backward retrieval.
- Full η_mem(T) at several a and depths, with heating. **Claim only what this shows.**
- Target that would be a real advance: **η ≳ 80 % with τ ≳ 1 s, no cavity,** which would beat Yang 2016 in both.
  A fair comparison must state the conditions, as these benchmarks differ in protocol and in what they call efficiency.
- Interferometry and clocks: a shorter section at most. T_QPN ≈ 0.4–0.5 s in case A; relevant only at long T.

### Higher-order modes

Worked out in `hom_memory_theory.md`. In short: LP₀m raises the on-axis OD by g = 2.3–6.3 for cooled atoms (a real link between
radial cooling and the memory), which helps efficiency only at low OD and with conversion ≳ 0.97 per pass, and helps bandwidth ∝ g.
A LP₁₁ signal on the node is useless. The largest possible storage-time gain is a blue dark trap in a donut HOM, which also lowers h by 1.94.
Compare always with the no-HOM alternatives: a smaller core and a lowered red lattice.

## Proposed figures

1. Scheme and regimes (keep existing Figs 1, 2).
2. Model hierarchy against numerics.
3. **Cooling map.**
4. Anomalous: collapse, tail, windows.
5. Best a; D1 vs D2.
6. Zeeman: with and without repumper.
7. Memory: η_mem(T) against benchmarks.

Fig. 4 of the current paper (robustness) can go in compact form with Fig. 5 or in Methods.

## Work list, in order

1. **Both radial axes:** alternating LP₁₁ lobes, a rotating pattern or a donut. Rate and P₀ per axis, heating of the other axis during each stage.
2. **Zeeman model,** D1 F′=1. It can overturn the choice of configuration.
3. Real fibre modes, purity and beating along z; a realistic core radius; definition of best a.
4. t_prep(P₀*) on the map; the deferred s-ramp steps.
5. Memory model with lattice scattering, lowering and inhomogeneity.
6. Fast-channel confirmations (Δ < 0, Δ scan); Rb imperfections at converged N.
7. Read: Roghani & Helm 2008; find the Waterloo paper.

## Where to publish

**Recommendation: PRX Quantum,** if both axes and Zeeman hold and the memory section is quantitative.

- **Cooling only, no memory result:** Physical Review Research or Physical Review A; **Physical Review Letters** is possible for a short paper on the concept plus the one-formula result, but it would need both axes done.
- **Nature Communications:** realistic only with a quantified memory advantage or an experimental collaborator. A pure proposal on cooling in one geometry is a hard sell.
- **Nature Physics / Nature Photonics:** not without an experiment.

## Open questions for the author

- How much of the present manuscript's two-fibre comparison to keep, given the move to the universal map?
- Fibre input for "best a": which fibres are realistic, and what loading is achievable?
- Is an experimental collaborator available? It strongly affects journal options.

## Reference status

Verified on Crossref (DOI, title): Yang 2016 Nat. Photon. 10, 381; Bao 2012 Nat. Phys. 8, 517; Dudin 2013 PRA 87, 031801;
Hsiao 2018 PRL 120, 183602; Wang 2019 Nat. Photon. 13, 346; Gorshkov 2007 PRL 98, 123601; Ozeri 2005 PRL 95, 030403;
Lechner 2016 PRA 93, 053401; 2025 PRL "Near-perfect broadband quantum memory enabled by intelligent spin-wave
compaction" (10.1103/kbwj-md9n).

Numbers taken from arXiv abstracts: Yang (1511.00407), Hsiao (1605.08519), Wang (2004.03123), Peters SNR (1906.05771);
Bao from an arXiv listing. **Still to verify:** Dudin efficiency and lifetime, the 2025 PRL figures, Berkeland & Boshier 2002.
