# Spin-wave quantum memory with ground-state atoms in a hollow-core fibre

Theory note, written so that a separate project, in a new Claude Code session, can start from it **without the cooling
code**. The cooling project only supplies the initial motional state (Section 1.3). All numbers are estimates for ⁸⁷Rb on
the D1 line, from two light scripts in this folder:

| File | What it does | Runtime |
|---|---|---|
| `memory_budget.py` | line strengths, retrieval kernel, efficiency against OD, column OD, phase mismatch, temporal capacity, read-out erasure, control photons, noise, lattice depth, gravity, vacuum, fibre vibration, $\eta(T)$ scenarios → `budget_output.txt`, `memory_budget.png` | ~8 min, 2 threads |
| `crosstalk.py` | spin-wave crosstalk model (Section 7) → `crosstalk_output.txt`, `crosstalk.png` | seconds |

Tags: **[script]** computed by these scripts; **[estimate]** back-of-envelope in this note; **[lit]** from a paper
(read as stated in Section 11); **[open]** not done. Companion notes: `../hom_memory_theory.md` (higher-order modes),
`../paper_draft_ideas.md` Section 6.

---

## 0. Summary

1. **The bottleneck is optical depth, not storage time.** With one atom per lattice site, OD $= N\sigma/A_{\rm eff}$.
   For a = 15 µm, a 1 cm column and half filling, OD ≈ 15 on a closed transition, and only **3.7 to 7.4 on a usable Λ**.
   That gives 23–40 % efficiency. Efficient storage (η ≥ 0.9) needs OD ≈ 90–100, so it needs a **small core**
   (a ≈ 3.5–5 µm gives OD 35–140), higher filling, a longer column, or a cavity. Small cores also lower the cooling parameter h, so the
   memory and the cooling point the same way. Vacuum and loading point the other way. [script]
2. **Achievable efficiency** (optimal control, backward retrieval) is $1-\eta \approx 11.6/{\rm OD}$ at large OD: 0.90 at OD = 93,
   0.95 at 201, 0.99 at 1100. In a cm-long column, backward retrieval loses a few per cent to the hyperfine phase mismatch
   (2.9 rad/cm). [script; Gorshkov 2007]
3. **Storage time.** Cool deep, store shallow. With depth lowered 100× and a spin echo, the model gives $T_{1/e}$ ≈ 2 s when limited by
   background gas at an assumed $\tau_{\rm bg}$ = 2 s. With $\tau_{\rm bg}$ = 100 s it gives **≈ 11 s, then limited by radial heating from lattice scattering.** Without echo,
   longitudinal lattice-intensity variations are the first limit: 1 % rms gives 21 ms at full depth and 1.3 s at 1 % depth. [script]
4. **What the radial ground state buys** is not mainly a longer $T_2$: a spin echo also removes static thermal dephasing.
   It buys (i) shallow traps without loss or anharmonic broadening, (ii) homogeneous coupling, which matters for efficiency,
   for LP₀ₘ signals and in small cores, (iii) $(2\bar n+1)$ ≈ 60× fewer motional jumps that defeat the echo (case A from 10 µK), and (iv) operation without echo
   up to the $P_0^4$ ceiling. A fair claim must compare against thermal atoms with an echo. [estimate]
5. **Multimode capacity is modest in one dimension.** Temporal (EIT/Raman): ≈ $0.18\sqrt d$ modes at η > 0.7, so about 4 at OD = 1000.
   Gradient echo (∝ d), mode pairs (LP₀ₘ signal / LP₀ₘ′ control, read in groups) and frequency-addressed segments with shelving are the
   candidates. The fibre's strength is efficiency × lifetime × fibre coupling, not mode number. [script/estimate]
6. **Crosstalk model (new, Section 7).** Collective bosonic operators give a crosstalk matrix built from a storage map and a read map.
   Independent per-atom noise transfers only $(1-|\chi|^2)/N$ into another orthogonal mode, which is negligible. Crosstalk comes from (a) phases
   **correlated** along z (Jacobi–Anger, $J_1^2(\phi_0)$; for example an LP₀₂ admixture in the lattice), (b) the **phase-matching tail**
   $X/X_0 \approx 1.6\,d/(\Delta k L)^2$, so channels must be separated by $\Delta kL \gg d$, not 2π, (c) **read-out erasure**: resonant EIT
   read-out destroys co-located spin waves, while far-detuned Raman leaves $\approx e^{-12.7/d}$ per access, and (d) random filling, $(1-f)/N$. [script]

---

## 1. Platform and parameters

### 1.1 Geometry

Atoms sit one per site (fraction f of sites filled) in a 1064 nm LP₀₁ standing wave on the fibre axis. The site spacing is 532 nm.
A 1 cm column has 18 800 sites. Signal and control are guided; for cold atoms on the axis the relevant intensity is the on-axis
one, $I(0) = P/A_{\rm eff}$ with $A_{\rm eff} = \pi a^2 J_1^2(u_{01})$ ≈ 0.85 a².

$${\rm OD} = \frac{N\,\sigma}{A_{\rm eff}},\qquad \sigma = r\,\frac{3\lambda^2}{2\pi},\qquad d = {\rm OD}/2 \ \ \text{(Gorshkov convention: intensity } e^{-2d}).$$

Here r is the branching ratio of the specific transition (r = 1 for a closed two-level transition).

### 1.2 Which Λ? Polarisation constraint in a fibre

Guided light is transverse. With the quantisation axis along the fibre (an axial bias field) only **σ±** exist, and there is **no π light**. Line strengths on D1 [script,
checked against `eitsim.atoms`]:

| g → s | via | r signal | r control | control couples g? (FWM channel) | Notes |
|---|---|---|---|---|---|
| \|1,0⟩ → \|2,0⟩ (clock) | F′=1, m′=±1 | 0.083 | 0.083 | 0.33 | weak |
| \|1,0⟩ → \|2,0⟩ (clock) | **F′=2, m′=±1** | **0.25** | 0.25 | 0.33 | reference scheme |
| \|2,1⟩ → \|1,−1⟩ (magic pair, 3.23 G) | F′=1 or 2, m′=0 | 0.25 | 0.083 | 0.17 | σ⁻ signal, σ⁺ control |
| **\|2,2⟩ → \|1,0⟩ (stretched)** | **F′=1, m′=1** | **0.50** | 0.083 | **0** | no Stokes channel; coherence field-sensitive |

The stretched scheme has twice the OD of the clock scheme and **no four-wave-mixing channel**: σ⁺ from \|2,2⟩ has no D1 partner.
Its g–s coherence shifts by 1.4 MHz/G. Store, then map it to the clock pair (for example with two microwave π pulses
\|2,2⟩ → \|1,1⟩ → \|2,0⟩), and map back before read-out. [estimate; to verify in S1]

### 1.3 Interface to the cooling project

Inputs a memory project needs from the cooling work:
- ground-state population per radial axis, $P_0$ (one axis: 0.95–0.98 in case A, depending on the scheme; both axes, $P_{00}$: 0.76–0.98 from the 2D chain, indicative);
- the shape of P(n), heavy-tailed, so ⟨n⟩ is not a good summary;
- the axial state (assumed Raman-sideband cooled);
- the preparation time (ms) and atom loss;
- the internal state after cooling (F = 1 or 2; the Zeeman repumper leaves atoms spread over sublevels, so optical pumping to g is a separate step).

### 1.4 Column optical depth [script]

| a (µm) | Column (mm) | f | N | OD (r = 1) | OD clock (r = ¼) | η, clock | OD stretched (r = ½) |
|---|---|---|---|---|---|---|---|
| 15 | 10 | 1 | 18 800 | 30 | 7.4 | 0.40 | 15 |
| 15 | 10 | ½ | 9 400 | 15 | 3.7 | 0.23 | 7.4 |
| 15 | 50 | ½ | 47 000 | 74 | 19 | 0.65 | 37 |
| 7.5 | 10 | ½ | 9 400 | 60 | 15 | 0.59 | 30 |
| 5 | 10 | ½ | 9 400 | 134 | 34 | 0.77 | 67 |
| 3.5 | 10 | ½ | 9 400 | 273 | 68 | 0.87 | 137 |
| 3.5 | 10 | 1 | 18 800 | 547 | 137 | 0.93 | 273 |

η is backward retrieval, optimal control. **Case A (15 µm) is not an efficient memory with single occupancy.** The manuscript's
sentence that "the large optical depth of a fibre strip allows efficient storage" holds for small cores only. [Paper text not changed;
needs the author's decision.]

Ways to raise OD: a smaller core (× $a^{-2}$); an LP₀ₘ signal (× g = 2.3–6.3, `hom_memory_theory.md`); a longer column (vacuum, loading);
f > 1 (collisions, light-assisted loss); a low-finesse cavity (fibre-end mirrors), which multiplies d by about the finesse (Gorshkov
2007: the cavity error is 1/C instead of 2.9/d).

---

## 2. Efficiency [script]

### 2.1 Optimal efficiency against OD

Retrieval kernel (Gorshkov et al., PRA 76, 033805, Eq. 15): $k_r(z,z') = \tfrac d2\, e^{-d(z+z')/2}\, I_0(d\sqrt{zz'})$. Its largest
eigenvalue is the optimal retrieval efficiency $\eta_r$. Optimal storage followed by backward retrieval gives $\eta_r^2$; forward retrieval gives the largest
eigenvalue of $K F K F$, where F flips z. Checks: flat spin wave against Eq. (30) agrees to 10⁻⁵; $(1-\eta_r)d$ → 2.72 at d = 400, approaching
Gorshkov's 2.9 from below.

| OD | 10 | 30 | 50 | 100 | 200 | 300 | 1000 |
|---|---|---|---|---|---|---|---|
| backward | 0.49 | 0.75 | 0.83 | 0.91 | 0.95 | 0.97 | 0.989 |
| forward | 0.39 | 0.60 | 0.70 | 0.80 | 0.88 | 0.91 | 0.968 |

Targets: η = 0.8 at OD = 41; 0.9 at 93; 0.95 at 201; 0.99 at 1096.

### 2.2 Backward retrieval and the hyperfine mismatch

Backward read-out of a forward-written spin wave leaves $\Delta k = 2\omega_{\rm hf}/c$, so $\Delta kL$ = 2.86 rad per cm.

| OD | L = 0.3 cm | 1 cm | 3 cm | no mismatch | forward |
|---|---|---|---|---|---|
| 50 | 0.825 | 0.768 | 0.510 | 0.831 | 0.695 |
| 100 | 0.903 | 0.872 | 0.710 | 0.906 | 0.801 |
| 300 | 0.964 | 0.954 | 0.892 | 0.965 | 0.911 |

At a fixed OD, a **short, dense** column is better. The optimal-control shape could recover part of the loss (not done here; Gorshkov
2007 Sec. VII treats re-optimisation). A backward spin wave has $k \approx 2k_{795}$, which matters for tunnelling (Section 3.5).

### 2.3 Control power

With $h = \int\Omega^2 dt/\gamma \approx 3d$, the photon number is $N_c \approx 3N\,r_s/r_c$, a few ×10⁴ per access, or **7–140 nW for 100 ns**.
Power is no limit. The question is control leakage into the detector and heating of the fibre, both small. [script]

---

## 3. Storage time [script unless tagged]

### 3.1 Model

$\eta(T) = \eta_0 \times$ motional dephasing (radial axes to the 4th power, the axial axis squared) × longitudinal inhomogeneity × magnetic term × $e^{-T/\tau_{\rm bg}}$.

- **Motional dephasing.** An atom in $|n\rangle$ gains phase $2\pi\epsilon n T$, with $\epsilon = \delta_{\rm diff}\,\nu/2U$. Lattice photon scattering (Rayleigh
  keeps the hyperfine coherence, Ozeri 2005) drives a ladder n → n ± 1 with rates $\gamma_\uparrow(n+1)$ and $\gamma_\uparrow n$. The coherent amplitude follows the ladder and rotates by $2\pi\epsilon n$.
  It starts in n = 0 (`axis_coherence`).
- **Longitudinal inhomogeneity.** A static Gaussian spread of the lattice intensity along the column, of rms dI, gives a differential-shift spread $\delta_{\rm diff}\,dI$. A spin echo removes it.
- Magnetic: not included (clock states). The background-gas lifetime $\tau_{\rm bg}$ is an input.

### 3.2 Lattice depth: case A geometry (a = 15 µm, ν_z = 300 kHz at u = 1)

| u | U (µK) | ν_z (kHz) | ν_r (Hz) | depth (lattice recoils) | J (Hz) | scattering (s⁻¹) | heating per radial axis (s⁻¹) | δ_diff (Hz) | ε_r (Hz) | ε_z (Hz) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 530 | 300 | 7150 | 5470 | 10⁻⁵⁸ | 4.2 | 0.47 | 760 | 0.24 | 10 |
| 0.1 | 53 | 95 | 2260 | 550 | 10⁻¹⁵ | 0.42 | 0.15 | 76 | 0.08 | 3.2 |
| 0.01 | 5.3 | 30 | 715 | 55 | 0.035 | 0.04 | 0.047 | 7.6 | 0.024 | 1.0 |
| 0.005 | 2.7 | 21 | 505 | 27 | 1.6 | 0.02 | 0.033 | 3.8 | 0.017 | 0.7 |

Heating in quanta falls only as $\sqrt u$, so lowering helps the shift ∝ u faster than the heating.
Gravity at u = 0.01: in a horizontal fibre the sag is 0.49 µm, and $mga$ = 1.5 µK against a 5.3 µK depth, so u ≈ 0.01 is close to the floor. In a vertical fibre the potential tilts by 1.1 kHz per site, a Wannier–Stark
ladder that suppresses tunnelling.

### 3.3 Scenarios (`memory_budget.png`, panel b)

| Scenario | $T_{1/e}$ | Limit |
|---|---|---|
| Full depth, 1 % intensity spread, no echo | **21 ms** | inhomogeneity |
| Full depth, echo | 0.9 s | radial heating (ε_r = 0.24 Hz, 0.47 quanta/s) |
| u = 0.01, 1 % spread, no echo | 1.3 s | inhomogeneity, background gas |
| u = 0.01, 0.1 % spread or echo, τ_bg = 2 s | 2.0 s | background gas |
| u = 0.01, echo, τ_bg = 10 s | 7 s | background gas and heating |
| u = 0.01, echo, τ_bg = 100 s | **11 s** | radial heating |

Literature for scale [lit]: 16 s storage with dynamical decoupling in a free-space 1D lattice (Dudin et al., PRA 87, 031801, 2013;
abstract and a DAMOP abstract only); 0.22 s at 76 % in a 3D lattice with a ring cavity (Yang et al. 2016); in hollow-core fibres, 20 ms 1/e
light storage at OD 40 (Leong et al., PRR 2, 033320, 2020), spin coherence 252 ms and trap lifetime 1.7–2.6 s at ≤ 10⁻⁸ Torr
(Xin et al., PRL 122, 163901, 2019), and Sr lifetimes of 350–500 ms along a fibre (Okaba et al., Nat. Commun. 5, 4096, 2014).

### 3.4 Longitudinal intensity: the first limit without an echo

Sources: fibre loss along the column, a focusing taper at the entrance, and **LP₀₂ admixture in the 1064 nm lattice**. The admixture beats with period
0.68 mm at a = 15 µm, with intensity amplitude $\pm 2\sqrt{pg}$: ±3 % for p = 10⁻⁴, that is ±23 Hz at full depth and ±0.23 Hz at u = 0.01.
Antiresonant fibres strip HOMs, which helps. This smooth, periodic phase is also the main **coherent crosstalk** source (Section 7.4).

### 3.5 Other terms

| Effect | Size | Comment |
|---|---|---|
| Tunnelling at u ≤ 0.005 | J ≈ 1.6 Hz | Harmless for a forward spin wave ($kd = 7.6\times10^{-5}$ rad). For a backward spin wave ($kd/2\pi$ = 0.34 mod 1), band dispersion dephases at ~J. Keep u ≳ 0.01 or use a vertical fibre. |
| Fibre bending modes | thermal amplitude 1–16 pm | Off resonance < 10⁻⁷ quanta/s; **if a mode coincides with ν_r**, 0.02–0.2 quanta/s. Avoid coincidence; this sets a clamping design rule. |
| Background gas | assumed $p = 3qL^2/8a^2\bar v$ | **The outgassing rate q is an assumption** (10⁻⁸ Pa m s⁻¹). Then τ_bg ∝ a²/L_fibre²: a = 3.5 µm and L = 10 mm give 0.15 s; L = 3 mm, 1.7 s; L = 1 mm, 15 s. OD²τ_bg ∝ 1/a² at a fixed q, so small cores still win, but **short fibre pieces** are needed. Needs measured data. |
| Magnetic noise | — | clock pair or magic field; the stretched scheme needs transfer to clock states for storage |
| Collisions (f > 1) | — | light-assisted loss during loading; density shifts |

### 3.6 What the radial ground state buys [estimate]

Without an echo, the ensemble coherence of motional dephasing saturates at $|\phi|^2 \to P_0^2$ per axis pair, so η has a ceiling $P_0^4$:
0.82, 0.89, 0.96 for $P_0$ = 0.95, 0.97, 0.99. **With an echo, static thermal dephasing is removed anyway,** so an honest comparison is with thermal atoms plus an echo:

- **Motional jumps defeat the echo.** Their rate per axis is $\gamma_\uparrow(2\bar n+1)$. Case A at 10 µK has $\bar n \approx 29$ (radial), so jumps are ≈ 60× more frequent than in the ground state.
  At u = 0.01 the phase per jump is small (ε_r = 0.024 Hz), so this matters at T ≳ 1 s and near the 11 s heating limit.
- **Shallow traps.** Adiabatic lowering keeps n. A 10 µK sample at u = 0.01 has $k_BT/U$ = 0.19: anharmonic, spilling out over seconds, and the trap must be
  deeper. A ground-state atom has $\hbar\nu/U \sim 10^{-2}$.
- **Coupling homogeneity.** At u = 0.01 thermal atoms spread over $x_{\rm rms}$ ≈ 2.2 µm, which reduces the on-axis OD by ≈ 20 % in a 15 µm core, more in small
  cores, and it destroys the g factor of an LP₀ₘ signal (`hom_memory_theory.md`).
- **Operation without an echo** removes pulse errors and makes random access simpler.

**Claim to test (S8):** at equal OD and depth, how much the ground state improves $\eta(T)$ over thermal atoms with an echo.

---

## 4. Noise and fidelity at the single-photon level [estimate]

- **Residual population in s** (imperfect optical pumping, p): about $d\,p$ incoherent photons per read in the guided mode. At d = 150 this needs p ≲ 10⁻⁵.
  The fibre helps: the guided fraction per atom is small, and the cold, trapped atoms can be pumped for a long time.
- **Four-wave mixing:** spin-wave probability ≈ $\tfrac32 d^2(\Gamma/\omega_{\rm hf})^2$: 2.4×10⁻² at d = 150, 0.27 at d = 500. It grows as d², so it limits high-OD operation.
  The **stretched scheme has no FWM channel**. Otherwise use a cavity, polarisation selection or phase mismatch (Lauk et al., PRA 88, 013823, 2013 for EIT fidelity with FWM).
- **Control leakage:** collinear geometry; filtering by polarisation (σ⁻ against σ⁺ in the magic-pair and stretched schemes), frequency (6.8 GHz etalons) and
  spatial mode (signal LP₀₂ against control LP₀₁: a mode filter; `hom_memory_theory.md`).
- Single-photon results in fibres for scale [lit]: SNR 59 at the single-photon level with cold atoms in HC-PCF (Peters et al., arXiv:1906.05771).

---

## 5. Bandwidth

Optimal adiabatic storage needs $T_{\rm pulse}\,d\gamma \gg 1$. With γ = Γ/2 = 1.8×10⁷ s⁻¹ and d = 35–70, $d\gamma$ ≈ 0.6–1.3×10⁹ s⁻¹, so pulses of ~10–20 ns, a
bandwidth of tens of MHz. No cavity linewidth applies. The time-bandwidth product of storage time × bandwidth would be ~10⁸ at 1–10 s. Raman storage at
large Δ keeps the efficiency and widens the bandwidth (Gorshkov 2007). [estimate]

---

## 6. Multimode operation

| Scheme | Capacity | Fibre implementation | Main issue |
|---|---|---|---|
| Temporal, EIT/Raman | ≈ $0.18\sqrt d$ modes at η > 0.7: 1, 2, 4, 7 at OD 100, 300, 1000, 3000 [script; Nunn 2008: ∝ √d] | kernel singular modes | small |
| Gradient echo (GEM/CRIB) | ∝ d (Nunn 2008) | axial magnetic gradient with the stretched scheme (1.4 MHz/G: 1 G/cm gives 1.4 MHz over 1 cm), then transfer to clock states | needs switching of the gradient; field-sensitive during write/read |
| Mode pairs: signal LP₀ₘ, control LP₀ₘ′ | M groups of M channels; spin wave $\kappa_{mm'} = (u_{m'}^2-u_m^2)/2ka^2$ | on-axis atoms couple to every LP₀ₘ | read in groups, crosstalk tail (Section 7.5), HOM guidance |
| Axial segments, frequency addressed | $N_{\rm seg}$, each with OD/$N_{\rm seg}$ | gradient plus microwave shelving into states the control does not couple | OD shared, so efficiency per channel falls |
| Polarisation | 2 | σ± dual rail | birefringence |

The fibre is one-dimensional: no transverse wavevector multiplexing as in free space (665 modes, Parniak et al., Nat. Commun. 8, 2140, 2017;
225 cells, Pu et al., Nat. Commun. 8, 15359, 2017). **The multimode case is not where this platform wins.**

**Random access and erasure.** All channels sit on the same atoms, and the control pumps any atom in s. Reading one spin wave therefore pumps the others
(Section 7.6): resonant EIT read-out erases them completely, and far-detuned Raman read-out keeps $\approx e^{-12.7/d}$ per access (about d/18 accesses
before half is lost). Shelving the waiting spin waves into a pair the control does not couple solves this at the cost of two transfer pulses per access.

---

## 7. Crosstalk between spin-wave modes: a bosonic model

**Literature status.** A short search found experiments that quantify crosstalk (Parniak 2017; Pu 2017; "Multiplexed storage and real-time manipulation based on a multiple degree-of-freedom quantum memory", Nat. Commun. 2018),
orthogonality as an overlap condition ("Spatial multiplexing in a cavity-enhanced quantum memory", arXiv:1903.02625), three-dimensional memory theory (Zeuthen et al., PRA 84, 043838, 2011;
Grodecka-Grad et al., PRL 109, 133601, 2012), singular-mode capacity (Nunn 2008), and phase mismatch for one mode (Gorshkov 2007). We did **not find a general
operator-level treatment that separates the crosstalk sources as below.** That is not a proof of novelty: a systematic search is item S5.0.

### 7.1 Operators

Atom j at $z_j$ with levels g, s. Spin-wave modes with mode functions $u_k$, normalised as $\sum_j |u_k(z_j)|^2 = 1$:

$$b_k = \sum_j u_k^*(z_j)\,\sigma_{gs}^{(j)},\qquad \sigma_{gs} = |g\rangle\langle s| .$$

At low excitation (Holstein–Primakoff), $[b_k, b_{k'}^\dagger] = \sum_j u_k^* u_{k'}\,(\sigma_{gg}-\sigma_{ss})^{(j)} \approx O_{kk'}$, the Gram matrix, with
corrections of order (excitations)/N. The modes are independent bosons only if O = 1.

**Write:** $b_k = \sum_l W_{kl}\,a_l^{\rm in} + \text{noise}$. **Storage:** $b_k(T) = \sum_{k'} M_{kk'}(T)\, b_{k'}(0) + F_k$, with

$$M_{kk'}(T) = \sum_j u_k^*(z_j)\,u_{k'}(z_j)\,\chi_j(T),\qquad \chi_j = e^{-i\phi_j(T)} ,$$

where $\phi_j$ is the phase of atom j from light shifts, fields and motion. The Langevin operator $F_k$ keeps the commutators; it is vacuum noise when no atom is in s (no added photons).
**Read** with control c (mode, frequency, timing): $a_{\rm out}^{(c)} = \sum_k R^{(c)}_k\, b_k + \text{vacuum}$. In the 1D Maxwell–Bloch model, $R^{(c)}$ comes
from the kernel: retrieving a spin wave $S(z)$ gives $\eta_r[S] = \langle S|K|S\rangle$.

**Crosstalk matrix:** $C_{ck} = |(R^{(c)} M)_k|^2$. Channel efficiency $\eta_c = C_{cc}$. The wrong-photon probability per read is $\sum_{k\neq c} C_{ck}\langle n_k\rangle$.
**Erasure** of spectators during a read: $b_k \to e^{-\Lambda_k} b_k$ + noise.

This separates four sources: O (geometry and filling), M (storage), R (read-out phase matching) and Λ (erasure).

### 7.2 Random filling: the Gram floor [script]

For plane waves $u_k = e^{ikz}/\sqrt N$ on the orthogonal grid $\Delta k = 2\pi m/L$: $O_{kk'} = 0$ exactly at f = 1. With random filling,
$E|O|^2 = (1-f)/N$ (Monte Carlo: 5.3×10⁻⁴ against 5.0×10⁻⁴ at N = 1000; 4.94×10⁻⁵ against 5.0×10⁻⁵ at N = 10⁴).
**Aliasing:** k and k + 2π/d are the same spin wave on the lattice. A forward spin wave sits at $kd/2\pi \approx 0$; a backward one at 0.34, so different channels
must not differ by a reciprocal lattice vector.

### 7.3 Independent per-atom phases [script]

If the $\phi_j$ are independent across atoms (motional state, local noise):

$$E|M_{kk'}|^2 = |\bar\chi|^2|O_{kk'}|^2 + (1-|\bar\chi|^2)\sum_j |u_k|^2|u_{k'}|^2 \;=\; (1-|\bar\chi|^2)/N \ \text{(orthogonal plane waves)}.$$

Monte Carlo with thermal-like $n_j$: $N|M_{kk'}|^2$ = 0.013, 0.172, 0.553 against the predicted 0.013, 0.180, 0.555. **The amplitude lost from a mode is spread over all N modes,
so the crosstalk into any one is negligible.** Motional dephasing is loss, not crosstalk.

### 7.4 Phases correlated along z [script]

For a smooth $\phi(z)$, M is the Fourier transform of $e^{-i\phi(z)}$ at $k-k'$. This is coherent transfer:
- a **periodic** phase $\phi_0\cos qz$ moves amplitude $J_m(\phi_0)$ into $k \pm mq$ (Jacobi–Anger);
- a **linear** gradient βz shifts every k by βT. A gradient echo memory uses this reversibly.

The lattice LP₀₂ beat (Section 3.4), with $\phi_0 = 2\pi\,\delta\nu\,T$, gives $J_1^2$ = 10⁻³ after (`crosstalk.png`, panel b):

| p | u = 1 | u = 0.01 |
|---|---|---|
| 10⁻³ | 0.14 ms | 14 ms |
| 10⁻⁴ | 0.43 ms | 43 ms |
| 10⁻⁵ | 1.4 ms | 0.14 s |

This becomes crosstalk only if channels are spaced by q (or mq). Otherwise it is loss, $1-J_0^2$. An echo removes it if it is static. Note that $q_{1064}$ for LP₀₁/LP₀₂ is 1.34× the
795 nm $\kappa$ of the same mode pair, so mode-pair channels and the lattice beat do not coincide by accident; check this for each design.

### 7.5 The read-out phase-matching tail [script]

Reading a spin wave with mismatch $\phi = \Delta kL$ into the output mode of channel c, with the optimal spin wave and the forward kernel:

| d | X/X₀ at 2π | 10π | φ = d | 2d | 4d | φ for 10⁻² | φ for 10⁻³ |
|---|---|---|---|---|---|---|---|
| 15 | 0.43 | 0.024 | 0.10 | 0.027 | 6.7×10⁻³ | 49 | 156 |
| 50 | 0.73 | 0.079 | 0.032 | 8.3×10⁻³ | 2.1×10⁻³ | 91 | 289 |
| 150 | 0.89 | 0.22 | 0.011 | 2.9×10⁻³ | 7.2×10⁻⁴ | 161 | 508 |

**Law, for φ ≳ d:** $X/X_0 \approx 1.6\,d/(\Delta kL)^2$, within 10 %. The EIT momentum window has a width of about d/L (Gorshkov 2007), so a dense medium reads out spin waves far from
phase matching. Channels need $\Delta kL \gg d$, not 2π.

**Mode pairs** (a = 15 µm, L = 1 cm, LP₀₁–LP₀₄): κL from 69 to 375 rad. A control LP₀ₘ′ reads every channel (x, m′) at once into LP₀ₓ, so these form groups. The nearest
channel of another group has a mismatch of 14 rad (for example control LP₀₁ and channel LP₀₄/LP₀₃ into LP₀₃), so X ≈ 1.6·d/200, about 3 % at d = 4.
Since $\kappa L \propto L/a^2$ and $d \propto L/a^2$, **X scales as $a^2/L$**. At a = 3.5 µm the minimum mismatch is 257 rad, and X ≈ 8×10⁻⁴ at d = 34. The κ = 0 spin wave (any m = m′) is read by every
control, so keep it empty.

### 7.6 Erasure during read-out [script]

A spectator spin wave not phase matched to the control sees only optical pumping: amplitude $e^{-h/(1+\delta^2)}$, δ = Δ/γ. Read-out of the target at
constant control (adiabatic Maxwell–Bloch, exact for piecewise-linear S; flat-spin-wave check against Eq. (30) agrees to 10⁻⁴) needs:

| OD | Δ/Γ | read efficiency | spectator survival |
|---|---|---|---|
| 30, 100, 300 | 0 (EIT) | 0.80, 0.89, 0.94 | ~0 (e⁻³⁰ or less) |
| 100 | 20 | 0.86 | 0.76 |
| 300 | 50 | 0.88 | 0.91 |
| 30, 100, 300 | ≫ d/2 (Raman) | 0.79, 0.87, 0.90 | 0.73, 0.80, 0.92 |

In the Raman limit, survival ≈ $e^{-12.7/d}$: the read needs $h \propto \delta^2/d$, so the loss per access is about const/d.
Hence the shelving of Section 6, or simultaneous read-out of all channels followed by sorting of the outputs (photonic lantern / multi-plane light converter).

### 7.7 Example budget [estimate]

a = 3.5 µm, 1 cm, f = ½, stretched scheme (OD 137), four mode-pair groups, Raman read-out with shelving. Gram floor 5×10⁻⁵; independent noise ≤ 10⁻⁴;
mode-pair tail ≈ 2×10⁻³ (d = 68, φ = 257 rad); LP₀₂ lattice beat ≲ 10⁻³ for T ≤ 40 ms at p = 10⁻⁴ and u = 0.01, unless an echo is used. **The read-out tail dominates.**

---

## 8. Against the state of the art

| Memory | Efficiency | Lifetime | Note |
|---|---|---|---|
| Hsiao 2018, cold Rb, EIT | 92 % | — | coherent pulses |
| Yang 2016, 3D lattice, ring cavity | 76 % | 0.22 s | best joint result found |
| Dudin 2013, 1D lattice | low (to check) | 16 s with DD | |
| Leong 2020, HC-PCF, cold Rb | — | 20 ms (1/e) | OD 40, 10 µK |
| Xin 2019, HC-PCF | coherence, not storage | 252 ms; trap 1.7–2.6 s | ≤ 10⁻⁸ Torr |
| **This platform, model** | 0.87–0.93 (a = 3.5 µm) | 2–11 s (τ_bg-dependent) | no cavity, fibre-coupled |

A real advance would be **η ≥ 0.9 with τ ≥ 1 s, no cavity, fibre-coupled.** The model allows it only with a small core (≤ 5 µm), a lowered lattice, an echo,
and τ_bg ≳ 10 s in a short fibre piece. The least certain input is the vacuum. The benchmarks use different protocols and definitions of efficiency.

---

## 9. Simulations for a future project (in priority order)

| # | Task | Method | Output | Effort |
|---|---|---|---|---|
| S1 | Efficiency with the real D1 structure: all F′ and m levels, off-resonant lines, FWM, optical pumping | 1D Maxwell–Bloch, multilevel, optimal control (Gorshkov iteration) | η(OD, Δ, scheme); stretched against clock | medium |
| S2 | Backward retrieval with mismatch, re-optimised shape | kernel with phase, iterative optimisation | η(L, OD) | light |
| S3 | Single-photon fidelity and SNR | Heisenberg–Langevin with FWM and residual population (Lauk 2013) | g², SNR against d | medium |
| S4 | η(T) Monte Carlo: Rayleigh and Raman scattering, recoil per axis, real lattice-intensity map from the fibre mode content, echo/DD with pulse errors, magnetic noise | trajectories per atom; input P(n) from the cooling project | η(T) curves | medium |
| S5.0 | **Literature search on crosstalk models** | — | novelty check | light |
| S5 | Crosstalk: discrete atoms with random filling; temporal, mode-pair and GEM channels; M from S4; R from S1 | the bosonic map of Section 7, numerically | crosstalk matrices; check the laws (1−f)/N, $J_m$, 1.6d/φ², $e^{-12.7/d}$ | medium |
| S6 | Vacuum in the core: molecular flow with measured outgassing; side holes | analytic plus literature data | τ_bg(a, L) | light |
| S7 | Design optimisation over a, column length, f, scheme, depth | combine S1–S6; cooling time and P₀ from the cooling project | η(T), η·τ, a repeater figure of merit | medium |
| S8 | **Value of the ground state:** thermal atoms with an echo against the ground state, at equal OD and depth | S4 with a thermal start; anharmonic trap; coupling inhomogeneity | the claim of Section 3.6 | light–medium |
| S9 | f > 1: collisions, light-assisted loss, density shifts | literature plus rate model | OD gain against loss | light |
| S10 | Low-finesse cavity on the fibre ends | cavity model (Gorshkov paper I) | η against finesse | light |

---

## 10. Starting a new project

- Copy this folder. The scripts need only numpy, scipy, sympy and matplotlib; constants come from `scipy.constants`, Rb data from Steck rev. 2.3.4.
- Treat every [estimate] as a hypothesis. The first checks are S5.0 and S2 (cheap), then S1 (the scheme choice and the stretched scheme).
- Inputs from the cooling project are listed in Section 1.3. Ask for the current values: they change as the cooling study converges.
- Conventions: d = OD/2 (Gorshkov); η is total (storage × retrieval); depth fraction u is relative to the depth that gives ν_z = 300 kHz.

## 11. References

Verified on Crossref (DOI, title) in this or earlier sessions:
Gorshkov et al., PRA 76, 033804 and 033805 (2007) [paper II read: Eqs. (15), (30), 2.9/d, phase mismatch]; Gorshkov et al., PRL 98, 123601 (2007);
Nunn et al., PRL 101, 260502 (2008); Zeuthen et al., PRA 84, 043838 (2011); Grodecka-Grad et al., PRL 109, 133601 (2012); Lauk et al., PRA 88, 013823 (2013);
Dudin et al., PRA 87, 031801 (2013); Dudin et al., PRA 81, 041805 (2010); Zhao et al., Nat. Phys. 5, 100 (2009); Yang et al., Nat. Photon. 10, 381 (2016);
Hsiao et al., PRL 120, 183602 (2018); Bao et al., Nat. Phys. 8, 517 (2012); Okaba et al., Nat. Commun. 5, 4096 (2014) [lifetime figure read];
Leong et al., PRR 2, 033320 (2020) [read: OD 40, 20 ms, 10 µK]; Xin et al., PRL 122, 163901 (2019) [read: 252 ms, 1.7–2.6 s, ≤ 10⁻⁸ Torr];
Langbecker et al., NJP 20, 083038 (2018); Parniak et al., Nat. Commun. 8, 2140 (2017); Pu et al., Nat. Commun. 8, 15359 (2017);
Mazelanik et al., npj QI 5, 22 (2019); Treutlein et al., PRL 92, 203005 (2004); Savard et al., PRA 56, R1095 (1997); Derevianko, PRL 105, 033002 (2010);
Sagi et al., PRL 105, 053201 (2010); Ozeri et al., PRL 95, 030403 (2005); Heshami et al., J. Mod. Opt. 63, 2005 (2016).

Not on Crossref under this title: Peters et al., arXiv:1906.05771 (read: η = 0.36 at 200 ns, d ≈ 80–100, SNR 59).
To verify (authors, details): the 2018 Nat. Commun. multi-degree-of-freedom memory (s41467-018-05669-5) and arXiv:1903.02625. Outgassing rate q and the fibre-clamp Q are assumptions.
