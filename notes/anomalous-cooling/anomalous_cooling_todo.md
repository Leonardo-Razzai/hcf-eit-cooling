# Anomalous (scale-free) Fock statistics in radial EIT cooling — research to-do

Possible separate paper. **Theory note: `anomalous_cooling.tex`** (figures from `make_plots.py`,
run from this folder). Status: **hypothesis**, from the case-A convergence problem
(branch `fix-truncation`, November 2026). Nothing here is established yet.

## Observation that started it

- Case A ($a = 15$ µm, $\eta_{\rm sp} = 0.71$) needs very large Fock truncations N; case B ($\eta_{\rm sp} = 0.34$) is
  converged at $N = 20$. The steady state has a long recoil-diffusion tail in P(n).
- Weak-probe N-scan (banded eliminated ME, B = 4, `results/convergence/nconv_weak.dat`):
  - $s = 0.3$: converged by $N \approx 250$ (< 0.1 %).
  - $s = 0.01$: +0.1 % per 50 states at $N = 300$ ($n_{\rm ss} \approx 0.0180$).
  - $s = 0.1$: still +0.65 % per 50 states at $N = 300$ ($n_{\rm ss} \approx 0.0324$ → extrapolated 0.0329);
    increments shrink ×0.67 per 50 states between $N \approx 225$ and 275.
  - W converges much faster than $n_{\rm ss}$ (W is measured from $n = 1$).
- Convergence is worst at intermediate s (~0.1), better at both smaller and larger s.

## Minimal model (to be checked)

One two-level system per Fock state: |$g_p$,n⟩ ↔ |+,n−1⟩ with Raman coupling $\Omega_R\sqrt n$, narrow-line
width $\Gamma_n$ ≈ ($\Gamma/\Delta$)$\nu_r$, effective saturation **$s_n = sn$**.

- Removal rate $R(n) \approx R_{\rm sat}\, sn/(1 + csn)$, $R_{\rm sat} \propto \Gamma_n/(2\kappa)$, $\kappa = 1/b_p$, c ≈ 2.25.
- Drift $v(n) = -(1 - h) R(n)$, with $h = \kappa\alpha\eta_{\rm sp}^2$ (recoil quanta per removed phonon).
- Recoil energy diffusion $D(n) \approx hnR(n)$.
- Fokker–Planck: $P(n) \propto \exp(\int v/D\,dn)/D$ → if $v/D \approx -(1-h)/(hn)$: $P(n) \propto n^{-1/h}$,
  exponent independent of s at this order. Case A $h \approx 0.40$ → P ~ $n^{-2.5}$; case B $h \approx 0.09 \to n^{-11}$.
- Crude: recoil jumps are large and discrete; the exponent may be modified.

## Hypothesis: crossover in s, transition in h

- **$sn = 1$ is a crossover** per Fock state (linear → saturated removal); s sets where the tail starts
  (n* ≈ 1/s) and its weight, which explains the worst convergence at intermediate s.
- **A genuine transition only in h, in the limit N → ∞** (the unbounded ladder is the thermodynamic
  limit, like the photon number at a laser threshold).

Features that would characterise it:
1. Scale-free (power-law) tail $P(n) \propto n^{-1/h}$ instead of exponential.
2. Successive divergence of moments: $h < 1/3$ all finite; 1/3 ≤ $h < 1/2$ Var(n) diverges;
   1/2 ≤ $h < 1$ normalisable, $P_0$ finite, $\langle n\rangle$ = ∞ ("anomalous cooling"); $h \geq 1$ no steady state,
   $P_0$ → 0. Two thresholds: **$h = 1/2$ (mean energy), $h = 1$ (existence of a steady state)** —
   stricter than the naive $h < 1$. Case A ($h \approx 0.40$) would be in the "finite $\langle n\rangle$, divergent Var" window.
3. Finite-size scaling: missing $\langle n\rangle \propto N^{2-1/h}$; N acts as system size; convergence exponent
   measures h. First rough estimate from $s = 0.1$: $h_{\rm eff} \approx 0.33$ (nominal 0.40).
4. Universality: exponent depends only on h; tails for different ($\nu_r/\omega_{\rm rec}$, $b_p$, s) at equal h
   should collapse.
5. Critical slowing down: non-exponential (power-law) relaxation of hot states; Liouvillian gap
   closing with N near threshold; W (from $n = 1$) stays finite.
6. Rounding by physical cut-offs (act as finite size), each bringing a new parameter:

   | scale | condition | case A, $s = 0.1$ | parameter |
   |---|---|---|---|
   | saturation | s $n = 1$ | $n \approx 10$ | s |
   | loss of sideband resolution (blue sideband resonant; high-intensity EIT breakdown) | s $n \approx 2 (\nu_r/\Gamma_n)^2$ ≈ 2($\Delta/\Gamma$)² | $n \approx 180$ | $\Delta/\Gamma$ |
   | Bessel non-linearity (capture range) | $u_{11}$·2$x_0\sqrt n$/a ≈ 1 | $n \approx 470$ | $x_0$/a (core size) |
   | trap anharmonicity shift ~ $\Gamma_n$ | (3/8)(n $\hbar\nu_r$/$U_0$)$\nu_r \approx \Gamma_n$ | $n \approx 1380$ | $U_0/\hbar\nu_r$ (lattice depth) |
   | trap depth (loss) | n $\hbar\nu_r$ = $U_0$ | $n \approx 1550$ | $U_0$ |

Consequences: W insensitive to the tail; $n_{\rm ss}$ = core ($\Gamma/4\Delta$)²(1+h)/(1−h)(1+s) + tail term set by h
(exponent) and s, $\Delta/\Gamma$, core size, lattice depth (cut-offs); $n_{\rm ss}$ very sensitive near h → 1/2.
$P_0$ is the robust figure of merit.

Physical analogy (verify on Crossref before citing): Lévy statistics in subrecoil laser cooling
(VSCPT, Raman cooling; Bardou, Bouchaud, Aspect, Cohen-Tannoudji) — dark state plus random recoil →
broad, non-ergodic statistics; random multiplicative ("Kesten") processes with power-law
stationary tails. Related high-intensity EIT breakdown: Bartolotta et al., arXiv:2411.18818 (2024);
Morillas-Rozas et al., arXiv:2607.22243 (2026, preprint).

## First diagnostics (7 Nov 2026, `results/convergence/tail_mechanism.png`)

Banded eliminated ME, case A (h = κα $\eta_{\rm sp}^2 = 0.40$), drift and diffusion from Fock states
(`rates_vs_n.dat`), tail at $s = 0.1$ (`tail_s0p1.dat`):

- **Scale-free tail confirmed at $s = 0.1$:** P(n) is a clean power law over $n = 30$–200,
  d ln P / d ln $n = -2.31$ … −2.39 (log-linear slope varies 4×, so not exponential). Prediction
  −1/h = −2.48. Steepens to −2.9 above $n \approx 200$. $N = 300$ and 400 give identical tails;
  $n_{\rm ss} = 0.03243$ → 0.03265 (+0.7 %), $W = 335$ → 338 s⁻¹.
- **Saturation crossover confirmed:** removal rate grows $\propto$ n up to $n \approx 1$/s and saturates at
  ≈ 2000 s⁻¹ ($s = 0.01$: linear to $n \approx 50$–80, plateau from $n \approx 120$; $s = 0.1$: plateau from $n \approx 10$–20;
  $s = 1$: plateau from $n = 1$). Diffusion per quantum saturates at ≈ 2800 s⁻¹, i.e. D $\propto$ n R(n) on
  the plateau, as assumed.
- **Unexpected: no collapse at loss of sideband resolution — a fast high-n cooling channel.**
  Beyond $sn \approx 20$ (≈ 2($\Delta/\Gamma$)² = 18, the predicted resolution-loss scale) the removal rate grows again,
  roughly $\propto$ n, far above the narrow-line bound ($s = 1$: 4.6×10³ s⁻¹ at $n = 30$, 4.4×10⁴ s⁻¹ at $n = 200$;
  $s = 0.1$: rising from $n \approx 120$). This channel cuts off the power-law tail. Hierarchy:
  linear ($n < 1$/s) → saturated plateau with power-law tail (1/s < n < ~20/s) → fast cooling
  (n > ~20/s). The plateau spans $n \approx 10$–200 at $s = 0.1$ (inside the computed range → slowest
  convergence), only $n \approx 1$–20 at $s = 1$ (tail cut early → fast convergence), and $n \approx 100$–2000 at
  $s = 0.01$ (tail weight tiny, but still drifting +0.1 %/50 states). Mechanism to be understood
  (hot atoms leave the narrow EIT resonance; Doppler-like cooling on the broad line?).
- **Fokker–Planck is only qualitative:** the measured local ratio −v n/D is 0.5–1.2 on the plateau,
  not $$(1-h)/h = 1.48$$; the large discrete recoil jumps matter. The tail exponent nevertheless is
  close to −1/h.

## The fast high-n channel: Sisyphus-like (gray-molasses) cooling? (7 Nov 2026)

Figure: `results/convergence/regimes.png` (from `rates_vs_n.dat`). Plotted against $sn_0$ the drift
curves **collapse**: $s = 0.01$ and 0.1 coincide in the resolved and plateau regions; $s = 0.1$ and 1
coincide at the onset of the fast channel, at $sn \approx 2 (\nu_r/\Gamma_n)^2$ = 18.0 ($\Gamma_n/\nu_r$ = 0.333 for
$\Delta = 3\Gamma$). Three regimes of one trapped atom:

| regime | condition | removal rate | limited by |
|---|---|---|---|
| resolved EIT (quantum, sideband) | $sn < 1$ | R $\propto$ n | narrow line, unsaturated |
| saturated plateau (power-law tail lives here) | 1 < $sn < 2 (\nu_r/\Gamma_n)^2$ | R ≈ const ≈ 2×10³ s⁻¹ | narrow line $\Gamma_n$ |
| fast, Sisyphus-like | $sn > 2 (\nu_r/\Gamma_n)^2$, i.e. $\Omega_R\sqrt n$ > $\nu_r$ | R $\propto$ n, far above $\Gamma_n$ bound | optical pumping / light shift |

**Onset, three equivalent forms (derived in the note, Eq. onset).** With $\tan\theta(x)$ =
($\Gamma_n/\nu_r$)$\sqrt{s/8}\,x/x_0$ along the LP₁₁ profile and the turning point $x_t = 2$$x_0\sqrt n$:
$\Omega_R\sqrt n$ = $\nu_r$ ⇔ $\tan\theta(x_t)$ = 1 (dark state rotated by π/4 over the orbit) ⇔ bright-state light-shift
hill $\delta_B(x_t) - \delta_B(0)$ = $\hbar\nu_r$ ⇔ $sn = 2(\nu_r/\Gamma_n)^2$. In general the hill is $\hbar\nu_r$ $\tan^2\theta(x_t)$ $\propto$ n, so
regime III removes $\propto$ n quanta per cycle ($s = 1$, $n = 200$: ~11 $\hbar\nu_r$ vs measured ~22× plateau rate).
The rate model with the calibrated c = 2.25 reproduces regimes I and II of the measured drift.

**Intuitive picture (hypothesis).** When the Raman coupling at level n exceeds the trap frequency
($\Omega_R\sqrt n$ > $\nu_r$) the internal state evolves faster than the motion: the atom follows the local dark
state $|D(x)\rangle \propto \Omega_c|g_p\rangle - \bar\Omega_p f_p(x)|g_c\rangle$, whose mixing angle $\theta(x) \approx \bar\Omega_p Gx/(a\Omega_c)$ varies along
the LP₁₁ gradient. Motion couples it non-adiabatically ($\propto$ velocity) to the bright state, which for
blue detuning ($\Delta > 0$) is light-shifted up: the atom climbs a potential hill before optical pumping
returns it to the dark state. That is dark-state Sisyphus / gray-molasses cooling (friction $\propto$ v², so
energy removal $\propto$ E $\propto$ n). It also explains why the GRM ($c_{\rm sat}$ calibrated for $n \leq 15$) underestimates
the ME drift at $n_{0} = 20$–30 in case A (GRM/ME = 0.42 at $n_{0} = 30$, $s = 1$).

**Literature (verified on Crossref/arXiv; full texts not yet read).** Known separately:
- gray molasses: Grynberg & Courtois, EPL 27, 41 (1994); Weidemüller et al., EPL 27, 109 (1994);
- unified theory of cooling of a trapped atom (resolved/unresolved sideband, gray molasses,
  $\Lambda$-enhanced gray molasses): Phatak et al., PRA 110, 043116 (2024), arXiv:2406.19153;
- EIT sub-Doppler friction for free atoms: He, Tengdin, Anderson & Rey, PRA 95, 053403 (2017);
- high-energy breakdown / capture range of EIT cooling for ions: Bartolotta et al., arXiv:2411.18818;
- $\Lambda$-enhanced gray molasses in tweezers: Janse van Rensburg et al., PRA 113, 033123 (2026).
Not found: the crossover resolved EIT → plateau → Sisyphus within one trapped atom vs n and s, in a
gradient-only (node) geometry. To check in full texts before claiming novelty.

**Fibre-specific aspects.** The mechanism is general $\Lambda$ physics; the regime is set by the fibre:
(i) the Sisyphus "lattice" is the transverse mode profile (scale $a \approx 15$ µm), not $\lambda$/2 — a smooth hill
over the whole radial excursion, large capture range, no corrugation; (ii) the soft trap puts the
crossover inside the thermal population (n̄ ≈ 30 at 10 µK): hot atoms are captured in the
Sisyphus regime and finished in the EIT regime (in tight tweezers/ion traps the crossover is at much
larger n); (iii) the node keeps the low-n regime carrier-free. Consequence: hot-atom capture may be
much faster than the GRM-based capture times of the paper.

## To do

- [x] Analyse the running diagnostics (`studies/convergence.py rates tailshape`) — see above:
      drift v(n), diffusion D(n) vs $n_0$ at $s = 0.01$, 0.1, 1 (`rates_vs_n.dat`); P(n) at $s = 0.1$,
      $N = 300/400$ (`tail_s0p1.dat`); figure `results/convergence/tail_mechanism.png`.
      Check: R(n) saturation at $n \approx 1$/s; D $\propto$ n R; local ratio −v n/D ≈ $$(1-h)/h$$, s-independent;
      log–log tail slope vs −1/h.
- [ ] Write the minimal analytic model properly (R(n), v, D, tail exponent, crossover scales,
      tail contribution to $n_{\rm ss}$); improve the Fokker–Planck treatment for large discrete jumps
      (e.g. master equation for P(n) with the exact recoil kernel).
- [ ] Scan h via $\nu_r/\omega_{\rm rec}$ and $b_p$ with the banded solver (reduced cases, branch `universal`):
      tail exponents, finite-size-scaling exponents, collapse in h, thresholds $h = 1/2$ and $h = 1$.
- [ ] Fast channel — test the Sisyphus picture:
  - [ ] map R(n, s) on a 2D grid ($n \approx 1$–300, $s = 0.01$–10; banded drift, ~6 min per s); collapse vs
        $sn$ and vs $\Omega_R\sqrt n/\nu_r$;
  - [ ] vary $\Delta/\Gamma$: onset should move as 2$(\nu_r/\Gamma_n)^2$/s ≈ 2($\Delta/\Gamma$)²/s;
  - [ ] flip the sign of $\Delta$: gray molasses needs blue detuning, so red detuning should turn the
        high-n regime into heating while the resolved regime survives;
  - [ ] scaling of the high-n rate with $\Omega_c$, $\Delta$, $\Gamma$ (Sisyphus friction, not $\Gamma_n$); bright-state
        population and dark-state following at high n;
  - [ ] semiclassical friction/diffusion of an atom moving through the LP₁₁ gradient (linear
        response of the optical Bloch equations) vs ME drift at high n;
  - [ ] consequence: capture from thermal states with banded time evolution (n̄ ≈ 30, $N \approx 400$)
        instead of the GRM — true capture time from 10 µK including this channel.
- [ ] Optional: Liouvillian gap vs N near threshold; time-domain relaxation from hot states
      (power-law vs exponential).
- [ ] Role of the cut-offs: include anharmonicity / trap depth / Bessel non-linearity and see how
      they round the transition; dependence on $\Delta/\Gamma$ and core size.
- [x] Literature check, done 9 Nov 2026: see below and note Sec. "Relation to earlier work". Lévy statistics in subrecoil cooling, Kesten processes,
      dissipative phase transitions in open quantum systems, high-intensity EIT cooling breakdown;
      full texts of Phatak 2024, Bartolotta 2024, He 2017 for the EIT/gray-molasses crossover.
- [ ] Decide: section of the fibre paper, or a separate paper on anomalous Fock statistics in
      laser cooling with soft traps.

## Observations from the D1/D2 convergence checks (November 2026, to follow up)

Three-level case A, `results/convergence/lines_nconv.dat`, `lines_nconv_delta.dat`.

- **$h > 1$ does not always mean a runaway.** At $b_p = 1/6$ ($h = 1.21$) $n_{\rm ss}$ runs away with N at $s = 0.1$
  (0.83, 1.88, 4.18 at $N = 60$, 120, 250), but it converges at $s = 1$ (2.69) and $s = 5$ (1.44) by $N = 120$.
  Possible reading, not yet tested: the fast channel above $sn = 2(\nu_r/\Gamma_n)^2$ ($n = 18$ at $s = 1$,
  180 at $s = 0.1$) cuts off the tail, so the thresholds in h would hold only where this channel is out of
  reach. To do: extend the $s = 0.1$ scan beyond $N = 250$, and compare the N where $n_{\rm ss}$ saturates with 18/s.
- **The tail reaches further at large $\Delta$.** $b_p = 1/2$, $s = 1$: $n_{\rm ss}$ and W are converged by $N = 120$ at
  $\Delta$ <= $3\Gamma$, but $N = 250$ -> 400 still changes W by 1.3 % and 2 %, $n_{\rm ss}$ by 0.7 % and 3 %, at 8 and
  $12\Gamma$. Consistent with the same cut-off moving out as $\Delta^2$ ($n = 288$ at $12\Gamma$, $s = 1$). To do:
  test this scaling in the $\Delta/\Gamma$ scan of the onset.

## Theory of the h windows and figures of merit (8 Nov 2026)

Written up in `anomalous_cooling.tex`, Sec. "Theory of the h regimes and figures of merit", figure
`figures/fom.pdf`. No new simulations: tests use stored data. Three thresholds, four windows.

Findings:
- **Correction.** The earlier "$-vn/D = 0.5$–1.2, Fokker-Planck only qualitative" used $D = d{\rm Var}/dt$.
  With the model's $D = \tfrac12\,d{\rm Var}/dt$ the ratio is 1.0-1.9 on the plateau, around $(1-h)/h = 1.48$, and the
  Fokker-Planck form built from the measured v and D reproduces the tail for n >~ 5 (slope 30 <= n <= 200:
  -2.17 vs ME -2.33 at $s = 0.1$). It fails below n ~ 5, where a single jump is comparable with n.
- **Regime III: more energy per photon, not more photons.** $s = 1$, $n = 20$ -> 200: drift x17, $D/n$ only x1.8.
  Model: $P(n) \sim n^{-1/h}e^{-n/n_x}$, $n_x = hn_c$. Right form, scale about 2x too small
  ($n_x \sim 16$ at $s = 1$, against 7). Sharp onset near $n_c$ not captured. To be derived semiclassically.
- **Thresholds are crossovers,** width $\delta h \sim h^2/\ln[2(\Delta/\Gamma)^2]$: 0.09 at $h = 1/2$, $3\Gamma$; 0.04 at
  $12\Gamma$. Variance check: ${\rm Var} \sim N^{0.54}$ at $s = 0.1$, predicted $N^{0.5}$.
- **$h > 1$: hot steady state, not a runaway (hypothesis).** $\langle n\rangle = (h-1)n_c$. $b_p = 1/6$: 3.8 vs 2.69 at $s = 1$,
  but 0.76 vs 1.44 at $s = 5$; $s = 0.1$ predicts ~38 (unconverged 4.18 at $N = 250$). Scaling not confirmed.
- **Effective h.** Neighbouring levels change $b_p$: D2 F'=2, gp = 1, $12\Gamma$ has beff = 0.37, h ~ 0.57.
- **W fails in a heavy tail.** From $n = 1$, $\langle n\rangle$ undershoots and rises as the tail fills, so
  $1/W \sim 1/W_c - n_{\rm tail}/\lambda$: pole, then $W < 0$. Explains the pole vs $\delta_{\rm e}$, $W < 0$ at $\Gamma$, and
  $W = 206$ -> 1254 /s with N at D2 F'=2, $12\Gamma$.
- **Ground-state rate $W_0$ (from $P_0$(t), `results/convergence/rate_p0.dat`) is stable.** Case A, $s = 1$:
  $W_0 = 2188$, 2220, 2225 /s at $N = 20$, 40, 120 (W: 1319, 1544, 1614). D2 F'=2, $12\Gamma$: 507-511 /s while W
  runs away. Smooth maximum vs $\delta_{\rm e}$ across the pole of W. **Correction:** I first expected the true rate
  below W (~1200 /s); $W_0$ is above W, by 1.2x even at weak probe, so the two measure different things. Open:
  origin of that factor; is $1/W_0$ the time to reach $P_0$ from a thermal state?
- **Figure of merit.** Coherence $\lvert\phi(T)\rvert^2$, phi = characteristic function of P(n). Long T: ~ $P_0^{2}$, stays.
  Memory: $P_0$ (threshold $h = 1$). Interferometer/clock: $\sigma_n$ at fixed $P_0$ ($h = 1/3$), $T_{\rm QPN} = 0.53$ s ($s = 0.1$),
  0.36 s ($s = 1$) for case A, $\delta_{\rm diff}$ = 800 Hz, no echo. Mean shifts: $\langle n\rangle$ ($h = 1/2$). Overall:
  $t_{\rm prep}$($$$P_0$$^*$). $\langle n\rangle$ or W alone mislead.
- **Best configuration (model).** $W(\varepsilon) \sim \nu_r\sqrt\varepsilon\, b_p (1-h)^{3/2}(1+h)^{-1/2}\, s(1+s)^{-3/2}$. Candidate:
  D1 F'=1, probe on F=2 ($b_p = 5/6$, $h = 0.24$), $\Delta$ ~ 5-$8\Gamma$, s ~ 1 in the last stage. Main risk: Zeeman
  dark states of F=2 -> F'=1, not in the model. Ramps of s and a spill of the tail are untested ideas.

Next, when simulations resume (cost order):
- [x] $W_0$ at case A, $s = 0.01$–5, $b_p = 5/6$, D2 points (`rate_p0.dat`).
- [x] Select $\delta_{\rm e}$ on $W_0$ in the D1/D2 study (`lines_ml_A`, `lines_delta_A`). $P_0 = 0.99$: D1 F'=1 at $5\Gamma$,
      $W_0$ = 1.9e3 /s; $b_p = 1/2$ near $6\Gamma$, ~1.0e3 /s. $W_0$ meaningful only for $P_0$ >~ 0.85.
- [ ] Report $\langle n\rangle$/(1-$P_0$) in all scans.
- [ ] Chain for P(n) with the exact recoil kernel and measured R(n), N up to several thousand, regime III
      on and off: sharp thresholds, then their rounding.
- [ ] ME at $h = 0.1$–1.2 by two routes ($\nu_r/\omega_{\rm rec}$; $b_p$), $s = 1$, $\Delta = 8\Gamma$, $N = 400$: universality.
      Recheck the band in regime III.
- [ ] Banded evolution for $t_{\rm prep}$($$$P_0$$^*$) with $P_0$(t) output: ($\Delta$, s, $b_p$) and ramps.
- [ ] Derive $n_x$ semiclassically.

## Optimal s ramp: deferred simulations (9 Nov 2026)

Theory in `s_ramp_theory.md`. **Do not start these unprompted: ask the author first,** after the map and the band checks.

- [ ] Drift and diffusion from Fock states at $s = 3$, 5, 10 (`convergence.py rates` extended, ~6 min per s): does the friction peak?
- [ ] Master equation at constant $s = 3$, 5 from $\bar n_0 = 30$: tests the 1/s estimate of the hot stage.
- [ ] Chain for P(n) with the exact kernel and tabulated R(n,s); solve the Pontryagin problem by forward-backward sweeps.
- [ ] Validate the ramp as piecewise-constant s in the banded ME; write it as a recipe in $\nu_r/\omega_{\rm rec}$, $b_p$, $\Delta/\Gamma$, $\bar n_0$, $$$P_0$$^*$.

## Literature check (9 Nov 2026): result

**Not new:** power-law energy tail, exponent 1/parameter, mean diverging at 1/2, no stationary state above 1, for Sisyphus cooling in a
harmonic trap: Dechant, Kessler, Barkai PRL 115, 173006 (2015); Dechant et al. PRE 94, 022151 (2016), D <-> h. Same reason: constant energy
removal, diffusion ~ E. Also carrier-free EIT at a node for ions: Xing & Mehta PR Applied 24, 014034 (2025), experiment arXiv:2512.10900.
Saturated n-independent rate: Morillas-Rozas et al., arXiv:2607.22243.

**Possibly new (targeted search only):** tail in resolved-sideband ground-state cooling, quantum Fock treatment, h set by the trap and $b_p$; coexists with
$P_0 > 0.95$, begins near n ~ 3 (Sisyphus: > $10^3E_r$); fast-channel cut-off, threshold widths, hot state above $h = 1$; failure of W, use of $W_0$; memory
consequences; LP11 radial EIT in a hollow-core fibre.

**Still to read:** Roghani & Helm 2008 full text (may contain Fock distributions); Waterloo transverse-cooling result in HCF (Anderson et al., DAMOP 2024);
He 2017 full text; last 13 % of Bartolotta. Consider contacting the Barkai group.

## Map in ($b_p$, $\nu_r/\omega_{\rm rec}$) (9 Nov 2026)

`results/universal_map.dat`, plots in `results/universal/` (`map`, `collapse`, `target`). $\nu_r/\omega_{\rm rec} = 1$–15, $b_p = 1/6$–5/6, $\Delta = 3$, $6\Gamma$, $s = 0.3$, 1, 3.

- **Partial collapse in h.** $\langle n\rangle/(1-P_0)$ at equal h by different routes agrees within about 1-30 % for h >~ 0.4 (3.33 vs 3.39 at $h = 0.40$, $3\Gamma$), less well at
  small h, where it measures the core (1.97 vs 1.36 at $h = 0.16$). $1-P_0$ and $W_0$ do not collapse (factors 2-4): they depend on $b_p$ and $\nu_r$ separately.
- **No sharp threshold,** a smooth rise through $h = 1/3$, 1/2, 1, steeper at $6\Gamma$. Consistent with crossovers; the map does not resolve the widths.
- Closed form with $n_x = 2$ h $n_c$ is low by about 1.5-2x at h ~ 0.3-0.8, and high for $b_p = 1/6$, $h > 1$.
- N-sensitivity of $n_{\rm ss}$ appears only at large h and $6\Gamma$ (up to +44 %), near zero at $h < 1/3$.
- Band: $P_0$ stable (B = 4 vs 8, 3e-4 or better for $b_p$ >= 1/2); $W_0$ within 0.6 %, except $b_p = 1/6$ corner (6 %, $P_0 < 0.55$). Thermal $\bar n_0 = 30$: $P_0$ within 0.1 %, $\langle n\rangle$ within 3 % (B = 4 vs 6).

## Fast channel: derived, not only hypothesised (9 Nov 2026)

Note, "Derivation of the fast channel". The eliminated H is gray molasses by construction. Golden rule on $\dot\theta$ gives exactly $ns\Gamma_n/2$ (the sideband rate);
strong coupling gives $R = R_{\rm sat}(1 + cn/n_c)$, measured c = 1.3-1.9 for $n = 50$–200, $s = 1$. So sign test, $\Delta$ scan and semiclassical calculation are now confirmations, not tests of the picture.
Sign flip predicts a capture range, $R_{\rm sat}(1 - cn/n_c)$.

Open: dip at $sn \sim 8$ (0.69-0.75 of the base rate, collapses in s n); value of c; onset shape; high-s limit. Cheap confirmation if wanted: one Fock-drift scan at $\Delta < 0$, one at another $\Delta$.

