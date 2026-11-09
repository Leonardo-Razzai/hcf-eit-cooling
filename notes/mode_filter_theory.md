# The atomic column as a mode filter for the probe: theory (9 Nov 2026)

Analytic model, using stored scattering rates. Nothing here has been simulated with propagation. Numbers: closed $\Lambda$, one radial axis, $\Delta = 3\Gamma$, $s = 1$, case A unless stated. $d$ is the two-level optical depth of the column in LP$_{01}$ on axis.

## Question

Atoms on the axis see the LP$_{01}$ part of the probe with the full optical depth, and LP$_{11}$ only through its gradient. Does the column clean the probe, $y(z) \equiv \lvert\epsilon(z)\rvert^2$, the LP$_{01}$ power fraction?

## Why EIT reduces the filter

The LP$_{01}$ part is the carrier, at two-photon detuning $\delta_{\rm e} \approx 0$, the transparency point. An atom in the dark state with the coupling beam does not absorb it. Photons are removed only when the atom has to re-form its dark state.

## Rate of absorption per atom

With the carrier, the dark state rotates by $\theta \simeq \epsilon\bar\Omega_p/\Omega_c$: a fraction $\theta^2$ of the population lies in $g_c$. Each scattering event at rate $\gamma_p$ (photons per atom per time, from any cause) resets the atom, and re-forming the dark state takes $\theta^2$ photons from LP$_{01}$ into the coupling beam. So

$$\dot N_{01} = \gamma_p\,\theta^2 = \gamma_p\,\frac{\epsilon^2\bar\Omega_p^2}{\Omega_c^2}.$$

With the photon flux $\Phi_{01}$ in LP$_{01}$ related to the Rabi frequency by $\epsilon^2\bar\Omega_p^2 = \Gamma\sigma_0\Phi_{01}/A$, the photons lost per atom per flux give

$$\frac{dy}{dz} = -\frac{1}{L}\,\frac{d\,\Gamma\,\gamma_p}{\Omega_c^2}\,y + \alpha_d\,y,$$

where $\alpha_d = \alpha_{11} - \alpha_{01}$ is the differential power loss of the fibre, which makes $y$ grow. This agrees with the standard EIT result when the coherence decays at $\gamma_p/2$. Using the design value $\Omega_c^2 \simeq 4\nu_r\Delta$,

$$\boxed{{\rm OD}_{\rm eff} = \frac{d}{4}\,\frac{\Gamma}{\Delta}\,\frac{\gamma_p}{\nu_r}},\qquad y(L) = y_0\,e^{-{\rm OD}_{\rm eff} + \alpha_d L}.$$

Note that $s$ and $\bar\Omega_p$ cancel. The filter depends on $\gamma_p$ only.

**Two limits for $\gamma_p$:**

- **Cold.** $\gamma_p = G_{\rm ss}$, the steady scattering per atom. Measured with the Zeeman solver code, three-level case: 0.0077 $\nu_r$ for $b_p = 1/2$, 0.0037 $\nu_r$ for $b_p = 5/6$; 0.0045 $\nu_r$ for D1 $F'=1$ with a repumper.
- **Hot.** The saturated cycle rate is $b_p\Gamma_n/c$ with $c = 2.25$, with $1/b_p$ photons per cycle, so $\gamma_p \to \Gamma_n/2.25$, independent of $b_p$. With $\Gamma_n = \Gamma\nu_r/\Delta$,

$$\boxed{{\rm OD}_{\rm eff}^{\rm hot} \simeq \frac{d}{9(\Delta/\Gamma)^2}},$$

which is 0.44 of the bare off-resonant value $d/[4(\Delta/\Gamma)^2]$, a useful check.

| $d$ (1 cm, one atom per site) | $a$ (µm), $\nu_z$ any | Cold, $b_p = 1/2$ | Cold, $b_p = 5/6$ | Hot |
|---|---|---|---|---|
| 30 | 15 | 0.08 dB | 0.04 dB | 1.6 dB |
| 70 | 10 | 0.2 dB | 0.09 dB | 3.8 dB |
| 140 | 7 | 0.4 dB | 0.2 dB | 7.5 dB |
| 550 | 3.5 | 1.5 dB | 0.7 dB | 30 dB |

**This corrects my earlier estimate,** which guessed $\gamma \approx 0.05\,\nu_r$ for cold atoms and gave −10 dB at $d = 550$. The measured rates are about 7–14 times lower. These $d$ are two-level values; the real $\Lambda$ line strength lowers them.

**Conclusion so far: a cold column is nearly transparent to the impurity,** by the same interference that makes the cooling work: the atoms that matter most, the ones in the ground state, do not absorb it.

## Nonlinear model, with the column's response to the impurity

Let $\gamma_p$ depend on the local $y$, through the local motional state. With $A \equiv d\,\Gamma/(\Omega_c^2 L) = d/(4 L\,(\Delta/\Gamma)\nu_r)$,

$$\frac{dy}{dz} = -A\,\gamma_p(y)\,y + \alpha_d\,y.$$

**Interpolation, an assumption to be tested:** $\gamma_p(y) = \gamma_0 + (\gamma_h - \gamma_0)\,y/(y + y_c)$, with $y_c = \epsilon_c^2 = (1.3\,\eta_\perp)^2$. Then:

1. **Clean input, $y_0 \ll y_c$:** $y$ follows the exponential with ${\rm OD}_{\rm eff}$ at the cold rate.
2. **Dirty input, $y_0 \gg y_c$:** the front is hot, and $y$ falls at the hot rate, over a length $z_f/L \simeq \ln(y_0/y_c)/{\rm OD}_{\rm eff}^{\rm hot}$. For $y_0 = 10\,y_c$: 0.34 at $d = 550$, more than the column at $d = 140$.
3. **Fixed point (clamp):** if $A\gamma_0 < \alpha_d < A\gamma_h$, then $y$ tends to $y^*$ with $A\gamma_p(y^*) = \alpha_d$:

$$\frac{y^*}{y_c} = \frac{\alpha_d L - {\rm OD}_{\rm eff}}{{\rm OD}_{\rm eff}^{\rm hot} - \alpha_d L}.$$

Then $y$ would not grow beyond $y^*$, whatever the length. Example: $d = 140$, $\alpha_d L = 3$ dB: $y^* \approx 0.6\,y_c$.

**A caution about this nonlinear part.** It rests on the carrier heating the atoms, and it may not do so. Stored scans show a carrier reducing the rate, with $n_{\rm ss}$ staying small: for example, $n_{\rm ss} \approx 0.05$ at $\epsilon = 0.05$, $\delta_{\rm e} \approx 0$. If $\gamma_p$ stays near $\gamma_0$ up to $y \gg y_c$, the fixed point moves far above $y_c$ and the filter stays negligible. **So $\gamma_p(y)$ must be computed before relying on items 2 and 3.**

## Time dependence [qualitative]

Cooling begins with a hot column, where the filter is strong, so downstream atoms see cleaner light and cool first. As they cool, the filter weakens and the impurity downstream rises to its steady value. So:

- The filter helps the hot stage, when the purity limit on $s_{\max}$ matters most. This may be its most useful role.
- It does not lower the final floor.
- Input atoms always see $y_0$, so the front sets a limit, unless a front section is discarded.

## What the model does not change

- **Probe power:** absorption of LP$_{11}$ is small; the report gives $P_{\min}$ of 0.02 nW for $N = 10^4$ (case A) against 8 nW.
- **Phase:** the beat between LP$_{01}$ and LP$_{11}$ along $z$ moves between a shift of the node and a true carrier. The effects are similar through $\epsilon_{\rm eff}$, so the model uses $\lvert\epsilon\rvert$. Period about 2.5 mm at $a = 15\ \mu$m.
- **Zeeman schemes:** optimum $\delta_{\rm e}$ moves to $+0.45$–$0.75\,\nu_r$ because of light shifts. If the EIT null moves with it, the result above holds; if not, absorption is much larger (the Fano factor is $\approx 0.05$ at $\delta = 0.6\,\nu_r$). Must be computed.

## How to test it (cheap; no propagation simulation)

1. **$\gamma_p(y)$:** steady 1D master equation at $\epsilon = 0$–$0.1$, at fixed settings. About 1 min per point.
2. **Absorption coefficient directly:** from the steady state, the atomic source of the LP$_{01}$ mode, $\langle u_{01}(\hat x)\,\sigma\rangle$, with $\sigma = c\,A\rho$ in the eliminated model. It should give ${\rm OD}_{\rm eff}$ without the photon-budget argument. **This is the key check of the boxed formula.**
3. Then integrate $y(z)$ with these tables, and $P_0(z)$ from $P_0(\epsilon)$; time dependence by a quasi-steady column of slices.

## Possible uses

- A hot-stage bonus for small cores.
- Not a replacement for probe purity: it fails where purity matters most.
- A deliberate filter section without the coupling beam would have ${\rm OD} = d/[4(\Delta/\Gamma)^2]$ off resonance, but it is hard to arrange in one fibre.
