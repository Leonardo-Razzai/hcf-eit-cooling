# Cooling both radial axes: theory (9 Nov 2026)

Analytic only; nothing here has been simulated. **All results so far are for one radial axis** with the other axis ignored, including its heating.

## 1. One probe cools one motional mode

Near the node a probe field is $f(x,y) \approx c_x x + c_y y$, so with $x \propto a + a^\dagger$, $y \propto b + b^\dagger$:

$$f \propto c_x a + c_y b + c_x a^\dagger + c_y b^\dagger.$$

On the red sideband (rotating-wave, $\Omega_R\sqrt n \ll \nu_r$) only the annihilation part is resonant, so **a single coherent probe couples to one motional mode,** $A = (c_x a + c_y b)/\sqrt{\lvert c_x\rvert^2+\lvert c_y\rvert^2}$. The orthogonal mode is dark in the resolved regime.

| Probe | $c_x$, $c_y$ | Mode cooled | Left uncooled |
|---|---|---|---|
| LP₁₁, lobes along x | 1, 0 | $a$ | $b$ |
| LP₁₁ rotated by $\phi$ | cos $\phi$, sin $\phi$ | $x'$ | $y'$ |
| Scalar OAM, $l = +1$, $x + iy$ | 1, i | $a_- = (a + ib)/\sqrt2$ | $a_+ = (a - ib)/\sqrt2$ |

**A single scalar OAM beam therefore does not solve it.** Angular momentum shows why: absorbing $l = +1$ and removing energy can only remove a quantum of the mode with $L = -1$. Removing an $L = +1$ quantum would need $l = -1$.
The creation part of $x + iy$ adds $L = +1$ quanta, but it is the blue sideband, suppressed by $(\Gamma/4\Delta)^2$.

## 2. Recoil couples the axes: the main change

One cycle removes one quantum from the addressed mode, but recoil heats **both** axes ($\alpha_x + \alpha_y = 0.8$ for the present dipole pattern). So per axis, with equal cycle rates:

$$v = -\frac{R}{2}(1 - 2h),\quad D = hnR \quad\Rightarrow\quad \frac{v}{D} = -\frac{1 - 2h}{2hn},$$

so **the effective value is $h_2 = 2h$** when both axes are saturated ($R$ is shared).

The other limit: x hot and saturated, y cold, where $y$ tracks $R(n_y) = R_x h/(1-h)$. This gives $h_2 = h/(1-h)$. The two estimates bracket the tail exponent $1/h_2$. They coincide at $h = 1/2$.

**Energy balance:** net removal per cycle is $1 - 2h$, so a cold steady state needs **$h < 1/2$, not $h < 1$.** (Fast channel not included, as in the 1D argument.)

**Core floor.** With rates $A_\mp \propto (n, n+1)$, equal axes:

  $n_{\rm core} = \xi$ (1 + 2h)/(1 − 2h), $\xi = (\Gamma/4\Delta)^2$,  against $\xi$(1 + h)/(1 − h) in 1D.

Table at $\nu_z = 300$ kHz:

| a (µm) | $b_p$ | h | 2h | $h/(1-h)$ | 2D/1D floor | $P_0$ both axes, $3\Gamma$ | $6\Gamma$ |
|---|---|---|---|---|---|---|---|
| 3.5 | 1/2 | 0.09 | 0.19 | 0.10 | 1.2 | 0.98 | 0.995 |
| 3.5 | 5/6 | 0.06 | 0.11 | 0.06 | 1.1 | 0.98 | 0.996 |
| 7 | 1/2 | 0.19 | 0.38 | 0.23 | 1.5 | 0.97 | 0.99 |
| 7 | 5/6 | 0.11 | 0.23 | 0.13 | 1.3 | 0.98 | 0.99 |
| 10 | 1/2 | 0.27 | 0.54 | 0.37 | 1.9 | 0.96 | 0.99 |
| 10 | 5/6 | 0.16 | 0.32 | 0.19 | 1.4 | 0.97 | 0.99 |
| 15 | 1/2 | 0.40 | 0.81 | 0.68 | 4.0 | 0.88 | 0.97 |
| 15 | 5/6 | 0.24 | 0.48 | 0.32 | 1.7 | 0.96 | 0.99 |
| 21 | 1/2 | 0.56 | 1.1 | 1.3 | no cold state | | |
| 21 | 5/6 | 0.34 | 0.68 | 0.51 | 2.6 | 0.93 | 0.98 |

Weak probe, core only, so upper bounds on $P_0$: no s dependence, no tail. **Case A with $b_p = 1/2$ (the current paper) moves into the window where the mean is set by the cut-off,** and the floor rises 4×.
With $2h$: $h_2 < 1/3$ needs $a < 10$ µm for $b_p = 5/6$ and 6 µm for $b_p = 1/2$.

**Consequences at $\nu_z = 300$ kHz:** the argument for small a and large $b_p$ gets stronger. See 2b for other depths.

**Polarisation matters.** $\alpha$ depends on the emission pattern: dipole along $z$ gives $\alpha_x = \alpha_y = 0.4$, radial sum 0.8; along x gives 0.2 and 0.4, sum 0.6; averaged, 2/3. Guided light is transverse, so the sum may be below 0.8. Must come from the Zeeman model.

## 2b. The limits are on $a/\nu_z$, not on $a$

$\nu_r = 0.357\,\nu_z/a(\mu{\rm m})$ at fixed lattice wavelength, so

$$h = 0.403\left(\frac{a}{15\ \mu{\rm m}}\right)\left(\frac{300\ {\rm kHz}}{\nu_z}\right)\left(\frac{1/2}{b_p}\right).$$

 The table above is for $\nu_z = 300$ kHz only. In fibres the depth can be
raised with lattice power, so the limits on a scale with $\nu_z$: $2h < 1/3$ needs $a < 10\ \mu$m $\times(\nu_z/300$ kHz$)$ for $b_p = 5/6$, and $6\ \mu$m for $b_p = 1/2$; $h < 1/2$ needs $a < 31\ \mu$m and $19\ \mu$m, times the same factor.

**Unchanged:** the map in $\nu_r/\omega_{\rm rec}$, one probe per mode, the doubling of h, the condition $h < 1/2$. Only the mapping to a changes.

**Costs of a deeper lattice (estimates).** Beam power $\propto \nu_z^2 a^2$ (≈ 0.14 W per beam at $a = 15$ µm, 300 kHz, from $\alpha \approx 690$ a.u.). Lattice scattering $\propto$ $\nu_z^2$ (4 s⁻¹ at 300 kHz); its radial
heating relative to W is unchanged, so it is harmless during cooling, but it limits storage. Purity gets harder: $\epsilon_c^2 \propto \eta_\perp^2 \propto 1/(a\nu_z)$, −31 dB at $a = 15$ µm, 300 kHz, −36 dB at 1 MHz.
The excited-state light shift grows with $U_0$: at 6 mK, $U_0/h \approx 120$ MHz $\approx 20\Gamma$, and it is **not in the model.** Not relevant to the map, but it limits how far this goes. The elimination needs $\nu_r \ll \Delta$, easily met.

**Smaller a helps both h and purity, and OD. A deeper lattice helps h but hurts purity.** So $\nu_z$ is a real design knob, not a free gain.

**For a memory: cool deep, store shallow.** Adiabatic lowering conserves n, so cooling at high $\nu_z$ (low h) and then lowering (low scattering, low $\delta_{\rm diff}$) keeps both. To be simulated.

## 3. Ways to cool both axes

| Scheme | Field | Works when | Rate per axis | Comment |
|---|---|---|---|---|
| **A. Alternating lobes** | LP₁₁a, then LP₁₁b | any switching | ≈ W/2 | simplest; needs mode switching ≲ 1/W ~ 0.1–1 ms |
| **B. Rotating pattern** | LP₁₁ at angle $\omega_{\rm rot}$ t | $W \ll \omega_{\rm rot} \ll \Gamma_n$ | ≈ W/2 | two frequencies $\delta f = 2\omega_{\rm rot}$ within the narrow line; continuous |
| **C. Trap anisotropy** | one LP₁₁ at 45° | $W \ll \lvert\nu_x - \nu_y\rvert \lesssim \Gamma_n$ | ≈ W/2 | beat between axes rotates the dark mode into the bright one; cw, single beam |
| **D. Two OAM** | $l = +1$ and $l = -1$, incoherent | no interference between them | ≈ W/2 | same as B if the frequencies differ |
| **E. Vector vortex** | TE$_{01}$ or TM₀₁ | the two circular polarisations drive independent Raman paths | possibly ≈ W | **most interesting** |

**C, numbers.** $W \approx 0.05$–0.08 $\nu_r$, $\Gamma_n = 0.33$ $\nu_r$ at $3\Gamma$: needs anisotropy of ≈ 10–30 %. At $6\Gamma$, $W \approx 0.03$, $\Gamma_n = 0.17$: ≈ 5–15 %. A fixed lobe aligned with a trap axis cools only that axis, so the angle matters, including any rotation along z.

**E, why.** TE$_{01}$ has $E \propto J_1(ur/a)(-\sin\phi, \cos\phi)$. In the circular basis, $\sigma^+$ carries $e^{+i\phi} \propto x + iy$ and $\sigma^-$ carries $e^{-i\phi} \propto x - iy$. Total angular momentum J = 0. So **$\sigma^+$ cools $a_-$ and $\sigma^-$ cools $a_+$: one beam addresses both modes,** provided the two polarisations drive distinct Raman paths with no interference. It is a true fibre eigenmode, with an exact axial zero, unlike LP₁₁, which is a superposition. **This needs the Zeeman model,** which is needed anyway. A magnetic field or atomic structure that mixes the paths would reduce it to one mode.

**Hot stage [hypothesis].** The one-mode rule is for resolved sidebands. For $\Omega_R\sqrt n \gg \nu_r$ the hill $\delta_B \propto 1 + \tan^2\theta(r)$ of a donut depends on r, so the Sisyphus channel should act on both axes, even with a scalar OAM beam, while the cold stage does not.
A possible protocol: donut or OAM when hot, scheme A, C or E when cold.

## 4. Numerical results (9 Nov 2026)

`simulation/studies/two_axes.py`, results in `simulation/results/two_axes_*.dat`.

**Rate chain on $P(n_x, n_y)$,** no parameter beyond those fixed in 1D. In 1D it gives $1-P_0$ of the master-equation map within 10 % for $h \lesssim 0.4$ (median 1.02); it overestimates $\langle n\rangle$ by 1.2–1.6× at $h \approx 0.4$ and fails above $h \approx 0.6$. So $P_0$ is reliable and $\langle n\rangle$ is an upper bound.

Steady states, both axes, $s = 1$:

| Point | $h$ | $P_0$, 1D | $P_{00}$, 2D, $3\Gamma$ | $\langle n_x\rangle$, 2D, $3\Gamma$ | $P_{00}$, 2D, $6\Gamma$ |
|---|---|---|---|---|---|
| Case A, $b_p = 1/2$ | 0.40 | 0.950 | 0.755 | 2.2 | 0.941 |
| Case A, $b_p = 5/6$ | 0.24 | 0.975 | 0.927 | 0.10 | 0.982 |
| $\nu_r/\omega_{\rm rec} = 3$, $b_p = 5/6$ | 0.16 | 0.979 | 0.946 | 0.045 | 0.987 |
| $\nu_r/\omega_{\rm rec} = 5$, $b_p = 5/6$ | 0.10 | 0.981 | 0.956 | 0.029 | 0.989 |
| Case B, $b_p = 5/6$ | 0.06 | 0.982 | 0.961 | 0.023 | 0.990 |

The floor ratio follows the formula of Section 2 (1.77 against 1.75 at $h = 0.24$, 1.41 at $h = 0.16$). Case A with $b_p = 1/2$ lies at $2h \approx 0.8$, near where the chain fails, so it is indicative. Above $h \approx 0.5$, $P_{00}$ falls strongly, but those runs have not converged.

**Exact check,** eliminated master equation at $N = 6$ per axis, against the chain at the same $N$ ($3\Gamma$, $s = 1$). One-axis limit: 0.97207 against 0.97211.

| Point | $h$ | $P_{00}$, exact | $P_{00}$, chain |
|---|---|---|---|
| $\nu_r/\omega_{\rm rec} = 8.43$, $b_p = 1/2$ | 0.09 | 0.933 | 0.938 |
| $\nu_r/\omega_{\rm rec} = 5$, $b_p = 1/2$ | 0.16 | 0.919 | 0.924 |
| Case A, $b_p = 5/6$ | 0.24 | 0.927 | 0.931 |
| $\nu_r/\omega_{\rm rec} = 3$, $b_p = 1/2$ | 0.26 | 0.893 | 0.898 |

Larger $h$ cannot be tested at this $N$.

**Switching must be slow.** At $h = 0.09$: $P_{0x} = 0.777$, 0.896, 0.955, 0.964, 0.964 for stage length $\tau\nu_r = 4$, 8, 16, 32, 64, so $\tau \gtrsim 10/\Gamma_n$ ($\approx 0.7$ ms in case A). Averaging the two Liouvillians is not alternation: it adds the probes coherently in $H$ and does not cool. One point; mechanism not analysed. This limits scheme A from below, and suggests checking schemes B and C.

**Time to 0.9 of the final $P_{00}$** (chain), both axes, $3\Gamma$, $s = 1$: case A, $b_p = 5/6$: 10, 22, 32 ms from $\bar n_0 = 3$, 10, 30; case B, $b_p = 1/2$: 2.5, 6.0, 9.3 ms.

## 5. Dimensionless control parameter

$$\zeta \equiv k_L a\,\frac{\omega_{\rm rec}}{\nu_z},\qquad \frac{\nu_r}{\omega_{\rm rec}} = \frac{2.11}{\zeta},\qquad h = 0.188\,\frac{\zeta}{b_p},\qquad \eta_\perp = \frac{G_{11}}{ka}\sqrt{\frac{\zeta}{2.11}}.$$

A cold state needs $\zeta < 2.66\,b_p$, a clean tail $\zeta < 0.89\,b_p$. Case A at 300 kHz has $\zeta = 1.07$, case B 0.25. At fixed $\zeta$, a larger core tightens mode purity as $1/(ka)$.

## 6. Next numerics

1. Schemes B, C and E; scheme E needs the Zeeman model.
2. Rerun the map in 2D if the configuration changes.
