# Optimal s ramp: theory (9 Nov 2026)

## Results of steps 1–3 (9 Nov 2026; `studies/ramp.py`, `convergence.py rates_hi`, `universal.py hotstage`)

Case A, $b_p = 1/2$, $\Delta = 3\Gamma$, one axis.

**Step 1: drift from Fock states at high s.** The removal rate keeps growing with s at moderate n (n0 = 30: 14, 20, 36×10³ s⁻¹ at s = 3, 5, 10) and saturates
or falls at n ≳ 100. **No friction peak is seen.** At n = 1 the drift changes sign for s ≥ 5 (+104, +437 s⁻¹): the mean rises while $P_0$ grows, because recoil builds
the tail. **Caveat:** at s = 5 the high-n drift is not converged in the coherence band (band 4 against 6 differs by up to 2× for n ≥ 30), while capture times agree within 1 %.

**Step 2: constant s from $\bar n_0 = 30$ (master equation).** t($P_0$ = 0.5) = 8.7, 7.1, 4.1, 3.0 ms at s = 0.3, 1, 3, 5. **Faster with s, but about as $s^{-0.5}$, not 1/s:**
the estimate of the table below is 1.2–1.5× high and its s scaling too strong. Higher s also lowers the final $P_0$ (0.89 at s = 3, 0.84 at s = 5).

**Step 3: chain and optimal control.** Population chain $\dot P = Q(s)P$: the GRM Franck–Condon recoil kernel, each Fock column scaled so that its drift equals the master-equation drift
for n ≥ 5, and the calibrated GRM saturation below (where the drift mixes removal and tail building). Checked against the master equation from $\bar n_0 = 30$ at constant s:

| s | t($P_0$ = 0.5, 0.7, 0.8), ME (ms) | chain (ms) | final $P_0$, ME / chain |
|---|---|---|---|
| 0.3 | 8.7, 16.5, 23.7 | 8.3, 13.5, 17.8 | 0.974 / 0.978 |
| 1 | 7.0, 11.3, 14.5 | 7.0, 10.8, 13.7 | 0.954 / 0.954 |
| 3 | 4.1, 6.1, 7.9 | 3.7, 5.4, 7.0 | 0.891 / 0.897 |
| 5 | 3.0, 4.6, 6.4 | 2.8, 4.3, 5.9 | 0.844 / 0.852 |

The chain is within 5 % at s = 1 and 10–25 % fast elsewhere. Minimum-time s(t) by forward–backward (Pontryagin) sweeps over 30 stages and 13 values of s:

| $\bar n_0$ | $P_0^*$ | best constant s: time | ramp, $s \le 10$ | $s \le 5$ | $s \le 3$ | two stages ($s_{\max}$, then s = 1), $s_{\max}$ = 10, 5, 3 |
|---|---|---|---|---|---|---|
| 30 | 0.5 | s = 10: 2.4 ms | 2.3 ms | 2.8 | 3.7 | — |
| 30 | 0.8 | s = 5: 6.0 ms | 4.3 ms (1.4×) | 5.1 | 6.6 | 4.3, 5.1, 6.6 ms |
| 30 | 0.9 | s = 2: 13.5 ms | 6.3 ms (2.1×) | 7.2 (1.9×) | 9.1 (1.5×) | 6.4, 7.2, 9.3 ms |
| 10 | 0.9 | s = 2: 10.5 ms | 5.1 ms (2.0×) | — | — | — |

**Findings.**
1. **The optimum is bang-bang:** $s_{\max}$ for the first 60–75 % of the time (hot stage), a short intermediate step, then s ≈ 1 → 0.7 (cold stage). This is the candidate protocol below.
2. **A two-stage protocol is within 2 % of the optimum.** Switch at 3.9, 5.2, 7.1 ms for $s_{\max}$ = 10, 5, 3 ($P_0^*$ = 0.9). A full control solution is not needed.
3. **The gain is large only for targets near or below the high-s floor:** 2× at $P_0^*$ = 0.9, 1.4× at 0.8, none at 0.5. It comes from using a high s that alone could not reach the target.
4. **It depends on $s_{\max}$**, and $s_{\max}$ is set by mode purity ($s \lesssim 4$ for $\epsilon = 0.03$) and by the validity of the chain at high s.

**Limits.** One axis; $b_p = 1/2$ only; chain calibrated on drift data that are not band-converged at high n and s; the chain is 10–25 % optimistic at constant s.
**Step 4** (validate the two-stage protocol in the banded master equation, one LU per stage) is cheap and decides whether the gain is real.


Analytic reasoning, tested only against stored results (`rate_p0.dat`, `rates_vs_n.dat`, `universal_thermal*.dat`).
No ramp has been simulated. Fits are for case A, $\Delta = 3\Gamma$.

## Exact statement

Populations obey $\dot P = Q(s)P$, linear in $P$. Minimise $T$ subject to $P_0(T) \geq P_0^*$, $s(t) \leq s_{\rm max}$. Pontryagin gives

$$\dot\lambda = -Q^T(s)\,\lambda,\qquad s^*(t) = \arg\max_s\, \lambda^T Q(s) P,$$

where $\lambda_n$ is the expected remaining time from level n. So the best s at any time depends on the *whole* P(n), not only
on $\langle n\rangle$ or $P_0$. Inputs: Q(s) for given $\nu_r/\omega_{\rm rec}$, $b_p$, $\Delta/\Gamma$. Initial condition: $\bar n_0$. Everything scales with $\nu_r$.

## Cold end: closed form

For $x = 1 - P_0$: $\dot x = -\lambda(s)[x - x_{\rm ss}(s)]$, with $\lambda = \Lambda_0 s/(1 + s/s_0)$, $x_{\rm ss} = a(1 + \beta s)$. Maximising −ẋ at each x:

$$x - x_{\rm ss}(s) = a\beta s\,(1 + s/s_0),\quad\text{i.e.}\quad s^*(x) = s_0\left[\sqrt{1 + \frac{x/a - 1}{\beta s_0}} - 1\right].$$

Then $\dot x = -\Lambda_0 a\beta s^2$, and

$$T = \frac{2}{\Lambda_0}\left[\frac{1}{s_f} - \frac{1}{s_i} + \frac{\ln(s_i/s_f)}{s_0}\right].$$

**Rule:** stay above the floor by $a\beta s(1+s/s_0)$, do not sit on it.

| Fit, $s \leq 1$ for rates | a | $\beta$ | $\Lambda_0$ (s⁻¹) | $s_0$ | $\Lambda_0$ in units of $\nu_r$ |
|---|---|---|---|---|---|
| $b_p = 1/2$, $h = 0.40$ | 0.0122 | 2.39 | 5125 | 0.75 | 0.114 |
| $b_p = 5/6$, $h = 0.24$ | 0.0086 | 1.66 | 5902 | 1.40 | 0.131 |

$\Lambda_0$ is close to $(1-h)\Gamma_n/2$ (0.100, 0.127), with no explicit $b_p$. The linear fit of $x_{\rm ss}$ is only fair: $a = 0.0122$ against $x = 0.0100$ at $s = 0.01$.
So targets near the floor are unreliable. The map should supply a, $\beta$, $\Lambda_0$, $s_0$ as functions of h and $b_p$.

Gain over the best constant s, from $x = 0.5$, in this one-variable model:

| $P_0^*$ | $b_p = 1/2$: ramp, constant | gain | $b_p = 5/6$: ramp, constant | gain |
|---|---|---|---|---|
| 0.90 | 0.87, 0.90 ms | 1.03 | 0.38, 0.40 ms | 1.03 |
| 0.95 | 1.57, 1.74 ms | 1.10 | 0.67, 0.73 ms | 1.09 |
| 0.98 | 4.62, 6.65 ms | 1.44 | 1.58, 2.08 ms | 1.32 |
| 0.985 | 10.2, 18.3 ms | 1.79 | 2.38, 3.55 ms | 1.49 |

## Adding $\Delta$ as a control

$x_{\rm ss} \propto \Delta^{-2}$ and rate $\propto \Delta^{-1}$ (stored D1 scans: $1-P_0 \propto \Delta^{-2.0}$, $W_0 \propto \Delta^{-1.0}$). With $\Delta$ free, the optimum is:

- **constant** $s^* = s_0\left[1 + \sqrt{1 + 8/(\beta s_0)}\right]/2$ (1.26 for $b_p = 1/2$, 2.17 for $b_p = 5/6$);
- **$\Delta$ rising as x falls,** so that $x = 3x_{\rm ss}(\Delta)$, that is $\Delta \propto x^{-1/2}$;
- then $\dot x \propto -x^{3/2}$, and $T \propto \epsilon^{-1/2}$.

Gain over the best fixed $\Delta$: 1.1 at $P_{0} = 0.95$, 1.35 at 0.99, 1.6 at 0.998. Keeping the resonance needs $\Omega_c^2 \propto \Delta$ and a retuned $\delta_{\rm e}$.
A fixed $\Delta$ also scales as $\epsilon^{-1/2}$ up to a logarithm, **so the gain is only logarithmic.**

## Hot start: where the time actually goes

From $\bar n_0 = 30$ the master equation needs ≈ 21 ms to $P_{0} = 0.9$, the one-variable model ≈ 1 ms. The late rate is only ≈ 0.08 $W_0$.
So the time is spent draining hot atoms, not in the last phonon. Estimate:

$$T_{\rm II} \approx \frac{n_c}{R_{\rm sat}},\quad n_c = \frac{2(\nu_r/\Gamma_n)^2}{s};\qquad T_{\rm III} \approx \frac{\ln(n_0/n_c)}{Ks},\quad K \approx \frac{cR_{\rm sat}}{18} \approx 200\ {\rm s^{-1}},$$

from $R \approx R_{\rm sat}(1 + cn/n_c)$, with $c \approx 1.3$–$1.9$ at $s = 1$, $n = 50$–$200$.

**Both fall as 1/s,** so the bottleneck is the plateau: at higher s it shrinks, and the fast channel takes over.

| s | $n_c$ | $T_{\rm III}$ + $T_{\rm II}$ (ms) | ME, t($P_{0} = 0.5$) (ms) |
|---|---|---|---|
| 0.3 | 60 | 13.3 | 8.7 |
| 1 | 18 | 10.8 | 7.1 |
| 3 | 6 | 5.0 | not run |
| 5 | 3.6 | 3.3 | not run |
| 10 | 1.8 | 1.7 | not run |

The estimate is 1.5× high but has the right ratio (1.23 against 1.22). **Above $s = 1$ it is an extrapolation:** drift data stop at $s = 1$,
and an earlier worry of a friction peak (adiabatic dark-state following) is not supported: passage through the node gets more diabatic as $\sqrt{n/n_c}$ (note, fast-channel derivation). Measure it anyway.

**A warning for naive ramps.** Lowering s raises $n_c$, so atoms still at n ~ 10–100 fall back onto the slow plateau. A ramp that follows
$P_0$ alone could strand the tail. This is why the problem needs P(n), not one variable.

## Candidate protocol [hypothesis]

1. **Hot stage:** s = $s_{\rm max}$ while the tail lasts. Possibly small $\Delta$, as $T_{\rm II} \propto \Delta^3/s$ if $R_{\rm sat} \propto 1/\Delta$.
2. **Cold stage:** s*(x) or constant s* with rising $\Delta$.
3. Switch: when the population above n ≈ $n_c$($s_{\rm final}$) is negligible, not when $P_0$ reaches a value.

**Limits on $s_{\rm max}$:** carrier from mode impurity, $s \lesssim 8(\Delta/3\Gamma)^2(\eta_\perp/\epsilon)^2$ ($\epsilon = 0.01$: ≈ 35; $\epsilon = 0.03$: ≈ 4 in case A); validity of the elimination
and band at large $sn$. **A ramp does not change h,** so it does not rescue $h > 1$: the final low s has large $n_c$.

## How to get the real optimum

1. Measure $v(n,s)$, $D(n,s)$ from Fock states at $s = 3$, 5, 10 (existing `rates` study, ≈ 6 min per s). Does the friction peak?
2. ME at constant $s = 3$, 5 from $\bar n_0 = 30$: tests the 1/s estimate directly.
3. Chain for P(n) with the exact kernel and tabulated $R(n,s)$. Solve the Pontryagin problem by forward–backward sweeps (cheap, $N \approx 400$).
4. Validate the resulting ramp, as piecewise-constant s, in the banded ME. One LU per stage.
5. Express as a recipe in $\nu_r/\omega_{\rm rec}$, $b_p$, $\Delta/\Gamma$, $\bar n_0$, $P_0^*$.

Also needed: band check at large $sn$.
