# Optimal s ramp: theory (9 Nov 2026)

Analytic reasoning, tested only against stored results (`rate_p0.dat`, `rates_vs_n.dat`, `universal_thermal*.dat`).
No ramp has been simulated. Fits are for case A, Δ = 3Γ.

## Exact statement

Populations obey Ṗ = Q(s)P, linear in P. Minimise T subject to P₀(T) ≥ P₀*, s(t) ≤ s_max. Pontryagin gives

  λ̇ = −Qᵀ(s)λ,  s*(t) = argmax_s λᵀ Q(s) P,

where λ_n is the expected remaining time from level n. So the best s at any time depends on the *whole* P(n), not only
on ⟨n⟩ or P₀. Inputs: Q(s) for given ν_r/ω_rec, b_p, Δ/Γ. Initial condition: n̄₀. Everything scales with ν_r.

## Cold end: closed form

For x = 1 − P₀: ẋ = −λ(s)[x − x_ss(s)], with λ = Λ₀s/(1 + s/s₀), x_ss = a(1 + βs). Maximising −ẋ at each x:

  x − x_ss(s) = aβ s (1 + s/s₀),  i.e.  s*(x) = s₀[√(1 + (x/a − 1)/(βs₀)) − 1].

Then ẋ = −Λ₀aβs², and

  T = (2/Λ₀)[1/s_f − 1/s_i + ln(s_i/s_f)/s₀].

**Rule:** stay above the floor by aβs(1+s/s₀), do not sit on it.

| Fit, s ≤ 1 for rates | a | β | Λ₀ (s⁻¹) | s₀ | Λ₀ in units of ν_r |
|---|---|---|---|---|---|
| b_p = 1/2, h = 0.40 | 0.0122 | 2.39 | 5125 | 0.75 | 0.114 |
| b_p = 5/6, h = 0.24 | 0.0086 | 1.66 | 5902 | 1.40 | 0.131 |

Λ₀ is close to (1−h)Γ_n/2 (0.100, 0.127), with no explicit b_p. The linear fit of x_ss is only fair: a = 0.0122 against x = 0.0100 at s = 0.01.
So targets near the floor are unreliable. The map should supply a, β, Λ₀, s₀ as functions of h and b_p.

Gain over the best constant s, from x = 0.5, in this one-variable model:

| P₀* | b_p = 1/2: ramp, constant | gain | b_p = 5/6: ramp, constant | gain |
|---|---|---|---|---|
| 0.90 | 0.87, 0.90 ms | 1.03 | 0.38, 0.40 ms | 1.03 |
| 0.95 | 1.57, 1.74 ms | 1.10 | 0.67, 0.73 ms | 1.09 |
| 0.98 | 4.62, 6.65 ms | 1.44 | 1.58, 2.08 ms | 1.32 |
| 0.985 | 10.2, 18.3 ms | 1.79 | 2.38, 3.55 ms | 1.49 |

## Adding Δ as a control

x_ss ∝ Δ⁻² and rate ∝ Δ⁻¹ (stored D1 scans: 1−P₀ ∝ Δ^−2.0, W₀ ∝ Δ^−1.0). With Δ free, the optimum is:

- **constant s*** = s₀[1 + √(1 + 8/(βs₀))]/2 (1.26 for b_p = 1/2, 2.17 for b_p = 5/6);
- **Δ rising as x falls,** so that x = 3x_ss(Δ), that is Δ ∝ x^(−1/2);
- then ẋ ∝ −x^(3/2), and T ∝ ε^(−1/2).

Gain over the best fixed Δ: 1.1 at P₀ = 0.95, 1.35 at 0.99, 1.6 at 0.998. Keeping the resonance needs Ω_c² ∝ Δ and a retuned δ_e.
A fixed Δ also scales as ε^(−1/2) up to a logarithm, **so the gain is only logarithmic.**

## Hot start: where the time actually goes

From n̄₀ = 30 the master equation needs ≈ 21 ms to P₀ = 0.9, the one-variable model ≈ 1 ms. The late rate is only ≈ 0.08 W₀.
So the time is spent draining hot atoms, not in the last phonon. Estimate:

  T_II ≈ n_c/R_sat, n_c = 2(ν_r/Γ_n)²/s;  T_III ≈ ln(n₀/n_c)/(K s), K ≈ c R_sat/18 ≈ 200 s⁻¹ (R ≈ R_sat(1 + c n/n_c), c ≈ 1.3–1.9 at s = 1, n = 50–200).

**Both fall as 1/s,** so the bottleneck is the plateau: at higher s it shrinks, and the fast channel takes over.

| s | n_c | T_III + T_II (ms) | ME, t(P₀ = 0.5) (ms) |
|---|---|---|---|
| 0.3 | 60 | 13.3 | 8.7 |
| 1 | 18 | 10.8 | 7.1 |
| 3 | 6 | 5.0 | not run |
| 5 | 3.6 | 3.3 | not run |
| 10 | 1.8 | 1.7 | not run |

The estimate is 1.5× high but has the right ratio (1.23 against 1.22). **Above s = 1 it is an extrapolation:** drift data stop at s = 1,
and an earlier worry of a friction peak (adiabatic dark-state following) is not supported: passage through the node gets more diabatic as √(n/n_c) (note, fast-channel derivation). Measure it anyway.

**A warning for naive ramps.** Lowering s raises n_c, so atoms still at n ~ 10–100 fall back onto the slow plateau. A ramp that follows
P₀ alone could strand the tail. This is why the problem needs P(n), not one variable.

## Candidate protocol [hypothesis]

1. **Hot stage:** s = s_max while the tail lasts. Possibly small Δ, as T_II ∝ Δ³/s if R_sat ∝ 1/Δ.
2. **Cold stage:** s*(x) or constant s* with rising Δ.
3. Switch: when the population above n ≈ n_c(s_final) is negligible, not when P₀ reaches a value.

**Limits on s_max:** carrier from mode impurity, s ≲ 8(Δ/3Γ)²(η⊥/ε)² (ε = 0.01: ≈ 35; ε = 0.03: ≈ 4 in case A); validity of the elimination
and band at large s·n. **A ramp does not change h,** so it does not rescue h > 1: the final low s has large n_c.

## How to get the real optimum

1. Measure v(n,s), D(n,s) from Fock states at s = 3, 5, 10 (existing `rates` study, ≈ 6 min per s). Does the friction peak?
2. ME at constant s = 3, 5 from n̄₀ = 30: tests the 1/s estimate directly.
3. Chain for P(n) with the exact kernel and tabulated R(n,s). Solve the Pontryagin problem by forward–backward sweeps (cheap, N ≈ 400).
4. Validate the resulting ramp, as piecewise-constant s, in the banded ME. One LU per stage.
5. Express as a recipe in ν_r/ω_rec, b_p, Δ/Γ, n̄₀, P₀*.

Also needed: band check at large s·n.
