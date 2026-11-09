# Zeeman structure of the cooling scheme: results (9 Nov 2026)

Eliminated, banded master equation on all 8 ground sublevels of $^{87}$Rb, one radial axis, D1, $\Delta = 3\Gamma$. Code `simulation/eitsim/zeeman.py`, studies `simulation/studies/zeeman_scan.py`, data `simulation/results/zeeman_*.dat`. Checks: reduces to the three-level solver, decay sum rule and branching ratios hold. The Zeeman $\Lambda$ excludes the far $F'=1$ level, and hyperfine schemes exclude the other excited level, here 142 $\Gamma$ away.

## Dark states are unavoidable

A fixed polarisation on $F \to F'$ leaves dark combinations whenever $2F+1$ exceeds the rank of the coupling, at every motional state:

| Configuration | Probe leg, dark of total | Coupling leg, dark of total |
|---|---|---|
| D1 $F'=1$, probe on $F=2$ | 2 of 5 | 1 of 3 |
| D1 $F'=2$, probe on $F=1$ | 0 of 3 | 1 of 5 |
| D1 $F'=2$, probe on $F=2$ | 1 of 5 | 0 of 3 |

Unmixed, there is no unique steady state (182 of the 1008 scan points, those with $L_z$ along a pure polarisation). So the closed $\Lambda$ results for $b_p = 5/6$ cannot be reached as they stand.

## Results

$s \le 3$, best over the scanned settings, $N = 12$, $B = 2$; rates in units of $\nu_r$. Case A trap, except the last row. Convergence checked up to $N = 28$, $B = 3$: $P_0$ and $W_0$ within 0.8 %; $\langle n\rangle$ not converged (up 25–45 %, same tail problem as before).

| Scheme | $P_0 \geq 0.95$ | $P_0 \geq 0.97$ | $P_0 \geq 0.99$ |
|---|---|---|---|
| Closed $\Lambda$, $b_p = 5/6$, no Zeeman (reference) | 0.077 (at $P_0 = 0.977$) | 0.077 | see `lines_delta_A.dat` |
| D1 $F'=1$, lin $\perp$ lin, repumper | 0.027 | 0.0225 | 0.0115 |
| D1 $F'=1$, lin $\perp$ lin, transverse field $L_x = 0.1\nu_r$ | 0.015 | 0.015 | 0.0085 |
| Zeeman $\Lambda$ in $F=2$, $\sigma^-$ probe, $\sigma^+$ coupling, $F=1$ repumper | 0.023 | 0.022 | 0.012 |
| Zeeman $\Lambda$, $\nu_r/\omega_{\rm rec} = 5$ | 0.039 | 0.039 | 0.025 |

The best hyperfine repumper rates are $\gamma \approx 0.03$–$0.1\,\nu_r$, with $\delta_{\rm e}$ moved to $+0.45$ to $+0.75\,\nu_r$ by light shifts. $F'=2$ configurations are worse: best $W_0 = 0.0096$ at $P_0 = 0.925$ (probe on $F=1$).

## Reading

- **The penalty is about 3× in rate relative to the closed $\Lambda$,** not the 15× of the first unoptimised point. In case A, $W_0 \approx 1.0\times10^3$ s$^{-1}$ at $P_0 = 0.97$, against $3.4\times10^3$ s$^{-1}$.
- **$F'=1$ remains the best hyperfine choice,** contrary to what the dark-state count alone suggests.
- **Hyperfine and Zeeman $\Lambda$ are comparable on one axis.** The Zeeman $\Lambda$ needs a field that resolves the Zeeman structure, and repumping at $\gamma \approx \nu_r$; it should be compared on recoil with both axes. Estimated photons per cycle: 4, so $h \approx 0.68$ in case A at 300 kHz by counting, which would rule it out with two axes there; that needs a 2D test.
- **Field remixing needs very low fields:** $L_x = 0.1\nu_r$ is about 1 mG in case A. A repumper looks more practical.

## Not included yet

Both radial axes; the other excited level as a neighbour; real polarisation structure of the modes, including TE$_{01}$/TM$_{01}$; quadratic Zeeman shifts; time evolution from thermal states; other $s$, $\Delta$ settings than those scanned. The grid is coarse, and the best points can improve.
