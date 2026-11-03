"""Numbers and plot data for the background report (pure Python, no numpy).

Self-contained scaling estimates (not using eitsim): capillary LP_lm modes,
87Rb D1 Lambda system, core-size / fibre-length trade-off, NANF loss budget.
Writes pgfplots tables to simulation/results/report/ (read by
paper/background_report/report.tex) and prints the tables used in the text.
"""
import math, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "report")

HBAR = 1.054571817e-34
M_RB = 1.44316060e-25          # 87Rb mass [kg]
GAMMA = 2 * math.pi * 5.746e6  # D1 natural linewidth [rad/s]
LAM_D1 = 794.98e-9
K_D1 = 2 * math.pi / LAM_D1
LAM_L = 1064e-9
K_L = 2 * math.pi / LAM_L
ISAT = 4.49e-3 / 1e-4          # D1 pi-polarised saturation intensity ~4.49 mW/cm^2 -> W/m^2
ALPHA = 0.4                    # recoil projection onto one radial axis (pi dipole)

os.makedirs(OUT, exist_ok=True)


def J(n, x, N=4000):
    """Bessel J_n via Bessel's integral (Simpson); accurate to ~1e-10 for |x| < 60."""
    h = math.pi / N
    s = 0.0
    for i in range(N + 1):
        t = i * h
        w = 1 if i in (0, N) else (4 if i % 2 else 2)
        s += w * math.cos(n * t - x * math.sin(t))
    return s * h / 3 / math.pi


def zeros(n, count):
    out, x, step = [], 0.5, 0.05
    f0 = J(n, x)
    while len(out) < count:
        x1 = x + step
        f1 = J(n, x1)
        if f0 * f1 < 0:
            lo, hi = x, x1
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                if J(n, lo) * J(n, mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            out.append(0.5 * (lo + hi))
        x, f0 = x1, f1
    return out


u01 = zeros(0, 1)[0]
u1 = zeros(1, 10)
J1u01 = J(1, u01)


def G(u):
    """On-axis gradient of power-normalised LP_1m, in units of LP01 peak amplitude per core radius a."""
    return (u / math.sqrt(2)) * abs(J1u01) / abs(J(2, u))


# Gaussian fit of LP01 (maximise overlap), w/a
best = (0, 0)
for i in range(400, 900):
    wa = i / 1000
    num = den1 = den2 = 0.0
    N = 400
    for k in range(1, N + 1):
        r = k / N
        f = J(0, u01 * r, 200)
        g = math.exp(-r * r / wa / wa)
        num += f * g * r; den1 += f * f * r; den2 += g * g * r
    ov = num * num / den1 / den2
    if ov > best[0]:
        best = (ov, wa)
W_OVER_A = best[1]
print(f"u01={u01:.4f}  Gaussian fit w/a={W_OVER_A:.3f} (overlap {best[0]:.4f})")

print("\nm   u_1m     G_1m    G/(0.46u^1.5)  loss ratio (u/u01)^2   n_lin(a=15um,nu=5kHz)")
rows = []
for m, u in enumerate(u1, 1):
    g = G(u)
    rows.append((m, u, g))
print_rows = []
nu_ref, a_ref = 2 * math.pi * 5e3, 15e-6
x0_ref = math.sqrt(HBAR / (2 * M_RB * nu_ref))
for m, u, g in rows:
    nlin = 0.4 * (a_ref / (u * x0_ref)) ** 2
    print(f"{m:2d} {u:8.4f} {g:8.3f} {g/(0.46*u**1.5):8.3f} {(u/u01)**2:10.1f} {nlin:10.1f}")
with open(os.path.join(OUT, "Gfactor.dat"), "w") as f:
    f.write("m u G loss nlin\n")
    for m, u, g in rows:
        f.write(f"{m} {u:.4f} {g:.4f} {(u/u01)**2:.3f} {0.4*(a_ref/(u*x0_ref))**2:.3f}\n")

# Transverse profiles along x at equal power, in units of LP01 peak amplitude
with open(os.path.join(OUT, "profiles.dat"), "w") as f:
    f.write("x LP01 LP11 LP13 LP15\n")
    for i in range(-200, 201):
        x = i / 200
        r = abs(x)
        vals = [J(0, u01 * r, 400)]
        for m in (1, 3, 5):
            u = u1[m - 1]
            amp = math.sqrt(2) * abs(J1u01) / abs(J(2, u))
            vals.append(amp * J(1, u * r, 400) * (1 if x >= 0 else -1))
        f.write(f"{x:.4f} " + " ".join(f"{v:.5f}" for v in vals) + "\n")


def R(delta, Delta, Oc, Op=1.0, gam=GAMMA):
    """Weak-probe absorption rate of the Lambda system; delta = two-photon detuning, probe detuning = Delta + delta."""
    dp = Delta + delta
    if delta == 0:
        return 0.0
    return Op**2 * gam * delta**2 / (gam**2 * delta**2 + 4 * (dp * delta - Oc**2 / 4) ** 2)


# Fano spectrum for nu = 2pi x 5 kHz, Delta = 3 Gamma
nu = 2 * math.pi * 5e3
Delta = 3 * GAMMA
Oc = math.sqrt(4 * Delta * nu)
Gn = GAMMA * Oc**2 / (4 * Delta**2)
print(f"\nFano: Delta=3Gamma, Oc=2pi x {Oc/2/math.pi/1e3:.1f} kHz, narrow width=2pi x {Gn/2/math.pi:.0f} Hz, "
      f"A+/A- = {R(-nu,Delta,Oc)/R(nu,Delta,Oc):.2e}")
with open(os.path.join(OUT, "fano.dat"), "w") as f:
    f.write("d R\n")
    peak = R(Oc**2 / (4 * Delta), Delta, Oc)
    for i in range(-1500, 1501):
        d = 2 * math.pi * (i * 10.0)  # -15..15 kHz
        f.write(f"{d/2/math.pi/1e3:.3f} {R(d, Delta, Oc)/peak:.6f}\n")


def scenario(a, nu_z=2 * math.pi * 300e3, Delta_over_G=3.0, m=1):
    w = W_OVER_A * a
    nu_r = nu_z * math.sqrt(2) / (K_L * w)
    x0 = math.sqrt(HBAR / (2 * M_RB * nu_r))
    g = G(u1[m - 1])
    eta = g * x0 / a
    eta_sp = K_D1 * x0
    Dl = Delta_over_G * GAMMA
    Oc = math.sqrt(4 * Dl * nu_r)
    Gn = GAMMA * Oc**2 / (4 * Dl**2)
    W_std = eta**2 * Oc**2 / 9 / GAMMA            # Omega_p = Omega_c / 3, carrier present
    W_std = min(W_std, Gn / 4)
    W_node = Gn / 4                                # saturation parameter s = 1
    Oeff = math.sqrt(GAMMA * Gn / 2)
    O01 = Oeff / (g * x0 / a)                      # LP01-equivalent peak Rabi frequency
    P = 0.847 * a**2 * 2 * ISAT * (O01 / GAMMA) ** 2  # power in LP_1m probe
    nss_fano = (GAMMA / (4 * Dl)) ** 2
    return dict(nu_r=nu_r, x0=x0, eta=eta, eta_sp=eta_sp, rec=ALPHA * eta_sp**2,
                Gn=Gn, W_std=W_std, W_node=W_node, P=P, nss=nss_fano, Oc=Oc)


print("\nScenario table (nu_z = 300 kHz, 1064 nm lattice, Delta = 3 Gamma, LP11 probe)")
for a in (15e-6, 3.5e-6):
    s = scenario(a)
    print(f"a={a*1e6:.1f}um nu_r=2pi x {s['nu_r']/2/math.pi/1e3:.2f} kHz x0={s['x0']*1e9:.0f} nm "
          f"eta_eff={s['eta']:.4f} eta_sp={s['eta_sp']:.2f} recoil/photon={s['rec']:.2f} "
          f"Gamma_n=2pi x {s['Gn']/2/math.pi:.0f} Hz  W_std={s['W_std']:.1f}/s  W_node={s['W_node']:.0f}/s "
          f"P_probe={s['P']*1e9:.2f} nW  Oc=2pi x {s['Oc']/2/math.pi/1e3:.0f} kHz")

# Trend: cooling-rate ceilings versus core radius for m = 1, 3, 5
with open(os.path.join(OUT, "trend_a.dat"), "w") as f:
    f.write("a nur Wnode Wstd1 Wstd3 Wstd5\n")
    for i in range(30):
        a = (3 + i * 1.0) * 1e-6
        s1, s3, s5 = scenario(a, m=1), scenario(a, m=3), scenario(a, m=5)
        f.write(f"{a*1e6:.1f} {s1['nu_r']/2/math.pi/1e3:.3f} {s1['W_node']:.3f} "
                f"{s1['W_std']:.4f} {s3['W_std']:.4f} {s5['W_std']:.4f}\n")

# Raman sideband cooling (node geometry), continuous: Gamma_op = nu/3, Omega_sb = Gamma_op/2
print("\nRSC (node geometry, continuous, Gamma_op = nu_r/3):")
for a in (15e-6, 3.5e-6):
    s = scenario(a)
    gop = s["nu_r"] / 3
    osb = gop / 2
    W = osb**2 / gop / (1 + 2 * osb**2 / gop**2)
    print(f"a={a*1e6:.1f}um  W_RSC ~ {W:.0f}/s   (sideband Rabi 2pi x {osb/2/math.pi:.0f} Hz)")


# ---------------------------------------------------------------------------
# Core-size trade-off: HOM loss over a fibre length L versus the gain in eta.
# Power loss alpha_lm = alpha01(a) (u_lm/u01)^2, alpha01(a) = A15 (15 um / a)^3 (Marcatili scaling).
# A residual LP01 amplitude eps0 at the input grows relative to LP_1m as exp(dalpha L/2).
# Unified ceiling (Delta = 3 Gamma, Gamma_n = nu_r/3): W = (1 - x) (Gamma_n/2) s/(1+s), s = min(1, 8 (eta/eps)^2),
# i.e. the carrier eps*Omega_p <= Omega_c/3 limits the sideband saturation; x = ALPHA eta_sp^2.
# ---------------------------------------------------------------------------
DB = 10 / math.log(10)
NU_SIL = 1.45
KAP_CAP = 2 * (NU_SIL**2 + 1) / (2 * math.sqrt(NU_SIL**2 - 1)) / (4 * math.pi**2)
A15_CAP = KAP_CAP * u01**2 * LAM_D1**2 / 15e-6**3     # bare silica capillary, ~350 dB/m
A15_IC = 0.070 / DB * (20 / 15) ** 3                  # kagome: 70 dB/km at 40 um core diameter
T_INIT = 10e-6        # radial temperature at the start of radial cooling
KB = 1.380649e-23


def best_mode(a, L, eps0, A15, umax=None, capture=True):
    """Return (W, m, eta/eps, x, nu_r) maximising the unified ceiling over LP_1m."""
    s = scenario(a)
    nu_r, x0, x = s["nu_r"], s["x0"], s["rec"]
    nth = KB * T_INIT / (HBAR * nu_r)
    k01 = A15 * (15e-6 / a) ** 3 / u01**2
    best = (0.0, 0, 0.0)
    for m, u in enumerate(u1, 1):
        if umax and u > umax:
            break
        if capture and 0.4 * (a / (u * x0)) ** 2 < 3 * nth:
            break
        eps = eps0 * math.exp(min(50.0, k01 * (u * u - u01 * u01) * L / 2))
        if eps > 1:
            continue
        r = G(u) * x0 / a / eps
        sat = min(1.0, 8 * r * r)
        W = max(0.0, 1 - x) * nu_r / 6 * sat / (1 + sat)
        if W > best[0]:
            best = (W, m, r)
    return best + (x, nu_r)


def agrid(n=160, lo=3.5e-6, hi=80e-6):  # 3.5 um: smallest hollow core in use (PBG)
    return [10 ** (math.log10(lo) + i * (math.log10(hi) - math.log10(lo)) / n) for i in range(n + 1)]


def optimum(L, eps0, A15, umax=None):
    return max((best_mode(a, L, eps0, A15, umax) + (a,) for a in agrid()), key=lambda r: r[0])


print("\nCore-size trade-off (nu_z = 300 kHz, Delta = 3 Gamma, T_init = 10 uK)")
print(f"LP01 loss at a = 15 um: capillary {A15_CAP*DB:.0f} dB/m, kagome-like {A15_IC*DB:.3f} dB/m")
for name, A15 in (("kagome-like", A15_IC), ("capillary", A15_CAP)):
    for L in (0.05, 0.3):
        for eps0 in (0.03, 0.1):
            out = []
            for a in (3.5e-6, 7e-6, 15e-6, 30e-6, 60e-6):
                W, m, r, x, nr = best_mode(a, L, eps0, A15)
                out.append(f"a={a*1e6:4.1f}: W={W:6.0f} m={m} r={r:4.2f}")
            W, m, r, x, nr, ao = optimum(L, eps0, A15)
            print(f"{name:11s} L={L*100:3.0f}cm eps0={eps0:.2f} | " + " | ".join(out) +
                  f" || a_opt={ao*1e6:.1f} W={W:.0f} m={m}")
print(f"recoil factor x = {ALPHA}*eta_sp^2 reaches 0.5 at a = "
      f"{[round(a*1e6,1) for a in agrid() if scenario(a)['rec'] >= 0.5][0]} um")

# number of guided modes with alpha L < 1 (Weyl count, both polarisations): N ~ u_max^2/2
for name, A15 in (("kagome-like", A15_IC), ("capillary", A15_CAP)):
    for a in (3.5e-6, 15e-6, 30e-6):
        for L in (0.05, 0.3):
            a01 = A15 * (15e-6 / a) ** 3
            uo = u01 * math.sqrt(1.5 / (a01 * L))
            print(f"{name:11s} a={a*1e6:4.1f}um L={L*100:3.0f}cm: N(alpha L<1) ~ {u01**2/(a01*L)/2:8.1f}"
                  f"  u_opt={uo:7.2f}")

with open(os.path.join(OUT, "tradeoff.dat"), "w") as f:
    f.write("a nur x Wic5a Wic30a Wic5b Wic30b Wcap5a Wcap30a Wcap5b Wcap30b\n")
    for a in agrid():
        rs = [best_mode(a, L, e, A15) for A15 in (A15_IC, A15_CAP) for e in (0.03, 0.1) for L in (0.05, 0.3)]
        rs = [rs[0], rs[1], rs[2], rs[3], rs[4], rs[5], rs[6], rs[7]]
        f.write(f"{a*1e6:.3f} {rs[0][4]/2/math.pi/1e3:.3f} {rs[0][3]:.3f} " +
                " ".join(f"{max(r[0], 1.0):.2f}" for r in rs) + "\n")

# optimum core radius versus fibre quality (LP01 loss at a = 15 um)
with open(os.path.join(OUT, "aopt.dat"), "w") as f:
    f.write("A15dB aopt5a W5a aopt30a W30a aopt5b W5b aopt30b W30b\n")
    for i in range(0, 41):
        A15 = 10 ** (-3 + i * 0.15) / DB
        cols = []
        for eps0 in (0.03, 0.1):
            for L in (0.05, 0.3):
                W, m, r, x, nr, ao = optimum(L, eps0, A15)
                cols += [ao * 1e6, W]
        f.write(f"{A15*DB:.4g} " + " ".join(f"{c:.3f}" for c in cols) + "\n")
print("aopt.dat:")
print(open(os.path.join(OUT, "aopt.dat")).read())

# NANF-type fibre: only LP11 guided; tolerable excess LP11 loss over L from eps(L) <= 2.83 eta (node ceiling kept)
print("\nNANF: max excess LP11 loss (dB over L) keeping eta/eps >= 0.354, and node rate")
for a in (12e-6, 15e-6, 17.5e-6):
    s = scenario(a)
    row = []
    for eps0 in (0.01, 0.03):
        budget = 2 * math.log(s["eta"] / 0.354 / eps0) * DB
        row.append(f"eps0={eps0}: {budget:5.1f} dB (={budget/0.1:5.0f} dB/m over 10 cm)")
    W = (1 - s["rec"]) * s["nu_r"] / 12
    print(f"a={a*1e6:4.1f}um nu_r=2pi x {s['nu_r']/2/math.pi/1e3:.1f} kHz eta={s['eta']:.4f} W={W:.0f}/s | " + " | ".join(row))
