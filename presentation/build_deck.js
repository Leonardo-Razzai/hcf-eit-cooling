// Build the talk "Radial EIT cooling in a hollow-core fibre" (16:9, 13.33 x 7.5 in).
// Plain layout: title, figure, a few statements, one reference line. No icons.
// Images come from make_assets.py (presentation/assets). Usage:
//   python make_assets.py && node build_deck.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const sharp = require("sharp");

const { applyTheme } = require("./apply_theme.js");
const OUT = path.join(__dirname, "radial_eit_cooling_talk.pptx");
const A = (f) => path.join(__dirname, "assets", f);

const THEME = {
  name: "Fibre Night",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1A2233", lt1: "FFFFFF", dk2: "0F1B2D", lt2: "EEF2F7",
    accent1: "0072B2", accent2: "D55E00", accent3: "E69F00",
    accent4: "009E73", accent5: "8C2D19", accent6: "56B4E9",
    hlink: "0072B2", folHlink: "6B4C9A",
  },
};
const HEX = THEME.colors;
const W = 13.333, H = 7.5, MX = 0.6;

async function main() {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
  pres.title = "Radial ground-state cooling of atoms in a hollow-core fibre";
  pres.author = "Leonardo Razzai";
  const C = pres.SchemeColor;

  // ------------------------------------------------------------ helpers
  const sizes = {};
  async function size(f) {
    if (!sizes[f]) { const m = await sharp(A(f)).metadata(); sizes[f] = [m.width, m.height]; }
    return sizes[f];
  }
  // place an image inside a box, keeping its aspect ratio
  async function img(slide, f, x, y, maxW, maxH, { align = "center", valign = "top", name } = {}) {
    const [pw, ph] = await size(f);
    let w = maxW, h = (maxW * ph) / pw;
    if (h > maxH) { h = maxH; w = (maxH * pw) / ph; }
    const dx = align === "center" ? (maxW - w) / 2 : align === "right" ? maxW - w : 0;
    const dy = valign === "middle" ? (maxH - h) / 2 : valign === "bottom" ? maxH - h : 0;
    slide.addImage({ path: A(f), x: x + dx, y: y + dy, w, h, objectName: name || f, altText: name || f });
    return { x: x + dx, y: y + dy, w, h };
  }
  // Physics notation in plain text: "n_ss", "Γ_n", "x_{s}" become subscript runs, "n^{-1/h}" superscript;
  // ν is set in Cambria italic (Calibri's ν looks like a v).
  function rich(str, base = {}) {
    const out = [];
    for (const part of String(str).split(/(\^\{[^}]+\}|_\{[^}]+\}|_[A-Za-z0-9]+|ν)/)) {
      if (!part) continue;
      if (part === "ν") out.push({ text: "ν", options: { ...base, fontFace: "Cambria", italic: true } });
      else if (part.startsWith("^")) out.push({ text: part.slice(2, -1), options: { ...base, superscript: true } });
      else if (part.startsWith("_")) out.push({ text: part.replace(/^_\{?|\}$/g, ""), options: { ...base, subscript: true } });
      else out.push({ text: part, options: { ...base } });
    }
    return out;
  }
  const expand = (runs) => runs.flatMap((r) => rich(r.text, r.options || {}));
  function text(slide, t, x, y, w, h, opts = {}) {
    slide.addText(typeof t === "string" ? rich(t) : t,
      { x, y, w, h, isTextBox: true, fontSize: 16, color: C.text1, valign: "top", margin: 0, ...opts });
  }
  // a bullet may be a string, or [bold lead, rest]
  function bullets(items, opts = {}) {
    return items.map((it, i) => {
      const runs = expand(Array.isArray(it) ? [{ text: it[0], options: { bold: true } }, { text: it[1] }] : [{ text: it }]);
      return runs.map((r, j) => ({ text: r.text, options: { ...(r.options || {}),
        ...(j === 0 ? { bullet: { indent: 16 } } : {}),
        ...(j === runs.length - 1 && i < items.length - 1 ? { breakLine: true } : {}),
        paraSpaceAfter: opts.space ?? 10 } }));
    }).flat();
  }
  function cite(slide, t) {
    text(slide, t, MX, 6.6, W - 2 * MX - 0.8, 0.3, { fontSize: 11, color: "5B6475", italic: true });
  }
  function table(slide, header, rows, x, y, colW, opts = {}) {
    const hdr = header.map((t) => ({ text: rich(t, { bold: true }), options: { color: C.background1, fill: { color: C.text2 } } }));
    const body = rows.map((r, i) => r.map((c) => ({ text: rich(c, { bold: false }),
      options: { fill: { color: i % 2 ? HEX.lt1 : HEX.lt2 } } })));
    slide.addTable([hdr, ...body], { x, y, w: colW.reduce((a, b) => a + b, 0), colW, fontSize: opts.fontSize || 14,
      color: C.text1, border: { type: "none" }, rowH: opts.rowH || 0.42, valign: "middle", margin: 0.06, objectName: opts.name || "table" });
  }

  // ------------------------------------------------------------ layouts
  pres.defineSlideMaster({
    title: "TITLE",
    background: { color: HEX.dk2 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 1.6, w: 7.6, h: 2.3, fontSize: 38, bold: true,
          color: C.background1, valign: "bottom", align: "left", margin: 0 }, text: "" } },
      { placeholder: { options: { name: "body", type: "body", x: MX, y: 4.1, w: 7.6, h: 1.2, fontSize: 20,
          color: C.accent6, valign: "top", align: "left", margin: 0 }, text: "" } },
    ],
  });
  pres.defineSlideMaster({
    title: "CONTENT",
    background: { color: HEX.lt1 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 0.32, w: W - 2 * MX, h: 0.85, fontSize: 30, bold: true,
          color: C.text2, valign: "middle", align: "left", margin: 0 }, text: "" } },
      { line: { x: MX, y: 1.2, w: W - 2 * MX, h: 0, line: { color: HEX.accent1, width: 1 } } },
      { text: { text: "Radial EIT cooling in a hollow-core fibre", options: { x: MX, y: 7.0, w: 6, h: 0.3, fontSize: 10, color: "5B6475", margin: 0 } } },
    ],
    slideNumber: { x: W - MX - 0.6, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: "5B6475", align: "right" },
  });
  pres.defineSlideMaster({
    title: "DARK",
    background: { color: HEX.dk2 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 0.32, w: W - 2 * MX, h: 0.85, fontSize: 30, bold: true,
          color: C.background1, valign: "middle", align: "left", margin: 0 }, text: "" } },
    ],
    slideNumber: { x: W - MX - 0.6, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: C.accent6, align: "right" },
  });
  const content = (sec) => pres.addSlide({ masterName: "CONTENT", sectionTitle: sec });
  // standard layouts: figure left (8.0 in) and statements right, or figure on top
  const FIG_W = 7.9, TX = MX + FIG_W + 0.35, TW = W - MX - TX;

  // ============================================================ 1 title
  pres.addSection({ title: "Introduction" });
  let s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Introduction" });
  s.addText("Radial ground-state cooling of atoms in a hollow-core fibre", { placeholder: "title" });
  s.addText("EIT cooling with an LP₁₁ probe: what is known, and what is new", { placeholder: "body" });
  text(s, "Leonardo Razzai", MX, 6.25, 6, 0.4, { fontSize: 16, color: C.background1 });
  await img(s, "mode_maps.png", 8.6, 2.0, 4.1, 3.4, { valign: "middle", name: "LP01 and LP11 cross-sections" });
  s.addNotes("The topic is cooling the radial motion of atoms trapped in a lattice inside a hollow-core fibre to the motional ground state. The key element is a probe in the LP11 mode, whose node sits exactly on the fibre axis, where the atoms are. I will first place this in the context of EIT cooling and of the closest related work, then say what is new, and then show the main results.");

  // ============================================================ 2 problem
  s = content("Introduction");
  s.addText("Atoms in a hollow-core fibre: the radial motion stays hot", { placeholder: "title" });
  await img(s, "lattice.png", MX, 1.5, FIG_W, 4.6, { valign: "middle", name: "LP01 lattice sketch" });
  text(s, bullets([
    "A cm-long atomic column in a guided 1064 nm lattice: high optical depth without a cavity",
    "Axial motion: ν_z = 300 kHz, ground-state cooled by Raman sideband cooling",
    "Radial motion: ν_r ≈ 7 kHz for a = 15 µm, ⟨n⟩ ≈ 29 at 10 µK",
    "Every guided beam pushes along z: LP₀₁ light alone cannot cool the radial motion",
  ], { space: 14 }), TX, 1.6, TW, 4.8, { fontSize: 17 });
  cite(s, "Fibre lattices: Okaba et al., Nat. Commun. 5, 4096 (2014); Leong et al., Nat. Commun. 11, 5295 (2020)");
  s.addNotes("In an LP01 lattice the axial confinement is strong and axial sideband cooling works. The radial trap is forty times softer, so a 10 microkelvin sample has about 30 radial quanta, and every guided beam carries momentum only along the axis. With both EIT beams in LP01 the radial coupling is only second order in x0 over the mode radius.");

  // ============================================================ 3 EIT cooling
  s = content("Introduction");
  s.addText("EIT cooling in one slide", { placeholder: "title" });
  await img(s, "levels.png", MX, 1.45, 6.6, 3.5, { name: "level scheme" });
  await img(s, "fano.png", 7.5, 1.45, 5.2, 3.3, { name: "Fano absorption spectrum" });
  text(s, bullets([
    "A strong coupling beam dresses the Λ system: the probe absorption gets a narrow peak next to a dark-state zero",
    "Put the narrow peak on the red sideband: absorption removes a phonon, the carrier and the blue sideband are suppressed",
    "Established for trapped ions and ion strings; recently used to load, image and cool ⁸⁷Rb tweezer arrays",
  ], { space: 10 }), MX, 5.15, W - 2 * MX, 1.4, { fontSize: 16 });
  cite(s, "Morigi, Eschner & Keitel, PRL 85, 4458 (2000); Roos et al., PRL 85, 5547 (2000); Lechner et al., PRA 93, 053401 (2016); Qiu et al., PRX Quantum 7, 020358 (2026)");
  s.addNotes("Standard EIT cooling. The coupling beam, detuned by Delta, creates a narrow dressed resonance at the light shift delta plus. Setting delta plus equal to the trap frequency places it on the red sideband. Theory by Morigi, Eschner and Keitel, first experiment by Roos and Blatt with a calcium ion, extended to long ion strings, and recently by the Lukin and Vuletic groups to load and image rubidium atoms in tweezer arrays while cooling axial and radial motion.");

  // ============================================================ 4 related work
  s = content("Introduction");
  s.addText("Closest related work", { placeholder: "title" });
  table(s, ["Work", "System", "What it showed"], [
    ["Cirac et al. 1992", "ion at a standing-wave node", "no carrier at the node: pure sideband coupling"],
    ["Morigi et al.; Roos et al. 2000", "trapped ion", "EIT ground-state cooling, theory and experiment"],
    ["Xing & Mehta 2025", "ion in structured light (theory)", "carrier-free EIT at a node; rate independent of the gradient"],
    ["Xing et al. 2025", "Ca⁺ in a phase-stable standing wave", "first EIT cooling at a node: n̄ ≈ 0.05, modes over 5 MHz in 150 µs"],
    ["Qiu, …, Lukin, Vuletić 2026", "⁸⁷Rb tweezer arrays", "EIT loading (68 %) and imaging (99.7 %), cools axial and radial motion"],
    ["Leong et al. 2020", "Rb lattice in a hollow-core fibre", "axial ground state by Raman sideband cooling"],
    ["Dechant et al. 2015, 2016", "Sisyphus cooling", "heavy-tailed, non-Boltzmann energy distributions"],
  ], MX, 1.5, [3.4, 3.7, 5.03], { fontSize: 14, rowH: 0.56, name: "related work table" });
  cite(s, "Cirac et al., PRA 46, 2668 (1992); Xing & Mehta, PR Applied 24, 014034 (2025); Xing et al., arXiv:2512.10900; Dechant et al., PRL 115, 173006 (2015), PRE 94, 022151 (2016)");
  s.addNotes("The idea of a node is old: Cirac and co-workers showed in 1992 that an ion at a standing-wave node sees no carrier. Xing and Mehta worked out carrier-free EIT cooling at a node of structured light for ions in integrated photonics, and their group demonstrated it in 2025 with a calcium ion in a phase-stable standing wave. On the neutral-atom side, the Harvard-MIT groups used EIT on rubidium in tweezers. In fibres, axial ground-state cooling has been shown. Heavy-tailed distributions are known for Sisyphus cooling. So the concept of node cooling is not ours; what is new is the regime and the platform.");

  // ============================================================ 5 scheme
  pres.addSection({ title: "Scheme" });
  s = content("Scheme");
  s.addText("This scheme: probe in LP₁₁, atoms on its node", { placeholder: "title" });
  await img(s, "setup.png", MX, 1.4, FIG_W, 5.1, { valign: "middle", name: "setup schematic" });
  await img(s, "eq_eta.png", TX, 1.6, TW, 0.7, { align: "left", valign: "middle", name: "effective Lamb-Dicke parameter" });
  text(s, bullets([
    "Lattice and coupling beam in LP₀₁, probe in LP₁₁, all guided from one end",
    "The LP₁₁ node lies on the axis by symmetry: no carrier, only the field gradient couples to radial motion",
    "η⊥ ≈ 0.02–0.04 is small, but without a carrier it only sets the probe power: nanowatts",
  ], { space: 14 }), TX, 2.6, TW, 3.9, { fontSize: 16 });
  s.addNotes("Coupling beam and lattice in LP01, the probe shaped into LP11 by an SLM or a phase plate and launched from the same end. Near the axis the LP11 field is linear in x, with slope G11 over a, so the probe couples only through its gradient, with effective Lamb-Dicke parameter eta perp equal to G11 x0 over a. This is the fibre version of the ion at a node. Being carrier-free, a small eta perp costs only probe power.");

  // ============================================================ 6 known vs new
  s = content("Scheme");
  s.addText("What is known, and what is new", { placeholder: "title" });
  table(s, ["Known (ions, tweezers)", "New here (atoms in a hollow-core fibre)"], [
    ["carrier-free EIT cooling at a node", "all-guided radial cooling with LP₁₁: self-aligned node, nW probe"],
    ["recoil is a small correction (h ≈ 0.01)", "recoil-dominated regime: h = 0.1–1 in a soft radial trap"],
    ["near-thermal steady states", "steady state with a power-law tail, P(n) ∝ n^{−1/h}, while P₀ > 0.95"],
    ["cooling rate from ⟨n⟩", "⟨n⟩-based rate fails; the ground-state rate W₀ is the robust measure"],
    ["heavy tails in Sisyphus cooling", "heavy tails in resolved-sideband cooling, from n ≈ 3"],
    ["one probe cools one motional mode", "both radial axes coupled by recoil: a cold 2D state needs h < 1/2"],
    ["—", "universal map in (ν_r/ω_rec, b_p); ⁸⁷Rb level and Zeeman structure; memory budget"],
  ], MX, 1.5, [4.6, 7.53], { fontSize: 15, rowH: 0.6, name: "known versus new table" });
  s.addNotes("Left, what we take from the ion work: the node concept and the fact that the rate is set by the narrow line, not by the gradient. Right, what is specific to atoms in a fibre and new to our knowledge. The radial trap is soft, so the recoil of each scattered photon is comparable to the phonon being removed. This produces a qualitatively different steady state, with a power-law tail, it breaks the usual definition of the cooling rate, and it couples the two radial axes. We map it out with one dimensionless parameter, h.");

  // ============================================================ 7 h
  pres.addSection({ title: "Results" });
  s = content("Results");
  s.addText("The soft trap: recoil per removed phonon, h", { placeholder: "title" });
  await img(s, "hscale.png", MX, 1.4, W - 2 * MX, 3.6, { name: "h scale" });
  await img(s, "eq_h.png", MX, 5.2, 4.6, 0.95, { align: "left", valign: "middle", name: "definition of h" });
  text(s, bullets([
    "Fibre atoms sit 10–50× above ions and tweezers, inside the range where the steady state changes",
    "Thresholds: h = 1/3 (variance), 1/2 (mean), 1 (ground-state cooling)",
    "h falls with a smaller core, a deeper lattice, or a Λ that returns to the probe state (b_p → 1)",
  ], { space: 8 }), 5.5, 5.15, W - MX - 5.5, 1.45, { fontSize: 15 });
  s.addNotes("h is the recoil energy fed back per removed phonon, in units of the phonon. alpha is the fraction of emission recoil along one radial axis, b_p the branching ratio back to the probe ground state, so 1 over b_p photons are scattered per cooling cycle. For a trapped ion h is about 0.01, for a rubidium tweezer about 0.03; in the fibre it is 0.1 to 0.4. Because nu_r scales as nu_z over a, h scales as a over b_p nu_z: a small core or a deep lattice helps.");

  // ============================================================ 8 rates
  s = content("Results");
  s.addText("Ground-state cooling with nanowatts of probe", { placeholder: "title" });
  await img(s, "saturation.png", MX, 1.4, FIG_W, 4.9, { valign: "middle", name: "rate and occupation versus saturation" });
  text(s, bullets([
    ["Case A, a = 15 µm, s = 1: ", "W = 1.6×10³ s⁻¹ (ground-state rate W₀ = 2.2×10³ s⁻¹), P₀ = 0.95, 8 nW"],
    ["Case B, a = 3.5 µm, s = 1: ", "W = 1.2×10⁴ s⁻¹ (W₀ = 1.4×10⁴ s⁻¹), ⟨n⟩ = 0.04, 0.4 nW"],
    ["Simple rate model: ", "valid only at weak probe and only if it counts 1/b_p ≈ 2 photons per cycle; at s ≈ 1 it overestimates W by up to 40 % and underestimates ⟨n⟩ up to 5×"],
  ], { space: 14 }), TX, 1.6, TW, 4.8, { fontSize: 16 });
  cite(s, "Converged master equation: excited state eliminated, banded coherences, up to 400 Fock states");
  s.addNotes("Symbols: master equation, converged in the number of Fock states, which needed adiabatic elimination of the excited state and a banded solver. Dashed: the standard rate model. s is the saturation of the narrow line, s = 2 times Omega_R over Gamma_n squared. The rate model fails at s near 1 because the Raman coupling splits the narrow line. W0 is the rate of the ground-state population; the rate from the mean occupation is not well defined, as the next slides show.");

  // ============================================================ 9 tail
  s = content("Results");
  s.addText("New: a heavy-tailed steady state", { placeholder: "title" });
  await img(s, "tail.png", MX, 1.4, 6.9, 5.1, { valign: "middle", name: "steady-state Fock distribution" });
  text(s, bullets([
    "P(n) falls as a power law, n^{−1/h}, over two decades; the exponent follows from h alone (fit −2.34, predicted −2.48)",
    "Case A, s = 1: P₀ = 0.95, yet 0.3 % of atoms above n = 10 carry a third of ⟨n⟩",
    "So ⟨n⟩ and the ⟨n⟩-based rate depend on the truncation; P₀ and W₀ do not",
    "Same law as in Sisyphus cooling, here in resolved-sideband ground-state cooling",
  ], { space: 14 }), 7.85, 1.6, W - MX - 7.85, 4.8, { fontSize: 16 });
  cite(s, "Sisyphus analogue: Dechant, Kessler & Barkai, PRL 115, 173006 (2015)");
  s.addNotes("This is the main new physics. Each removed phonon costs about two photons, and their recoil adds energy proportional to n on average, while the removal rate saturates. The balance gives a power law with exponent minus 1 over h. A thermal state with the same mean occupation, dashed, falls off exponentially. The tail ends where the narrow line no longer resolves the sideband. In practice the ground-state population is high, but the mean is carried by a few hot atoms, so it converges only with hundreds of Fock states.");

  // ============================================================ 10 regimes
  s = content("Results");
  s.addText("New: three cooling regimes on one curve", { placeholder: "title" });
  await img(s, "regimes.png", MX, 1.4, FIG_W, 5.0, { valign: "middle", name: "removal rate regimes" });
  text(s, bullets([
    "The phonon-removal rate depends on s and n only through s·n",
    "R ∝ n, then a saturated plateau, then a fast channel above s·n ≈ 2(ν_r/Γ_n)² ≈ 18",
    "The plateau makes hot atoms cool linearly in time, not exponentially",
    "Fast channel: drift grows 17× from n = 20 to 200 while diffusion grows 1.8×; Sisyphus-like reading still to be tested",
  ], { space: 14 }), TX, 1.6, TW, 4.8, { fontSize: 16 });
  s.addNotes("Drift of the mean occupation from Fock states, master equation, case A. Curves for different s collapse onto one curve against s times n, the saturation of the n to n minus 1 sideband. At large s n the Raman splitting exceeds the trap frequency and a much faster channel opens; our working hypothesis is dark-state following over the smooth mode profile, a Sisyphus-like mechanism. The sign-of-detuning test is still to be done.");

  // ============================================================ 11 map
  s = content("Results");
  s.addText("New: one cooling map for any atom and trap", { placeholder: "title" });
  await img(s, "map.png", MX, 1.35, W - 2 * MX, 4.4, { name: "cooling map" });
  text(s, bullets([
    "Only ν_r/ω_rec, b_p, Δ/Γ and s matter; Γ/ν_r and x₀/a change results by < 0.3 %",
    "h = 1/3 separates clean ground states from heavy tails; Δ = 3Γ, s = 1 shown",
    "For ⁸⁷Rb, D1 F′=1 with the probe on F=2 (b_p = 5/6): P₀ = 0.975 and 1.5× the rate of D1 F′=2",
  ], { space: 8 }), MX, 5.8, W - 2 * MX, 1.0, { fontSize: 15 });
  s.addNotes("Converged master-equation map in the plane of trap stiffness in recoil units and branching ratio. Left, ground-state population; right, the ground-state rate in units of nu_r, shown only where the state is cold. Lines of constant h follow the colour contours. The open symbols are the configurations we considered for rubidium. Ions and tweezers would sit far to the right; we plan to add them as points.");

  // ============================================================ 12 two axes
  s = content("Results");
  s.addText("New: both radial axes", { placeholder: "title" });
  await img(s, "two_axes.png", MX, 1.4, 6.9, 5.1, { valign: "middle", name: "two-axis ground state" });
  text(s, bullets([
    "One probe pattern cools one motional mode; recoil heats both axes",
    "The second axis doubles the effective recoil: a cold 2D state needs h < 1/2",
    "Case A, Δ = 3Γ: P₀₀ = 0.76 with b_p = 1/2, 0.93 with b_p = 5/6; 0.98 at Δ = 6Γ",
    "Since h ∝ a/ν_z, a deeper lattice relaxes the limit",
  ], { space: 14 }), 7.85, 1.6, W - MX - 7.85, 4.8, { fontSize: 16 });
  cite(s, "Rate chain on P(n_x, n_y), checked against the exact 2D master equation (0.927 vs 0.931)");
  s.addNotes("Open symbols: one axis; filled: both axes, cooled by alternating the LP11 orientation. Below h of one third the second axis costs a factor of about two in the excited population; above one half the two-dimensional ground state is lost at three Gamma. A larger detuning helps. A single scalar vortex beam does not cool both modes; a vector mode might, which is still open.");

  // ============================================================ 13 capture
  s = content("Results");
  s.addText("Cooling from 10 µK: time grows with temperature", { placeholder: "title" });
  await img(s, "capture.png", MX, 1.4, 6.9, 5.1, { valign: "middle", name: "cooling time versus temperature" });
  text(s, bullets([
    "Above n̄ ≈ 3 the sideband saturates: the cooling time grows linearly with the initial temperature",
    "90 % in n = 0: 2.6–255 ms in case A, below 8 ms in case B, for 0.5–53 µK",
    "From 10 µK: 97.5 % in n = 0 after 150 ms (A); ⟨n⟩ = 0.04 after 2.4 ms (B)",
    "The rate model is conservative: the master equation is about 2× faster from n̄ = 30",
  ], { space: 14 }), 7.85, 1.6, W - MX - 7.85, 4.8, { fontSize: 16 });
  s.addNotes("Calibrated rate model with the full trap depth in Fock space. The time to reach 90 percent in the ground state grows linearly with temperature because the phonon removal rate per atom is constant in the plateau. The small core is about thirty times faster. A master-equation check from thirty quanta gives 21 milliseconds against 47 from the rate model in case A, so these times are upper bounds.");

  // ============================================================ 14 Rb
  pres.addSection({ title: "Rubidium" });
  s = content("Rubidium");
  s.addText("⁸⁷Rb: which transition, and what the Zeeman structure costs", { placeholder: "title" });
  table(s, ["Λ configuration (case A, Δ = 3Γ, s = 1)", "W₀ (s⁻¹)", "P₀"], [
    ["D1 F′=2, b_p = 1/2", "2.2×10³", "0.953"],
    ["D1 F′=1, probe on F=2, b_p = 5/6", "3.2–3.6×10³", "0.974–0.979"],
    ["D2 F′=2", "2.1–2.2×10³", "0.94–0.955"],
    ["D2 F′=1, probe on F=1", "2.0–4.2×10³ (sign)", "fails above 8Γ"],
  ], MX, 1.5, [3.3, 2.3, 1.7], { fontSize: 13, rowH: 0.5, name: "transition table" });
  await img(s, "rb.png", 8.15, 1.35, W - MX - 8.15, 5.3, { valign: "middle", name: "Rb transitions and Zeeman cost" });
  text(s, bullets([
    ["All excited hyperfine levels: ", "D1 is robust (within 1 %); D2 F′=1 depends on the sign of the interference"],
    ["All 8 ground sublevels: ", "a fixed polarisation leaves dark states, so population must be remixed; a weak repumper or a Zeeman Λ inside F=2 (as in Qiu et al.) cools at about 1/3 of the closed-Λ rate, W₀ ≈ 10³ s⁻¹ at P₀ ≥ 0.97"],
    ["Open: ", "Zeeman structure with both axes, and the real polarisation of the fibre modes"],
  ], { space: 10 }), MX, 4.2, 7.3, 2.4, { fontSize: 14 });
  s.addNotes("Figure a: ground-state rate against the population left outside the ground state, as the detuning goes from 2 to 12 Gamma; up and left is better. Figure b: with all Zeeman sublevels, the best rate reachable at each ground-state population for each remixing scheme, against the closed lambda. First the excited hyperfine levels. On D1 the second excited level is 814 MHz away and changes nothing; on D2 the levels are close and the F'=1 configuration depends on the interference sign. D1 F'=1 with the probe on F=2 is the best option, because five sixths of the decays return to the probe state. With all Zeeman sublevels, any fixed polarisation leaves dark states independent of motion; remixing by a repumper costs about a factor of three in rate. The Zeeman lambda used by Qiu and co-workers in tweezers performs about the same.");

  // ============================================================ 15 requirements
  s = content("Rubidium");
  s.addText("Practical requirements", { placeholder: "title" });
  table(s, ["Requirement", "Value", "Reason"], [
    ["Probe mode purity", "LP₀₁ power < 10⁻³", "rate halves at an LP₀₁ amplitude ε ≈ 1.3 η⊥ ≈ 0.03"],
    ["Trap centring", "offset below ≈ x₀ (0.1 µm in case A)", "an offset x_s acts as an admixture ε = η⊥ x_s/x₀"],
    ["Radial frequency uniformity", "± 10 % along the column", "keeps W above 70 % of its maximum"],
    ["Probe power", "0.4–8 nW", "set by s, not by the small gradient"],
    ["Mode order", "LP₁₁ or LP₁₂", "same rate at equal s; a higher mode only lowers the power"],
  ], MX, 1.5, [2.2, 2.25, 2.85], { fontSize: 13, rowH: 0.7, name: "requirements table" });
  await img(s, "robust.png", 8.15, 1.35, W - MX - 8.15, 5.3, { valign: "middle", name: "purity and trap-frequency tolerance" });
  s.addNotes("Figure a: the maximum rate against the LP01 admixture in units of eta perp; squares are trap offsets converted with epsilon = eta perp x_s / x0, which fall on the same curve; dashed is the simple model, a factor two too optimistic. Figure b: lasers set for the design trap frequency while the actual one varies, as it would along a lossy fibre. Mode purity is the critical experimental requirement, as for ions: a small LP01 admixture restores a carrier. The simple estimate of the tolerance is a factor two too optimistic. A displaced trap is equivalent to an admixture. The trap frequency must be uniform along the column to within ten percent, which sets the tolerable lattice-beam loss.");

  // ============================================================ 16 memory
  pres.addSection({ title: "Application" });
  s = content("Application");
  s.addText("Application: a spin-wave memory in the fibre", { placeholder: "title" });
  await img(s, "memory.png", MX, 1.35, W - 2 * MX, 4.05, { name: "memory optical depth and storage time" });
  text(s, bullets([
    ["Optical depth is the bottleneck: ", "OD = 0.302 r µm² × ρ_s; 90 % efficiency needs ρ_s ≈ 1200 µm⁻² on the clock Λ: a ≲ 4 µm for a 1 cm column at one atom per site, or several atoms per site"],
    ["Storage: ", "cool deep, store shallow, with an echo: 2–11 s, limited by vacuum, then by lattice heating"],
    ["Ground state buys ", "shallow traps, uniform coupling and 60× fewer motional jumps; an echo alone removes static thermal dephasing"],
  ], { space: 6 }), MX, 5.5, W - 2 * MX, 1.35, { fontSize: 14 });
  s.addNotes("Theory only. Left: optical depth per centimetre of column for one atom per site, as a function of core radius; rho_s is the number of atoms per square micrometre of mode area. A 15 micrometre core gives only 7 to 15 per centimetre on a usable lambda scheme, so efficient storage needs a small core, a longer column or several atoms per site. Right: storage model with lattice scattering and a one percent intensity spread along the column. Lowering the lattice a hundred times after cooling and using an echo gives seconds, limited by the background gas inside the fibre, which is the least certain input. The radial ground state is not mainly a coherence-time gain; it allows shallow traps and uniform coupling.");

  // ============================================================ 17 summary
  pres.addSection({ title: "Summary" });
  s = pres.addSlide({ masterName: "DARK", sectionTitle: "Summary" });
  s.addText("Summary", { placeholder: "title" });
  const sumItems = [
    ["All-guided radial ground-state cooling: ", "an LP₁₁ probe puts the atoms at a node; 10³–10⁴ s⁻¹ with nanowatts"],
    ["Recoil-dominated regime: ", "the steady state has a power-law tail, n^{−1/h}; P₀ and W₀, not ⟨n⟩, measure the cooling"],
    ["One map for any atom and trap: ", "h < 1/3 for a clean state, h < 1/2 for both radial axes, h ∝ a/(b_p ν_z)"],
    ["⁸⁷Rb: ", "D1 F′=1 with the probe on F=2; the Zeeman structure costs about 3× in rate"],
    ["Memory: ", "small cores for optical depth; seconds of storage after lowering the lattice"],
  ];
  text(s, bullets(sumItems.map(([a, b]) => [a, b]), { space: 16 }).map((r) => ({ ...r, options: { ...r.options, color: r.options.bold ? HEX.lt1 : HEX.accent6 } })),
    MX, 1.5, W - 2 * MX, 3.9, { fontSize: 19 });
  text(s, "Next: optimal saturation ramp from thermal states, real fibre-mode polarisation, Zeeman structure with both axes, memory simulations.",
    MX, 5.6, W - 2 * MX, 0.8, { fontSize: 16, italic: true, color: C.background1 });
  s.addNotes("Take-home messages. The concept of node cooling comes from ions; what is new is the platform, all beams guided in a hollow-core fibre, and the regime, where recoil is comparable to the trap quantum and changes the nature of the steady state. Next steps are listed at the bottom.");

  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);          // write THEME colours/fonts into the deck theme
  console.log("wrote", OUT);
}

main().catch((e) => { console.error(e); process.exit(1); });
