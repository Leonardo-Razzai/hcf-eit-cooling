// Build the talk "Radial EIT cooling in a hollow-core fibre" (16:9, 13.33 x 7.5 in).
// Images come from make_assets.py (presentation/assets). Usage:
//   python make_assets.py && node build_deck.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const sharp = require("sharp");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const fa = require("react-icons/fa6");

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
  pres.title = "Radial EIT cooling of lattice-trapped atoms in a hollow-core fibre";
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
  async function icon(Comp, color) {
    const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(Comp, { color, size: 256 }));
    const buf = await sharp(Buffer.from(svg)).png().toBuffer();
    return "image/png;base64," + buf.toString("base64");
  }
  async function iconDot(slide, Comp, x, y, d, fill, name) {
    slide.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill }, line: { color: fill }, objectName: name + " circle" });
    const p = d * 0.22;
    slide.addImage({ data: await icon(Comp, "#FFFFFF"), x: x + p, y: y + p, w: d - 2 * p, h: d - 2 * p, objectName: name + " icon", altText: name });
  }
  function card(slide, x, y, w, h, name, fill = C.background2) {
    slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: fill }, line: { color: fill }, objectName: name });
  }
  // Physics notation in plain text: "n_ss", "Γ_n", "x_{s}" become subscript runs;
  // ν is set in Cambria italic (Calibri's ν looks like a v).
  function rich(str, base = {}) {
    const out = [];
    for (const part of String(str).split(/(_\{[^}]+\}|_[A-Za-z0-9]+|ν)/)) {
      if (!part) continue;
      if (part === "ν") out.push({ text: "ν", options: { ...base, fontFace: "Cambria", italic: true } });
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
  function bullets(items, opts = {}) {
    return items.map((it, i) => {
      const runs = expand(Array.isArray(it) ? it : [{ text: it }]);
      return runs.map((r, j) => ({ text: r.text, options: { ...(r.options || {}),
        ...(j === 0 ? { bullet: { indent: 18 } } : {}),
        ...(j === runs.length - 1 && i < items.length - 1 ? { breakLine: true } : {}),
        paraSpaceAfter: opts.space ?? 8 } }));
    }).flat();
  }
  function cite(slide, t) {
    text(slide, t, MX, 6.62, 9.5, 0.3, { fontSize: 10, color: C.accent5, italic: true });
  }

  // ------------------------------------------------------------ layouts
  pres.defineSlideMaster({
    title: "TITLE",
    background: { color: HEX.dk2 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 1.7, w: 7.4, h: 2.2, fontSize: 40, bold: true,
          color: C.background1, valign: "bottom", align: "left", margin: 0 }, text: "" } },
      { placeholder: { options: { name: "body", type: "body", x: MX, y: 4.1, w: 7.4, h: 1.2, fontSize: 20,
          color: C.accent6, valign: "top", align: "left", margin: 0 }, text: "" } },
    ],
  });
  pres.defineSlideMaster({
    title: "CONTENT",
    background: { color: HEX.lt1 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 0.32, w: W - 2 * MX, h: 0.85, fontSize: 32, bold: true,
          color: C.text2, valign: "middle", align: "left", margin: 0 }, text: "" } },
      { text: { text: "Radial EIT cooling in a hollow-core fibre", options: { x: MX, y: 7.0, w: 6, h: 0.3, fontSize: 10, color: C.accent1, margin: 0 } } },
    ],
    slideNumber: { x: W - MX - 0.6, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: C.accent1, align: "right" },
  });
  pres.defineSlideMaster({
    title: "DARK",
    background: { color: HEX.dk2 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: MX, y: 0.32, w: W - 2 * MX, h: 0.85, fontSize: 32, bold: true,
          color: C.background1, valign: "middle", align: "left", margin: 0 }, text: "" } },
    ],
    slideNumber: { x: W - MX - 0.6, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: C.accent6, align: "right" },
  });
  const content = (sec) => pres.addSlide({ masterName: "CONTENT", sectionTitle: sec });

  // ============================================================ 1 title
  pres.addSection({ title: "Introduction" });
  let s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Introduction" });
  s.addText("Radial EIT cooling of atoms in a hollow-core fibre", { placeholder: "title" });
  s.addText("An LP₁₁ probe puts the atoms at a node — towards long-lived spin-wave memories", { placeholder: "body" });
  text(s, "Leonardo Razzai", MX, 6.25, 6, 0.4, { fontSize: 16, color: C.background1 });
  {
    const [pw, ph] = await size("mode_maps.png");
    const iw = 3.9, ih = (iw * ph) / pw, cy = (H - ih) / 2 - 0.2;
    card(s, 8.45, cy - 0.25, iw + 0.4, ih + 0.5, "title figure card", C.background1);
    await img(s, "mode_maps.png", 8.65, cy, iw, ih, { name: "LP01 and LP11 cross-sections" });
  }
  s.addNotes("Topic: cooling the radial motion of atoms trapped in a hollow-core fibre lattice. The key element is a probe in the LP11 mode, whose node sits exactly where the atoms are. We end with what this enables for spin-wave quantum memories.");

  // ============================================================ 2 problem
  s = content("Introduction");
  s.addText("Radial motion stays hot in a fibre lattice", { placeholder: "title" });
  await img(s, "lattice.png", MX, 1.45, 6.8, 3.5, { name: "LP01 lattice sketch" });
  text(s, bullets([
    "All guided beams propagate along z: no transverse photon momentum",
    "Both EIT beams in LP₀₁: radial coupling only ~ (x₀/w)² ≈ 10⁻⁴",
    "Axial motion already ground-state cooled in fibres (Raman sideband cooling)",
  ]), MX, 5.1, 7.6, 1.35, { fontSize: 15 });
  const stats = [["300 kHz", "axial ν_z (1064 nm lattice)"], ["7 kHz", "radial ν_r, 40× softer (2a = 30 µm)"], ["29", "⟨n_r⟩ at 10 µK, per radial axis"]];
  stats.forEach(([big, small], i) => {
    const y = 1.5 + i * 1.65;
    card(s, 8.6, y, 4.1, 1.4, "stat card " + (i + 1));
    text(s, big, 8.85, y + 0.15, 3.7, 0.7, { fontSize: 32, bold: true, color: i === 1 ? C.accent1 : C.accent5, fontFace: "Cambria" });
    text(s, small, 8.85, y + 0.85, 3.7, 0.4, { fontSize: 14 });
  });
  cite(s, "Axial RSC in a fibre lattice: Leong et al., Nat. Commun. 11, 5295 (2020)");
  s.addNotes("In an LP01 lattice the axial confinement is strong (300 kHz) and axial sideband cooling works. The radial trap is 40 times softer and every guided beam carries momentum only along z. With both EIT beams in LP01 the ratio Omega_p/Omega_c is independent of r, so the radial coupling is only second order in x0/w.");

  // ============================================================ 3 setup
  pres.addSection({ title: "Mechanism" });
  s = content("Mechanism");
  s.addText("Proposed setup", { placeholder: "title" });
  await img(s, "setup.png", MX, 1.3, W - 2 * MX, 5.0, { name: "experimental setup schematic" });
  text(s, "Vertical antiresonant fibre (2a ≈ 30 µm) loaded from a MOT; 1064 nm LP₀₁ standing wave; D1 coupling beam in LP₀₁ and phase-locked probe (+6.8 GHz) shaped into LP₁₁ by an SLM, co-propagating from the bottom end.",
    MX, 6.35, W - 2 * MX, 0.55, { fontSize: 14, color: C.text2 });
  s.addNotes("Implementation I. MOT above a vertical antiresonant fibre, transfer into a 1064 nm LP01 standing wave from two counter-propagating beams. The two EIT beams come from one D1 laser system, 6.8 GHz apart. The probe is converted to LP11 with an SLM or a 0/pi phase plate, combined on a PBS with the LP01 coupling beam and launched from the same end, so the axial wavevector difference is negligible (axial Lamb-Dicke parameter < 1e-3). Inset: the coupling/lattice mode peaks on axis, the probe has its node there.");

  // ============================================================ 4 key idea
  s = content("Mechanism");
  s.addText("Key idea: put the atoms at a probe node", { placeholder: "title" });
  await img(s, "mode_maps.png", MX, 1.45, 6.0, 4.6, { valign: "middle", name: "mode maps" });
  card(s, 7.0, 1.5, 5.73, 1.9, "equation card");
  await img(s, "eq_gradient.png", 7.25, 1.62, 5.2, 0.85, { align: "left", valign: "middle", name: "gradient coupling equation" });
  await img(s, "eq_eta.png", 7.25, 2.5, 5.2, 0.8, { align: "left", valign: "middle", name: "effective Lamb-Dicke parameter" });
  text(s, bullets([
    [{ text: "No carrier: ", options: { bold: true } }, { text: "the probe drives only Δn = ±1 — the fibre analogue of an ion at a standing-wave node" }],
    [{ text: "Self-aligned: ", options: { bold: true } }, { text: "the LP₁₁ node sits on the axis by symmetry, at the lattice centre" }],
    [{ text: "Small but harmless: ", options: { bold: true } }, { text: "η⊥ ≈ 0.02–0.04; without a carrier it costs only probe power (nW)" }],
  ], { space: 12 }), 7.0, 3.65, 5.73, 2.9, { fontSize: 16 });
  cite(s, "Cirac, Blatt, Zoller & Phillips, Phys. Rev. A 46, 2668 (1992)");
  s.addNotes("Near the axis the LP11 field is linear in x with slope G11/a at equal power, G11 = 3.49. Expanding x = x0(b + b^dagger) gives a pure sideband coupling with effective Lamb-Dicke parameter eta_perp = G11 x0/a. This is exactly the node configuration of Cirac et al.: carrier and its recoil diffusion vanish. The small eta_perp would cripple a carrier-present scheme, but here there is nothing that limits the probe intensity except the saturation of the narrow EIT resonance.");

  // ============================================================ 5 EIT
  s = content("Mechanism");
  s.addText("EIT cooling on the radial sideband", { placeholder: "title" });
  await img(s, "levels.png", MX, 1.4, 6.6, 3.6, { name: "level scheme" });
  text(s, bullets([
    "Coupling (Δ = 3Γ) dresses |g_c⟩ into a narrow state |+⟩: δ₊ ≈ Ω_c²/4Δ, Γ_n ≈ ΓΩ_c²/4Δ²",
    "δ₊ = ν_r: red sideband on the narrow peak, blue sideband suppressed by (Γ/4Δ)²",
  ]), MX, 5.15, 6.6, 1.4, { fontSize: 15 });
  await img(s, "fano.png", 7.5, 1.4, 5.2, 3.3, { name: "Fano absorption spectrum" });
  card(s, 7.5, 4.95, 5.23, 1.3, "Fano equation card");
  await img(s, "eq_fano.png", 7.65, 5.05, 4.93, 1.1, { valign: "middle", name: "Fano profile equation" });
  cite(s, "Morigi, Eschner & Keitel, Phys. Rev. Lett. 85, 4458 (2000)");
  s.addNotes("Standard EIT cooling on the D1 line, but with the motional coupling provided by the probe gradient. Choosing Omega_c^2 = 4 nu_r (Delta + nu_r) places the narrow dressed state exactly on the red sideband. Because the probe has no carrier, the Fano zero at delta = 0 is not even needed for the carrier; it suppresses heating only through the blue sideband.");

  // ============================================================ 6 simple model
  s = content("Mechanism");
  s.addText("Simple model: the narrow line sets the rate, not η⊥", { placeholder: "title" });
  card(s, MX, 1.5, 6.3, 3.2, "rate equations card");
  await img(s, "eq_W.png", MX + 0.25, 1.75, 5.8, 0.85, { align: "left", valign: "middle", name: "rate equation" });
  await img(s, "eq_n.png", MX + 0.25, 2.85, 5.8, 0.85, { align: "left", valign: "middle", name: "steady-state equation" });
  await img(s, "eq_s.png", MX + 0.25, 3.8, 5.8, 0.75, { align: "left", valign: "middle", name: "saturation parameter" });
  text(s, "Rate bounded by Γ_n/2 ≈ Γν_r/2Δ for any η⊥; the gradient only sets the probe power, P ∝ (a/G₁₁x₀)² s.",
    MX, 4.95, 6.3, 0.9, { fontSize: 15, italic: true, color: C.text2 });
  const rows6 = [
    [fa.FaArrowsRotate, C.accent1, "κ ≈ 2 photons per cycle", "Decays into |g_c⟩ are re-pumped by the coupling beam: recoil per cycle h = καη_sp² (0.40 in a 30 µm core)"],
    [fa.FaGaugeHigh, C.accent2, "Only the red sideband saturates", "The off-resonant blue sideband does not, so n_ss grows as (1 + s)"],
    [fa.FaBolt, C.accent4, "25–50× faster than with a carrier", "A carrier-present probe must stay perturbative: W ∝ η⊥², only ≈ 26 s⁻¹ in case A"],
  ];
  for (let i = 0; i < rows6.length; i++) {
    const [ic, col, head, body] = rows6[i];
    const y = 1.5 + i * 1.6;
    await iconDot(s, ic, 7.3, y + 0.05, 0.62, col, head);
    text(s, head, 8.1, y, 4.63, 0.4, { fontSize: 17, bold: true, color: C.text2 });
    text(s, body, 8.1, y + 0.45, 4.63, 1.0, { fontSize: 14 });
  }
  s.addNotes("Two corrections to the textbook rate equation matter in a fibre. First, every decay into g_c is followed by further scattering of the coupling beam, so a cycle costs kappa = 1/b_p photons, about 2; with the soft radial trap eta_sp approaches 0.7 and the recoil term h = 0.4 in case A. Second, saturation reduces only the resonant red sideband. The saturation parameter is s = 2(Omega_R/Gamma_n)^2 with Omega_R the effective Raman coupling; this is what breaks the rate picture at s ~ 1.");

  // ============================================================ 7 platforms
  pres.addSection({ title: "Results" });
  s = content("Results");
  s.addText("Two representative fibres", { placeholder: "title" });
  const hdr = (t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } });
  const rowsT = [
    ["Core radius a", "15 µm", "3.5 µm"],
    ["Radial frequency ν_r", "7.15 kHz", "30.6 kHz"],
    ["x₀ ; η⊥ = G₁₁x₀/a", "90 nm ; 0.021", "44 nm ; 0.043"],
    ["η_sp ; recoil per cycle h", "0.71 ; 0.40", "0.34 ; 0.094"],
    ["Ω_c/2π ; Γ_n/2π", "0.70 MHz ; 2.4 kHz", "1.45 MHz ; 10.2 kHz"],
    ["Probe power at s = 1", "8 nW", "0.4 nW"],
    ["⟨n_x⟩ at 10 µK", "29", "6.3"],
  ];
  s.addTable([[hdr("Parameter"), hdr("A: antiresonant"), hdr("B: bandgap")],
    ...rowsT.map((r, i) => r.map((c, j) => ({ text: rich(c, { bold: j === 0 }), options: { fill: { color: i % 2 ? C.background1 : C.background2 } } })))],
    { x: MX, y: 1.5, w: 8.0, colW: [3.3, 2.35, 2.35], fontSize: 15, color: C.text1, border: { type: "none" }, rowH: 0.55, valign: "middle", objectName: "case table" });
  card(s, 9.1, 1.5, 3.63, 4.4, "assumptions card");
  text(s, "Common settings", 9.35, 1.7, 3.2, 0.4, { fontSize: 17, bold: true, color: C.text2 });
  text(s, bullets(["⁸⁷Rb, D1 line, closed Λ, Δ = 3Γ", "1064 nm LP₀₁ lattice, ν_z = 300 kHz (U₀ = 530 µK)", "One radial axis, axial motion frozen", "Exact Bessel profiles, dipole-pattern recoil"], { space: 10 }),
    9.35, 2.2, 3.2, 3.6, { fontSize: 14 });
  s.addNotes("Case A is the mature platform: a large-core antiresonant fibre. Case B is a short small-core bandgap section, which raises nu_r four-fold and reduces the recoil per cycle. Note the nW probe powers. The master equation includes exact Bessel mode profiles, the LP01 admixture and trap offset as imperfections, and spontaneous-emission recoil with the dipole pattern.");

  // ============================================================ 8 rate
  s = content("Results");
  s.addText("Ground-state cooling with nanowatts of probe", { placeholder: "title" });
  await img(s, "saturation.png", MX, 1.45, 8.3, 4.6, { name: "rate and occupation versus saturation" });
  const st8 = [["1.3×10³ s⁻¹", "case A, s = 1: n_ss = 0.13 (8 nW)", C.accent1], ["1.2×10⁴ s⁻¹", "case B, s = 1: n_ss = 0.04 (0.4 nW)", C.accent2], ["0.017 / 0.009", "weak-probe floor n_ss (A / B)", C.text2]];
  st8.forEach(([big, small, col], i) => {
    const y = 1.5 + i * 1.6;
    card(s, 9.3, y, 3.43, 1.35, "result card " + (i + 1));
    text(s, big, 9.5, y + 0.15, 3.1, 0.6, { fontSize: 26, bold: true, color: col, fontFace: "Cambria" });
    text(s, small, 9.5, y + 0.8, 3.1, 0.45, { fontSize: 13 });
  });
  text(s, "At weak probe the simple model matches the master equation within 5 % (with κ = 2). Above s ≈ 0.1 it overestimates W by up to 65 % and underestimates n_ss by up to 4×.",
    MX, 6.0, 8.3, 0.6, { fontSize: 14, italic: true, color: C.text2 });
  s.addNotes("Symbols master equation, dashed simple model. The maximum rate, 3.3e3 and 1.8e4 per second, is reached at s about 5 but costs a higher occupation (0.4 and 0.2); s about 1 is a sensible working point. The rates are 30-45 percent of the resolution ceiling (1-h)Gamma_n/2.");

  // ============================================================ 9 robustness
  s = content("Results");
  s.addText("Robustness: probe purity is the critical requirement", { placeholder: "title" });
  await img(s, "robust.png", MX, 1.35, W - 2 * MX, 3.95, { name: "robustness plots" });
  const tk = [
    ["Mode impurity", "Rate halves at ε ≈ 1.4η⊥ ≈ 0.03: need < 10⁻³ LP₀₁ power. A trap offset acts as ε_eff = η⊥x_s/x₀."],
    ["Two-photon error", "With a clean node, positive errors cost rate but not temperature; with ε = 0.05 the carrier heats."],
    ["Trap-frequency mismatch", "W > 70 % of its maximum for ±10 % in ν_r: sets the uniformity along the strip."],
  ];
  tk.forEach(([h1, b1], i) => {
    const x = MX + i * 4.1;
    text(s, h1, x, 5.45, 3.85, 0.35, { fontSize: 16, bold: true, color: C.text2 });
    text(s, b1, x, 5.85, 3.85, 0.8, { fontSize: 13 });
  });
  s.addNotes("Left: the simple-model threshold epsilon < 2.8 eta_perp is a factor two too optimistic; squares show that a trap offset is exactly equivalent to a mode admixture. Middle: this is the node advantage of Cirac et al.: with no carrier, a two-photon error only shifts the sidebands on the Fano profile. Right: lasers set for a design frequency; plus or minus 10 percent in nu_r corresponds to about 8 dB of lattice-beam loss or 0.8 dB of coupling-beam loss across the strip.");

  // ============================================================ 10 dynamics
  s = content("Results");
  s.addText("Cooling from 10 µK: the thermal tail sets the time", { placeholder: "title" });
  await img(s, "dynamics.png", MX, 1.45, 8.3, 4.6, { name: "cooling dynamics" });
  card(s, MX, 5.6, 8.3, 0.95, "dynamics takeaway");
  text(s, "Each decay into |g_c⟩ lengthens the cycle: phonon removal per atom saturates, so a hot thermal tail cools linearly in time, not exponentially.",
    MX + 0.25, 5.68, 7.8, 0.8, { fontSize: 15, italic: true, color: C.text2, valign: "middle" });
  text(s, bullets([
    [{ text: "Saturated removal: ", options: { bold: true } }, { text: "1.3–1.9 phonons/ms (A), 10–12 (B) — only 30–45 % of (1−h)Γ_n/2" }],
    [{ text: "Calibrated rate model ", options: { bold: true } }, { text: "(c_sat = 2.25) tracks the master equation within 8 %" }],
    [{ text: "Staged Δ = Γ → 3Γ: ", options: { bold: true } }, { text: "A: 97.5 % in n = 0 after 150 ms (⟨n⟩ = 1.4, carried by a 0.4 % hot tail); B: ⟨n⟩ = 0.04 in 2.4 ms" }],
    [{ text: "LP₁₂ = LP₁₁: ", options: { bold: true } }, { text: "the rate depends on s, not on the mode order" }],
  ], { space: 12 }), 9.2, 1.5, 3.53, 5.0, { fontSize: 14 });
  s.addNotes("Left: from n = 2 the master equation relaxes far more slowly than the single-exponential simple model. Each decay into g_c lengthens the cooling cycle, so the phonon-removal rate per atom saturates well below the naive two-level bound. A rate model with the full Franck-Condon factors and a calibrated saturation reproduces the dynamics. Right: from 10 microkelvin the hot thermal tail dominates in case A; a broad first stage helps.");

  // ============================================================ 11 initial temperature
  s = content("Results");
  s.addText("Cooling time grows linearly with the initial temperature", { placeholder: "title" });
  await img(s, "temperature.png", MX, 1.45, 8.3, 4.0, { name: "cooling versus initial temperature" });
  card(s, MX, 5.6, 8.3, 0.95, "temperature takeaway");
  text(s, "Above ⟨n⟩ ≈ 3 the sideband saturates: the phonon-removal rate is constant, so the cooling time scales as k_{B}T_r/U_0 — and a small core gains twice.",
    MX + 0.25, 5.68, 7.8, 0.8, { fontSize: 15, italic: true, color: C.text2, valign: "middle" });
  text(s, bullets([
    [{ text: "Scan: ", options: { bold: true } }, { text: "k_{B}T_r/U_0 = 10⁻³–10⁻¹ (0.5–53 µK), 20 log steps; harmonic trap within 4 %" }],
    [{ text: "90 % in n = 0: ", options: { bold: true } }, { text: "≈ 2.5 s × k_{B}T_r/U_0 (A), ≈ 80 ms × k_{B}T_r/U_0 (B)" }],
    [{ text: "Saturated removal: ", options: { bold: true } }, { text: "1.9 phonons/ms (A), 12 (B)" }],
    [{ text: "Hot tail: ", options: { bold: true } }, { text: "recoil spreads energy by η_sp√n per photon; ⟨n⟩ lags behind P(n = 0)" }],
    [{ text: "Rate model ", options: { bold: true } }, { text: "underestimates the master-equation drift above n ≈ 15: times are conservative" }],
  ], { space: 10 }), 9.2, 1.5, 3.53, 5.0, { fontSize: 14 });
  s.addNotes("We scan the initial radial temperature in units of the lattice depth, from one thousandth to one tenth of U0, which is the range where the radial trap stays harmonic and the calibrated rate model applies. At the same k_B T over U0 the small core is four times colder in units of its trap quantum, and it removes six times more phonons per millisecond. Above a few quanta the red sideband saturates, the phonon-removal rate is constant and the cooling time grows linearly with temperature. The ground-state fraction converges faster than the mean, because photon recoil lets a small fraction of atoms diffuse to high energy. Beyond the calibration range the rate model is conservative: at 20 to 30 quanta the master equation removes phonons 1.2 to 2.4 times faster.");

  // ============================================================ 12 D1 versus D2
  s = content("Results");
  s.addText("D1 or D2? The excited hyperfine structure decides", { placeholder: "title" });
  await img(s, "lines.png", MX, 1.45, 8.3, 4.0, { name: "D1 versus D2 cooling" });
  card(s, MX, 5.6, 8.3, 0.95, "lines takeaway");
  text(s, "D2's excited levels lie 12–44 Γ from the Λ, D1's 140 Γ: their light shifts barely touch D1 but can wash out D2's narrow resonance.",
    MX + 0.25, 5.68, 7.8, 0.8, { fontSize: 15, italic: true, color: C.text2, valign: "middle" });
  text(s, bullets([
    [{ text: "Model: ", options: { bold: true } }, { text: "all excited hyperfine levels in the master equation; two-photon detuning re-optimised; bars span both signs of the Raman interference" }],
    [{ text: "D1 is robust: ", options: { bold: true } }, { text: "< 5 % change after a 0.1 ν_r retune" }],
    [{ text: "Best Λ: D1 F′=1, probe on F=2 ", options: { bold: true } }, { text: "(b_p = 5/6): 2.1–2.3× the paper rate, n_ss ≈ 0.04 (A)" }],
    [{ text: "D2 F′=1 is fragile: ", options: { bold: true } }, { text: "×2 sign dependence, no cooling for Δ ≥ 8Γ" }],
    [{ text: "Next: ", options: { bold: true } }, { text: "Zeeman model — dark states of F=2→F′=1" }],
  ], { space: 9 }), 9.2, 1.5, 3.53, 5.0, { fontSize: 14 });
  s.addNotes("We compared all closed Lambda systems on both lines. In the three-level picture the only difference is the branching ratio: when most decays return to the probe ground state, fewer photons are scattered per removed phonon and the cycle is shorter, so the rate doubles. Including all excited hyperfine levels in the master equation, the extra light shift moves the narrow resonance; we re-optimise the two-photon detuning as an experiment would. On D1 the other level is 814 MHz away and changes the results by a few percent. On D2 the levels are 72 to 267 MHz apart: the F'=2 Lambda still works, but the F'=1 Lambda, the D2 option with the favourable branching ratio, depends strongly on the sign of the interference through the other levels, which only a Zeeman-resolved model can fix, and stops cooling at larger detuning. Our recommendation is D1 via F'=1 with the probe on F=2, pending a check of Zeeman dark states.");

  // ============================================================ 13 model
  pres.addSection({ title: "Model" });
  s = content("Model");
  s.addText("Simple model versus master equation", { placeholder: "title" });
  const colw = (W - 2 * MX - 0.4) / 2;
  const cols = [
    [fa.FaCircleCheck, C.accent4, "What we added to the simple model", [
      "κ = 1/b_p photons per cycle: without it n_ss is 40 % too low (case A)",
      "Saturation of the red sideband only → n_ss ∝ (1 + s)",
      "Rate model with exact Franck–Condon factors and calibrated saturation (c_sat = 2.25): dynamics within 8 %",
      "Full master equation: Bessel profiles, dipole-pattern recoil, impurity, offset, detuning errors",
    ]],
    [fa.FaTriangleExclamation, C.accent2, "Where the simple model fails", [
      "s ≳ 0.1: Raman coupling Ω_R ~ Γ_n splits the narrow line — W up to 65 % too high, n_ss up to 4× too low",
      "Purity threshold ε ≲ 2.8η⊥ is 2× too optimistic",
      "Thermal capture: predicts a single exponential, far too fast",
      "Two-photon error ≲ −0.45ν_r: predicts heating where the master equation stays finite",
    ]],
  ];
  for (let i = 0; i < cols.length; i++) {
    const [ic, col, head, items] = cols[i];
    const x = MX + i * (colw + 0.4);
    card(s, x, 1.45, colw, 3.75, head);
    await iconDot(s, ic, x + 0.3, 1.65, 0.6, col, head);
    text(s, head, x + 1.05, 1.72, colw - 1.3, 0.5, { fontSize: 19, bold: true, color: C.text2 });
    text(s, bullets(items, { space: 10 }), x + 0.3, 2.45, colw - 0.6, 2.7, { fontSize: 15 });
  }
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: MX, y: 5.45, w: W - 2 * MX, h: 1.05, rectRadius: 0.12,
    fill: { color: C.background1 }, line: { color: C.accent2, width: 1.25 }, objectName: "Raman coupling outline" });
  await img(s, "eq_s.png", MX + 0.3, 5.55, 5.2, 0.85, { align: "left", valign: "middle", name: "saturation parameter equation" });
  text(s, "The rate picture breaks when the sideband Raman coupling Ω_R reaches the narrow linewidth Γ_n.",
    6.4, 5.62, 6.1, 0.75, { fontSize: 15, italic: true, color: C.text2, valign: "middle" });
  s.addNotes("The simple model is a good design tool at weak probe: Fano placement, scaling with core size and detuning. At the working points that maximise the rate it is not reliable, because the coherent sideband Raman coupling becomes comparable to the narrow linewidth. Use the master equation, or the calibrated rate model for large phonon numbers.");

  // ============================================================ 14 missing physics
  s = content("Model");
  s.addText("Beyond both models: what is still missing", { placeholder: "title" });
  const grid = [
    [fa.FaAtom, C.accent1, "D1 Zeeman structure", "Leakage out of the Λ raises κ; case A (h = 0.4) is pushed towards h → 1"],
    [fa.FaArrowsLeftRight, C.accent2, "Second radial axis", "One lobe cools one axis: alternate LP₁₁ orientations or use lobe rotation"],
    [fa.FaWaveSquare, C.accent4, "Vector modes", "TE₀₁/TM₀₁/HE₂₁ beating rotates the lobes and the polarisation along the fibre"],
    [fa.FaLayerGroup, C.accent5, "Axial–radial interplay", "One Fano peak cannot serve ν_r and ν_z ≈ 40ν_r: alternate with axial RSC"],
    [fa.FaUsers, C.accent6, "Collective effects", "Optical depth, rescattering, probe depletion in case B (absorbed ≤ 0.8 nW)"],
    [fa.FaTemperatureHigh, C.accent3, "Technical heating", "Lattice intensity and pointing noise, higher-order-mode content of the lattice beam"],
  ];
  const gw = (W - 2 * MX - 0.6) / 3;
  for (let i = 0; i < grid.length; i++) {
    const [ic, col, head, body] = grid[i];
    const x = MX + (i % 3) * (gw + 0.3), y = 1.5 + Math.floor(i / 3) * 2.55;
    card(s, x, y, gw, 2.3, head);
    await iconDot(s, ic, x + 0.25, y + 0.25, 0.55, col, head);
    text(s, head, x + 0.95, y + 0.3, gw - 1.15, 0.45, { fontSize: 17, bold: true, color: C.text2 });
    text(s, body, x + 0.25, y + 1.0, gw - 0.5, 1.2, { fontSize: 14 });
  }
  s.addNotes("These are not captured by either the simple model or our master equation, and set the agenda for the next round of simulations: the full D1 Zeeman structure first, then both radial axes.");

  // ============================================================ 15 implementations
  pres.addSection({ title: "Experiment" });
  s = content("Experiment");
  s.addText("Two implementations", { placeholder: "title" });
  const implStat = [["10³ s⁻¹", "cooling rate, 150 ms from 10 µK"], ["10⁴ s⁻¹", "cooling rate, 3 ms from 10 µK"]];
  const impl = [
    ["I", "Large-core antiresonant fibre", C.accent1, [
      "2a ≈ 30 µm NANF, low loss at 795 and 1064 nm",
      "MOT loading above a vertical fibre",
      "LP₁₁ probe from an SLM, purity checked by S² imaging",
      "2–8 nW probe, W ≈ 10³ s⁻¹, 97 % in n = 0 after 150 ms",
      "Most mature; weakness: slow capture of the thermal tail",
    ]],
    ["II", "Small-core bandgap section", C.accent2, [
      "1–5 cm of 2a ≈ 7 µm fibre, moving-lattice loading",
      "Probe injected from the near end (lossy LP₁₁ mode)",
      "4× higher ν_r, 4× smaller recoil per cycle",
      "W ≈ 10⁴ s⁻¹, ground state from 10 µK in 3 ms",
      "Costs: loading volume, modal loss, wall proximity",
    ]],
  ];
  for (let i = 0; i < impl.length; i++) {
    const [num, head, col, items] = impl[i];
    const x = MX + i * (colw + 0.4);
    card(s, x, 1.45, colw, 4.6, "implementation " + num);
    s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 1.65, w: 0.65, h: 0.65, fill: { color: col }, line: { color: col }, objectName: "number " + num });
    text(s, num, x + 0.3, 1.65, 0.65, 0.65, { fontSize: 20, bold: true, color: C.background1, align: "center", valign: "middle", fontFace: "Cambria" });
    text(s, head, x + 1.15, 1.75, colw - 1.4, 0.5, { fontSize: 20, bold: true, color: C.text2 });
    text(s, bullets(items, { space: 10 }), x + 0.3, 2.55, colw - 0.6, 2.5, { fontSize: 15 });
    text(s, implStat[i][0], x + 0.3, 5.05, 2.6, 0.7, { fontSize: 30, bold: true, color: col, fontFace: "Cambria", valign: "middle" });
    text(s, implStat[i][1], x + 2.9, 5.05, colw - 3.2, 0.7, { fontSize: 14, color: C.text2, valign: "middle" });
  }
  text(s, "Key requirement for both: LP₀₁ power content of the probe below 10⁻³ along the atom section (differential-loss budget ≈ 9 dB for a 1 % launch admixture).",
    MX, 6.2, W - 2 * MX, 0.6, { fontSize: 15, italic: true, color: C.text2 });
  s.addNotes("Implementation I follows existing Rb and Sr fibre-lattice experiments (Bajcsy, Blatt, Okaba, Leong). Implementation II trades loading volume and wall distance for a ten-fold faster cooling. Thermometry: Raman sideband spectroscopy with one LP11 and one LP01 Raman beam gives the same carrier-free gradient coupling.");

  // ============================================================ 16 memory mechanism
  pres.addSection({ title: "Quantum memory" });
  s = content("Quantum memory");
  s.addText("Spin-wave memory: radial temperature sets the dephasing", { placeholder: "title" });
  await img(s, "memory.png", MX, 1.4, 6.3, 4.6, { name: "light-shift spread versus radial occupation" });
  card(s, 7.2, 1.45, 5.53, 3.25, "memory equations card");
  await img(s, "eq_spinwave.png", 7.4, 1.6, 5.13, 0.85, { align: "left", valign: "middle", name: "spin-wave equation" });
  await img(s, "eq_lightshift.png", 7.4, 2.55, 5.13, 0.95, { align: "left", valign: "middle", name: "differential light shift" });
  await img(s, "eq_spread.png", 7.4, 3.6, 5.13, 0.9, { align: "left", valign: "middle", name: "light-shift spread" });
  text(s, bullets([
    "Collinear storage in the fibre: λ_sw ≈ 4.4 cm, insensitive to atomic motion",
    "Dephasing from the lattice differential light shift sampled by radial motion",
    "10 µK: σ_δ ≈ 15 Hz (T₂* ~ 20 ms) → after EIT: 0.01–0.1 Hz (T₂* ~ seconds)",
  ]), 7.2, 4.9, 5.53, 1.7, { fontSize: 14 });
  cite(s, "Kuhr et al., PRA 72, 023406 (2005); Zhao et al., Nat. Phys. 5, 100 (2009); Dudin et al., PRA 87, 031801 (2013)");
  s.addNotes("Order-of-magnitude estimate. The hyperfine coherence sees the differential light shift of the 1064 nm lattice, delta_0 about 0.8 kHz at U0/h = 11 MHz. An atom with radial potential energy V sees a shift reduced by V/U0, so the inhomogeneous width is delta_0 times the thermal radial energy over U0. At 10 microkelvin this is about 15 Hz; after radial ground-state cooling the thermal part drops by two to three orders of magnitude, and the zero-point part is common to all atoms. T2* = 1/(pi sigma) follows the Kuhr et al. analysis. Other limits then take over.");

  // ============================================================ 17 memory benefits
  s = content("Quantum memory");
  s.addText("What a cold atom strip brings to a fibre memory", { placeholder: "title" });
  const mem = [
    [fa.FaRuler, C.accent1, "Motion-insensitive spin wave", "Collinear signal and control: λ_sw ≈ 4.4 cm; axial motion already in the Lamb–Dicke regime"],
    [fa.FaCircleNodes, C.accent4, "Uniform coupling, huge optical depth", "All atoms on the mode axis; OD ~ 10³ in fibres → efficiency 1 − O(1/d)"],
    [fa.FaClock, C.accent2, "Compensation that works", "With radial thermal spread removed, light-shift compensation can reach the seconds-to-minute regime of free-space lattice memories"],
    [fa.FaMagnet, C.accent5, "Remaining limits", "Magnetic noise (use m_F = 0 clock states), lattice heating, collisions, retrieval mode matching"],
  ];
  const mw = (W - 2 * MX - 0.4) / 2;
  for (let i = 0; i < mem.length; i++) {
    const [ic, col, head, body] = mem[i];
    const x = MX + (i % 2) * (mw + 0.4), y = 1.45 + Math.floor(i / 2) * 2.45;
    card(s, x, y, mw, 2.2, head);
    await iconDot(s, ic, x + 0.3, y + 0.3, 0.6, col, head);
    text(s, head, x + 1.1, y + 0.38, mw - 1.35, 0.45, { fontSize: 18, bold: true, color: C.text2 });
    text(s, body, x + 1.1, y + 0.9, mw - 1.35, 1.2, { fontSize: 15 });
  }
  cite(s, "Blatt et al., Opt. Lett. 39, 446 (2014); Gorshkov et al., PRL 98, 123601 (2007); Peters et al., Opt. Express 28, 5340 (2020)");
  s.addNotes("Hollow-core memories so far used freely expanding cold atoms or warm vapour. A three-dimensionally ground-state-cooled lattice strip combines the optical depth of the fibre with the coherence of free-space lattice memories. Open questions: the magic or compensated configuration in the fibre, and heating from lattice intensity noise over second-long storage.");

  // ============================================================ 18 summary
  pres.addSection({ title: "Summary" });
  s = pres.addSlide({ masterName: "DARK", sectionTitle: "Summary" });
  s.addText("Summary and next steps", { placeholder: "title" });
  const sum = [
    [fa.FaBullseye, "LP₁₁ probe = carrier-free node", "Gradient coupling η⊥ = G₁₁x₀/a, self-aligned with the lattice"],
    [fa.FaBolt, "10³–10⁴ s⁻¹ with nW", "n_ss = 0.02–0.2 for an ideal Λ; small cores strongly favoured"],
    [fa.FaTriangleExclamation, "Purity and model limits", "LP₀₁ content < 10⁻³; rate models must be checked against the master equation"],
    [fa.FaMemory, "Long-lived fibre memories", "Radial ground state removes the thermal light-shift dephasing of the spin wave"],
  ];
  for (let i = 0; i < sum.length; i++) {
    const [ic, head, body] = sum[i];
    const x = MX + (i % 2) * (mw + 0.4), y = 1.5 + Math.floor(i / 2) * 1.85;
    await iconDot(s, ic, x, y, 0.65, C.accent1, head);
    text(s, head, x + 0.9, y + 0.02, mw - 1.0, 0.45, { fontSize: 19, bold: true, color: C.background1 });
    text(s, body, x + 0.9, y + 0.55, mw - 1.0, 0.9, { fontSize: 15, color: C.accent6 });
  }
  text(s, "Next: full D1 Zeeman structure and both radial axes in the simulations; a first experiment in a large-core antiresonant fibre with Raman-sideband thermometry.",
    MX, 5.45, W - 2 * MX, 0.9, { fontSize: 16, italic: true, color: C.background1 });
  s.addNotes("Take-home messages and the next steps.");

  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);          // write THEME colours/fonts into the deck theme
  console.log("wrote", OUT);
}

main().catch((e) => { console.error(e); process.exit(1); });
