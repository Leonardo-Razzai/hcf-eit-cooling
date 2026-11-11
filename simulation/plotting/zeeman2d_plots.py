"""Analysis of the 2D Zeeman runs (studies/zeeman2d.py run, on the cluster). Reads the run JSON files only.

Writes into results/zeeman2d/ (or --out DIR/..):
  summary.txt         flags (drift, frame), best runs per scheme, comparison with the closed-Lambda
                      two-axis model (vector_twoaxes.dat, case A), convergence in N
  overview.pdf/.png   a: P00 against W00 for every run, by scheme (upper right is better);
                      b: best W00 (P00 >= 0.85) per scheme and s
  detunings.pdf/.png  W00 (P00 >= 0.85) against the two-photon detunings: TE01 (d1, d2) maps at the best
                      s and repumper, LP11 against d
  convergence.pdf     P00 and W00 against N (if convergence runs exist)
  traces.pdf          P00(t) of the best run of each scheme
Usage: python zeeman2d_plots.py [--runs DIR]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import style
from style import MM, WDOUBLE, CA, CB, CG, CL, CO, letter

R = Path(__file__).resolve().parents[1] / 'results'
style.apply()
COL = {'lp11_hf': CG, 'lp11_zl': CL, 'te01_bz': CB, 'te01_bx': CA}
LAB = {'lp11_hf': 'LP$_{11}$ switched, hyperfine $\\Lambda$', 'lp11_zl': 'LP$_{11}$ switched, Zeeman $\\Lambda$',
       'te01_bz': 'TE$_{01}$, $B\\parallel z$', 'te01_bx': 'TE$_{01}$, $B\\parallel x$'}
P_COLD = 0.85
DRIFT_MAX = 1e-3


def load_runs(d):
    runs = []
    for f in sorted(Path(d).glob('*.json')):
        try:
            runs.append(json.loads(f.read_text()))
        except json.JSONDecodeError:
            pass
    return runs


def ok(r):
    return abs(r['drift']) < DRIFT_MAX and r['frame_mismatch'] < 1e-6


def best(runs, key='W00'):
    good = [r for r in runs if ok(r) and r['P00'] >= P_COLD]
    return max(good, key=lambda r: r[key]) if good else None


def summary(runs, out):
    L = []
    L.append('2D Zeeman runs: %d files' % len(runs))
    bad = [r for r in runs if not ok(r)]
    L.append('flagged (|drift| >= %.0e or frame mismatch): %d' % (DRIFT_MAX, len(bad)))
    for r in bad[:20]:
        L.append('   %s  drift %.1e  mismatch %.1e' % (r['tag'], r['drift'], r['frame_mismatch']))
    L.append('')
    L.append('Best run per scheme (P00 >= %.2f, unflagged), main resolution:' % P_COLD)
    L.append('%-9s %4s %5s %5s %5s %6s %6s %8s %9s %3s' % ('scheme', 'case', 's', 'd1', 'd2', 'gamma', 'P00', 'W00/nu', 'W00 (1/s)', 'N'))
    ref = None
    f = R / 'vector_twoaxes.dat'
    if f.exists():
        ref = np.genfromtxt(f, names=True)
    main_N = min((r['N'] for r in runs), default=6)
    for sch in COL:
        for case in ('A', 'B'):
            sel = [r for r in runs if r['scheme'] == sch and r['case'] == case and r['N'] == main_N]
            if not sel:
                continue
            b = best(sel)
            if b is None:
                bP = max(sel, key=lambda r: r['P00'])
                L.append('%-9s %4s   no cold run; highest P00 %.3f (%s)' % (sch, case, bP['P00'], bP['tag']))
                continue
            L.append('%-9s %4s %5.2f %5.2f %5.2f %6.2f %6.3f %8.4f %9.0f %3d' % (
                sch, case, b['s'], b['d1'], b['d2'], b['gamma'], b['P00'], b['W00'], b['W00_per_s'], b['N']))
    if ref is not None:
        r = ref[(np.abs(ref['nu_wrec'] - 1.97) < 0.05) & (ref['s'] == 1.0)]
        if r.size:
            names = {0: 'switched LP11', 3: 'TE01 separate g_c'}
            L.append('')
            L.append('Closed-Lambda two-axis reference (minimal model, case A, b_p = 5/6, s = 1, N = 6):')
            for k, n in names.items():
                q = r[r['scheme'] == k]
                if q.size:
                    L.append('   %-18s P00 %.3f  W00 %.4f nu_r  (%.0f s^-1)' % (n, q['P00'][0], q['W00'][0], q['W00'][0] * 2 * np.pi * 7148.78))
    conv = {}
    for r in runs:
        conv.setdefault((r['scheme'], r['case'], r['s'], r['d1'], r['d2'], r['gamma'], r['tau']), []).append(r)
    conv = {k: sorted(v, key=lambda r: r['N']) for k, v in conv.items() if len({r['N'] for r in v}) > 1}
    if conv:
        L.append('')
        L.append('Convergence in N (same parameters):')
        for k, v in conv.items():
            L.append('   %s %s s %g d %g/%g gamma %g: ' % (k[0], k[1], k[2], k[3], k[4], k[5]) +
                     '; '.join('N %d: P00 %.4f W00 %.4f' % (r['N'], r['P00'], r['W00']) for r in v))
    txt = '\n'.join(L)
    (out / 'summary.txt').write_text(txt + '\n')
    print(txt)
    return conv


def overview(runs, out):
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE, 65 * MM), constrained_layout=True)
    a = ax[0]
    for sch, c in COL.items():
        sel = [r for r in runs if r['scheme'] == sch]
        if not sel:
            continue
        x = [r['W00'] for r in sel]
        y = [r['P00'] for r in sel]
        mk = ['o' if ok(r) else 'x' for r in sel]
        for xi, yi, m in zip(x, y, mk):
            a.plot(xi, yi, m, color=c, ms=3, mfc=c if m == 'o' else 'none')
        a.plot([], [], 'o', color=c, ms=3, label=LAB[sch])
    a.axhline(P_COLD, color=CG, lw=0.5, ls=':')
    a.set_xlabel(r'$W_{00}/\nu_r$')
    a.set_ylabel(r'$P_{00}$')
    a.set_xlim(left=0)
    a.legend(loc='lower right', fontsize=5.5)
    b = ax[1]
    ss = sorted({r['s'] for r in runs})
    width = 0.8 / max(1, len(ss))
    for j, s in enumerate(ss):
        vals = []
        for sch in COL:
            bb = best([r for r in runs if r['scheme'] == sch and r['s'] == s])
            vals.append(bb['W00'] if bb else np.nan)
        b.bar(np.arange(len(COL)) + (j - (len(ss) - 1) / 2) * width, np.nan_to_num(vals), width * 0.95,
              color=[COL[k] for k in COL], alpha=0.5 + 0.5 * j / max(1, len(ss) - 1), label='$s=%g$' % s)
    b.set_xticks(range(len(COL)), [LAB[k].replace(', ', '\n') for k in COL], fontsize=5.5)
    b.set_ylabel(r'best $W_{00}/\nu_r$ ($P_{00}\geq%.2f$)' % P_COLD)
    b.legend(loc='upper left', fontsize=5.5)
    letter(a, 'a')
    letter(b, 'b')
    fig.savefig(out / 'overview.pdf')
    fig.savefig(out / 'overview.png', dpi=200)
    plt.close(fig)


def detunings(runs, out):
    fig, ax = plt.subplots(1, 3, figsize=(WDOUBLE, 55 * MM), constrained_layout=True)
    for a, sch in zip(ax[:2], ('te01_bz', 'te01_bx')):
        sel = [r for r in runs if r['scheme'] == sch and ok(r)]
        b = best(sel)
        a.set_title(LAB[sch], fontsize=6.5)
        if b is None:
            a.text(0.5, 0.5, 'no cold run', ha='center', transform=a.transAxes)
            continue
        sel = [r for r in sel if r['s'] == b['s'] and r['gamma'] == b['gamma'] and r['case'] == b['case'] and r['N'] == b['N']]
        d1s, d2s = sorted({r['d1'] for r in sel}), sorted({r['d2'] for r in sel})
        Z = np.full((len(d2s), len(d1s)), np.nan)
        for r in sel:
            Z[d2s.index(r['d2']), d1s.index(r['d1'])] = r['W00'] if r['P00'] >= P_COLD else np.nan
        m = a.imshow(Z, origin='lower', cmap='Blues', aspect='auto',
                     extent=(min(d1s) - 0.2, max(d1s) + 0.2, min(d2s) - 0.2, max(d2s) + 0.2))
        fig.colorbar(m, ax=a, label=r'$W_{00}/\nu_r$')
        a.set_xlabel(r'$d_1/\nu_r$')
        a.set_ylabel(r'$d_2/\nu_r$')
        a.text(0.02, 1.02, '$s=%g$, $\\gamma=%g$' % (b['s'], b['gamma']), transform=a.transAxes, fontsize=5.5, va='bottom')
    a = ax[2]
    for sch in ('lp11_hf', 'lp11_zl'):
        sel = [r for r in runs if r['scheme'] == sch and ok(r)]
        for key in sorted({(r['s'], r['gamma'], r['tau']) for r in sel}):
            q = sorted([r for r in sel if (r['s'], r['gamma'], r['tau']) == key], key=lambda r: r['d1'])
            W = [r['W00'] if r['P00'] >= P_COLD else np.nan for r in q]
            a.plot([r['d1'] for r in q], W, 'o-', color=COL[sch], ms=2.5, lw=0.7, alpha=0.4 + 0.6 * (key[0] == 1.0))
        a.plot([], [], 'o-', color=COL[sch], label=LAB[sch])
    a.set_xlabel(r'$d/\nu_r$')
    a.set_ylabel(r'$W_{00}/\nu_r$ ($P_{00}\geq%.2f$)' % P_COLD)
    a.legend(fontsize=5)
    for x, l in zip(ax, 'abc'):
        letter(x, l)
    fig.savefig(out / 'detunings.pdf')
    fig.savefig(out / 'detunings.png', dpi=200)
    plt.close(fig)


def convergence(conv, out):
    if not conv:
        return
    fig, ax = plt.subplots(1, 2, figsize=(WDOUBLE * 0.7, 50 * MM), constrained_layout=True)
    for k, v in conv.items():
        N = [r['N'] for r in v]
        ax[0].plot(N, [r['P00'] for r in v], 'o-', color=COL[k[0]], ms=3, label='%s, %s' % (k[0], k[1]))
        ax[1].plot(N, [r['W00'] for r in v], 'o-', color=COL[k[0]], ms=3)
    ax[0].set_ylabel(r'$P_{00}$')
    ax[1].set_ylabel(r'$W_{00}/\nu_r$')
    for a in ax:
        a.set_xlabel('$N$ per axis')
    ax[0].legend(fontsize=5)
    fig.savefig(out / 'convergence.pdf')
    plt.close(fig)


def traces(runs, out):
    fig, ax = plt.subplots(figsize=(WDOUBLE * 0.5, 50 * MM), constrained_layout=True)
    for sch, c in COL.items():
        b = best([r for r in runs if r['scheme'] == sch])
        if b:
            ax.plot(b['t'], b['P00_t'], color=c, lw=0.8, label=LAB[sch])
    ax.set_xlabel(r'$t\nu_r$')
    ax.set_ylabel(r'$P_{00}$ from $|1,1\rangle$')
    ax.legend(fontsize=5)
    fig.savefig(out / 'traces.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', default=str(R / 'zeeman2d' / 'runs'))
    a = ap.parse_args()
    runs = load_runs(a.runs)
    out = Path(a.runs).parent
    out.mkdir(parents=True, exist_ok=True)
    if not runs:
        print('no runs in', a.runs)
        return
    conv = summary(runs, out)
    overview(runs, out)
    detunings(runs, out)
    convergence(conv, out)
    traces(runs, out)
    print('figures and summary in', out)


if __name__ == '__main__':
    main()
