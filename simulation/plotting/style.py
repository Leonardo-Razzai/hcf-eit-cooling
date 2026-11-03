"""Nature-journal figure style shared by all plotting scripts.

Arial 7 pt text, 8 pt bold lowercase panel letters, double-column width
183 mm, outward ticks, no top/right spines, Okabe-Ito colour-blind-safe
palette, DejaVu Sans for math (Arial's italic nu looks like a v).
Rule: legends and panel letters go outside the axes, or inside regions the
data provably leave empty (fix the axis limits), never over the curves.
"""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

MM = 1 / 25.4
WDOUBLE = 183 * MM            # Nature double-column width
WSINGLE = 89 * MM             # Nature single-column width

# Okabe-Ito palette: case A, case B, grey (simple model), light blue, orange
CA, CB, CG, CL, CO = '#0072B2', '#D55E00', '#595959', '#56B4E9', '#E69F00'

RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'mathtext.fontset': 'dejavusans',
    'font.size': 7, 'axes.labelsize': 7, 'legend.fontsize': 6.5,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'axes.linewidth': 0.6, 'lines.linewidth': 1.0, 'lines.markersize': 3,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'legend.handlelength': 2.0,
    'pdf.fonttype': 42, 'savefig.dpi': 600,
}


def apply():
    plt.rcParams.update(RC)


def letter(ax, s, x=-0.2, y=1.02):
    """Bold lowercase panel letter outside the top-left corner of ax."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight='bold',
            va='bottom', ha='left')
