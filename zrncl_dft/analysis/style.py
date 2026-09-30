"""Shared matplotlib style: validated categorical palette, recessive axes."""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

# categorical slots in fixed order (validated palette, light surface)
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
MAGENTA, GREEN, VIOLET, RED = '#e87ba4', '#008300', '#4a3aa7', '#e34948'
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
SURFACE = '#fcfcfb'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE, 'savefig.dpi': 180,
    'axes.edgecolor': MUTED, 'axes.linewidth': 0.8,
    'axes.labelcolor': INK, 'axes.titlecolor': INK,
    'axes.titlesize': 11, 'axes.labelsize': 10,
    'xtick.color': INK2, 'ytick.color': INK2,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 8.5, 'legend.frameon': False,
    'font.family': 'DejaVu Sans', 'lines.linewidth': 1.6,
    'grid.color': GRID, 'grid.linewidth': 0.6,
})
