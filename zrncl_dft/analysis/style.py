"""Shared matplotlib style: validated categorical palette, recessive axes."""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

# categorical slots in fixed order (validated palette, light surface)
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
MAGENTA, GREEN, VIOLET, RED = '#e87ba4', '#008300', '#4a3aa7', '#e34948'
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
SURFACE = '#fcfcfb'

def fermi_level(tag, npz=None):
    """Dense-mesh Fermi level for a run: GPAW's fixed_density() keeps the SCF
    Fermi level, so the doped/intercalate analyses store a tetrahedron
    electron-counting value in their JSON; fall back to the npz 'ef'."""
    import json
    import os
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
    key = {'intercalate': ('intercalate_electronic.json', 'ef_dos')}.get(
        tag, (f'{tag}_electronic.json', 'ef_tetra'))
    try:
        v = json.load(open(os.path.join(res, key[0]))).get(key[1])
        if v is not None:
            return float(v)
    except FileNotFoundError:
        pass
    return float(npz['ef']) if npz is not None else None


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
