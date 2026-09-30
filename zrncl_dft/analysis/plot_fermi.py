"""Fermi-surface cross-sections (kz-projected) of electron-doped beta-ZrNCl.

Uses the dense-mesh conduction-band energies stored by electronic.py: the
self-consistent jellium x = 0.10 run and the pristine run with the rigid-band
Fermi level for x = 0.10.  The IBZ (time-reversal only) is unfolded to the
full mesh; the first kz plane of the compact-cell mesh is contoured (the CB
disperses by only a few meV along kz).

Usage: python analysis/plot_fermi.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import plt, BLUE, ORANGE, INK2, MUTED  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(HERE, '..', 'figures')
A = 3.6046


def full_grid(npz, band=24):
    d = np.load(os.path.join(RES, npz), allow_pickle=True)
    ibzk, eig = d['ibzk'], d['eig_dense'][:, band]
    size = tuple(int(round(1 / np.min(np.abs(ibzk[:, i][np.abs(ibzk[:, i]) > 1e-6]))))
                 for i in range(3))
    g = np.full(size, np.nan)
    for k, e in zip(ibzk, eig):
        for s in (1, -1):
            idx = tuple(int(v) for v in np.rint(np.mod(s * k * size, size)).astype(int) % size)
            g[idx] = e
    return g, size, d


def tile_plane(g, size):
    """kz-plane l=0 on in-plane scaled coords, tiled 3x3 for contouring."""
    n1, n2 = size[0], size[1]
    plane = g[:, :, 0]
    i = np.arange(-n1, 2 * n1 + 1)
    j = np.arange(-n2, 2 * n2 + 1)
    I, J = np.meshgrid(i, j, indexing='ij')
    E = plane[I % n1, J % n2]
    b1 = 2 * np.pi / A * np.array([1, 1 / np.sqrt(3)])
    b2 = 2 * np.pi / A * np.array([0, 2 / np.sqrt(3)])
    KX = (I / n1)[..., None] * b1 + (J / n2)[..., None] * b2
    return KX[..., 0], KX[..., 1], E


def main():
    rb = json.load(open(os.path.join(RES, 'rigid_band.json')))
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    kK = 4 * np.pi / (3 * A)
    hexagon = np.array([[kK * np.cos(np.radians(60 * t)), kK * np.sin(np.radians(60 * t))]
                        for t in range(7)])
    ax.plot(hexagon[:, 0], hexagon[:, 1], color=MUTED, lw=1.0)
    g, size, d = full_grid('pristine_electronic.npz')
    x10 = [r for r in rb['rigid_band'] if abs(r['x'] - 0.10) < 1e-9][0]
    ef_rb = rb['cbm'] + x10['ef_minus_cbm']
    X, Y, E = tile_plane(g, size)
    ax.contour(X, Y, E, levels=[ef_rb], colors=[BLUE], linewidths=1.6)
    lines = [plt.Line2D([], [], color=BLUE, lw=1.6, label='rigid band, x = 0.10')]
    if os.path.exists(os.path.join(RES, 'x0.10_electronic.npz')):
        g2, size2, d2 = full_grid('x0.10_electronic.npz')
        X2, Y2, E2 = tile_plane(g2, size2)
        ax.contour(X2, Y2, E2, levels=[float(d2['ef'])], colors=[ORANGE],
                   linewidths=1.6, linestyles='--')
        lines.append(plt.Line2D([], [], color=ORANGE, lw=1.6, ls='--',
                                label='extra e⁻ + jellium, x = 0.10'))
    for t in range(6):
        ang = np.radians(60 * t)      # K = (b1+b2)/3 sits at 60 deg
        ax.text(1.12 * kK * np.cos(ang), 1.12 * kK * np.sin(ang), "K" if t % 2 == 1 else "K′",
                ha='center', va='center', fontsize=9, color=INK2)
    ax.text(0, 0, 'Γ', ha='center', va='center', fontsize=10, color=INK2)
    ax.set_aspect('equal')
    lim = 1.3 * kK
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel('k$_x$ (Å$^{-1}$)')
    ax.set_ylabel('k$_y$ (Å$^{-1}$)')
    ax.legend(handles=lines, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=1)
    fig.savefig(os.path.join(FIG, 'fermi_surface_x0.10.png'), bbox_inches='tight')
    print('wrote fermi surface figure')


if __name__ == '__main__':
    main()
