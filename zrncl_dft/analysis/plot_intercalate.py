"""Figures for the ZrNCl{Co(Cp)2}0.10 model: PDOS, unfolded bands and the
plane-averaged density difference (structure panels: plot_structures.py).
Each panel is skipped if its data file is missing.

Usage: python analysis/plot_intercalate.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import (plt, BLUE, ORANGE, AQUA, YELLOW, MAGENTA, INK, INK2,  # noqa: E402
                   MUTED, GRID)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
RES = os.path.join(ROOT, 'results')
FIG = os.path.join(ROOT, 'figures')
RUN = os.path.join(ROOT, 'runs', 'intercalate')
os.makedirs(FIG, exist_ok=True)


def smooth(y, e, width=0.05):
    de = e[1] - e[0]
    n = int(4 * width / de)
    k = np.exp(-0.5 * (np.arange(-n, n + 1) * de / width) ** 2)
    return np.convolve(y, k / k.sum(), mode='same')


def pdos():
    f = os.path.join(RES, 'intercalate_dos.npz')
    if not os.path.exists(f):
        return
    d = np.load(f, allow_pickle=True)
    e = d['energies'] - d['ef']
    pd = dict(zip(list(d['keys']), d['pdos']))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.2), sharex=True,
                                 gridspec_kw={'hspace': 0.08})
    a1.fill_between(e, 0, smooth(d['dos'], e), color=GRID, lw=0, label='total')
    for key, col, lab in [('Zr-d', BLUE, 'Zr 4d'), ('N-p', ORANGE, 'N 2p'),
                          ('Cl-p', AQUA, 'Cl 3p')]:
        a1.plot(e, smooth(pd[key], e), color=col, lw=1.3, label=lab)
    a1.set_ylabel('DOS (states/eV/cell)')
    a1.set_title('ZrNCl layer', loc='left', fontsize=10, color=INK2)
    a1.legend(loc='upper left', ncol=4)
    for key, col, lab in [('Co-d', BLUE, 'Co 3d'), ('C-p', ORANGE, 'C 2p')]:
        a2.plot(e, smooth(pd[key], e, 0.03), color=col, lw=1.3, label=lab)
    a2.set_ylabel('PDOS (states/eV)')
    a2.set_title('Co(Cp)₂ guest', loc='left', fontsize=10, color=INK2)
    a2.legend(loc='upper left')
    for a in (a1, a2):
        a.axvline(0, color=MUTED, lw=0.8, ls='--')
        a.set_xlim(-6, 3)
        a.set_ylim(0, None)
    a2.set_xlabel('E − E$_F$ (eV)')
    fig.savefig(os.path.join(FIG, 'intercalate_pdos.png'), bbox_inches='tight')
    # zoom around the gap / E_F
    fig, a = plt.subplots(figsize=(6.4, 2.8))
    a.fill_between(e, 0, smooth(d['dos'], e, 0.03), color=GRID, lw=0, label='total')
    a.plot(e, smooth(pd['Zr-d'], e, 0.03), color=BLUE, lw=1.3, label='Zr 4d (layer)')
    a.plot(e, smooth(pd['Co-d'], e, 0.03), color=ORANGE, lw=1.3, label='Co 3d (guest)')
    a.axvline(0, color=MUTED, lw=0.8, ls='--')
    a.set_xlim(-2.5, 2.0)
    a.set_ylim(0, None)
    a.set_xlabel('E − E$_F$ (eV)')
    a.set_ylabel('DOS (states/eV/cell)')
    a.legend(loc='upper left')
    fig.savefig(os.path.join(FIG, 'intercalate_pdos_zoom.png'), bbox_inches='tight')
    print('wrote pdos')


def unfold():
    f = os.path.join(RES, 'intercalate_unfold.npz')
    if not os.path.exists(f):
        return
    d = np.load(f, allow_pickle=True)
    x, e, P, F = d['x'], d['e_kn'] - d['ef'], d['P_kn'], d['Fmol_kn']
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ref = os.path.join(RES, 'x0.10_electronic.npz')
    if os.path.exists(ref):
        r = np.load(ref, allow_pickle=True)
        xr, er = r['x'], r['eps_band'] - r['ef']
        keep = xr <= r['ticks'][3] + 1e-9                       # G-M-K-G part
        for n in range(er.shape[1]):
            ax.plot(xr[keep], er[keep, n], color=MUTED, lw=0.9, zorder=1,
                    label='β-ZrNCl + 0.10 e⁻ (jellium)' if n == 0 else None)
    X = np.repeat(x[:, None], e.shape[1], axis=1)
    lay = F < 0.5
    s = 28 * P
    ax.scatter(X[lay], e[lay], s=s[lay], color=BLUE, lw=0, zorder=2,
               label='intercalate, layer states')
    ax.scatter(X[~lay], e[~lay], s=np.maximum(s[~lay], 4), color=ORANGE, lw=0, zorder=3,
               label='intercalate, Co(Cp)₂ states')
    for t in d['ticks']:
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)
    ax.axhline(0, color=MUTED, lw=0.8, ls='--')
    ax.set_xticks(d['ticks'])
    ax.set_xticklabels(['Γ', 'M', 'K', 'Γ'])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-4, 3)
    ax.set_ylabel('E − E$_F$ (eV)')
    ax.legend(loc='lower left', fontsize=7.5)
    ax.set_title('Unfolded bands of ZrNCl{Co(Cp)₂}₀.₁₀ (1×1 ZrNCl zone)', loc='left',
                 fontsize=10)
    fig.savefig(os.path.join(FIG, 'intercalate_unfolded_bands.png'), bbox_inches='tight')
    print('wrote unfolded bands')


def drho():
    f = os.path.join(RES, 'intercalate_drho.npz')
    if not os.path.exists(f):
        return
    d = np.load(f, allow_pickle=True)
    z, prof, cdc = d['z'], d['drho_z'], d['cdc']
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 4.6), sharex=True,
                                 gridspec_kw={'hspace': 0.08})
    a1.axhline(0, color=MUTED, lw=0.8)
    a1.fill_between(z, 0, prof, where=prof > 0, color=BLUE, lw=0, alpha=0.8,
                    label='electron gain')
    a1.fill_between(z, 0, prof, where=prof < 0, color=ORANGE, lw=0, alpha=0.8,
                    label='electron loss')
    a1.set_ylabel('Δρ̄(z) (e/Å)')
    a1.legend(loc='upper right')
    a2.plot(z, cdc, color=INK, lw=1.4)
    a2.axhline(0, color=MUTED, lw=0.8)
    a2.set_ylabel('ΔQ(z) (e)')
    a2.set_xlabel('z (Å)')
    zs, sy = d['zatoms'], d['symbols']
    for a in (a1, a2):
        for s, col in [('Zr', BLUE), ('Cl', AQUA), ('Co', YELLOW)]:
            for zz in sorted(set(np.round(zs[sy == s], 1))):
                a.axvline(zz, color=col, lw=0.8, ls=':', zorder=0)
    a1.text(0.01, 0.92, 'dotted: Zr (blue), Cl (aqua), Co (yellow) planes',
            transform=a1.transAxes, fontsize=7.5, color=INK2)
    fig.savefig(os.path.join(FIG, 'intercalate_drho.png'), bbox_inches='tight')
    print('wrote drho')


if __name__ == '__main__':
    pdos()
    unfold()
    drho()
