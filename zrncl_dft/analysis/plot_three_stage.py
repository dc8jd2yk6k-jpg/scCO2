"""One-figure summary: bands along G-M-K-G for
  (a) pristine beta-ZrNCl  (rigid-band E_F for x = 0.10 marked),
  (b) beta-ZrNCl + 0.10 e-/ZrNCl on a jellium background (self-consistent),
  (c) ZrNCl{Co(Cp)2}0.10, supercell bands unfolded onto the 1x1 zone.
All panels are aligned at the conduction-band minimum (K) of the ZrNCl layer.

Usage: python analysis/plot_three_stage.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import plt, fermi_level, BLUE, ORANGE, INK, INK2, MUTED, GRID  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(HERE, '..', 'figures')
EMIN, EMAX = -2.6, 1.6


def host_panel(ax, tag, ef_rel, title, ef_label):
    d = np.load(os.path.join(RES, f'{tag}_electronic.npz'), allow_pickle=True)
    x, eps = d['x'], d['eps_band']
    keep = x <= d['ticks'][3] + 1e-9
    cbm = eps[keep, 24].min()
    for n in range(eps.shape[1]):
        ax.plot(x[keep], eps[keep, n] - cbm, color=BLUE if n >= 24 else INK2,
                lw=1.3 if n >= 24 else 1.0)
    ax.axhline(ef_rel, color=INK, lw=1.0, ls='--')
    ax.text(x[keep][-1] * 0.02, ef_rel + 0.05, ef_label, ha='left', va='bottom',
            fontsize=8, color=INK)
    ticks = d['ticks'][:4]
    for t in ticks:
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(['Γ', 'M', 'K', 'Γ'])
    ax.set_xlim(x[0], ticks[-1])
    ax.set_title(title, loc='left', fontsize=9.5, color=INK2)
    return x[keep][-1]


def main():
    rb = json.load(open(os.path.join(RES, 'rigid_band.json')))
    x10 = [r for r in rb['rigid_band'] if abs(r['x'] - 0.10) < 1e-9][0]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.2), sharey=True,
                             gridspec_kw={'wspace': 0.06})
    host_panel(axes[0], 'pristine', x10['ef_minus_cbm'], '(a) β-ZrNCl, rigid band',
               f"E$_F$(x=0.10) = CBM + {1e3 * x10['ef_minus_cbm']:.0f} meV")
    if os.path.exists(os.path.join(RES, 'x0.10_electronic.npz')):
        d = np.load(os.path.join(RES, 'x0.10_electronic.npz'), allow_pickle=True)
        keep = d['x'] <= d['ticks'][3] + 1e-9
        ef_rel = fermi_level('x0.10', d) - d['eps_band'][keep, 24].min()
        host_panel(axes[1], 'x0.10', ef_rel, '(b) + 0.10 e⁻/ZrNCl, jellium (SCF)',
                   f'E$_F$ = CBM + {1e3 * ef_rel:.0f} meV')
    f = os.path.join(RES, 'intercalate_unfold.npz')
    ax = axes[2]
    if os.path.exists(f):
        d = np.load(f, allow_pickle=True)
        x, e, P, F = d['x'], d['e_kn'], d['P_kn'], d['Fmol_kn']
        lay = (F < 0.5) & (P > 0.05)
        ef_i = fermi_level('intercalate', d)
        cbm = e[lay & (e > ef_i - 1.0)].min()
        e = e - cbm
        ef_rel = ef_i - cbm
        X = np.repeat(x[:, None], e.shape[1], axis=1)
        s = 26 * P
        ax.scatter(X[lay], e[lay], s=s[lay], color=BLUE, lw=0, zorder=2)
        mol = F >= 0.5
        ax.scatter(X[mol], e[mol], s=6, color=ORANGE, lw=0, zorder=3, alpha=0.9)
        ax.axhline(ef_rel, color=INK, lw=1.0, ls='--')
        ax.text(x[-1] * 0.02, ef_rel - 0.08, f'E$_F$ = CBM + {1e3 * ef_rel:.0f} meV',
                ha='left', va='top', fontsize=8, color=INK)
        for t in d['ticks']:
            ax.axvline(t, color=GRID, lw=0.8, zorder=0)
        ax.set_xticks(d['ticks'])
        ax.set_xticklabels(['Γ', 'M', 'K', 'Γ'])
        ax.set_xlim(x[0], x[-1])
        ax.scatter([], [], s=20, color=BLUE, label='ZrNCl layer (unfolded weight)')
        ax.scatter([], [], s=20, color=ORANGE, label='Co(Cp)₂ guest states')
        ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.45), fontsize=7.5)
    ax.set_title('(c) ZrNCl{Co(Cp)₂}₀.₁₀, unfolded', loc='left', fontsize=9.5, color=INK2)
    axes[0].set_ylabel('E − E$_{CBM}$(ZrNCl) (eV)')
    axes[0].set_ylim(EMIN, EMAX)
    for a in axes:
        a.axhline(0, color=MUTED, lw=0.6, ls=':')
    fig.savefig(os.path.join(FIG, 'three_stage_bands.png'), bbox_inches='tight')
    print('wrote figures/three_stage_bands.png')


if __name__ == '__main__':
    main()
