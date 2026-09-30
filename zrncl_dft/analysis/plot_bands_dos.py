"""Band structure + projected DOS figure for one electronic.py run.

Usage: python analysis/plot_bands_dos.py pristine [x0.10 ...]
Energies are referenced to the VBM (pristine) or E_F (doped).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import plt, BLUE, ORANGE, AQUA, INK, INK2, MUTED, GRID  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)

LABEL = {'G': 'Γ', 'M': 'M', 'K': 'K', 'Z': 'Z'}


def smooth(y, e, width=0.04):
    """Gaussian smoothing of a tetrahedron DOS for display only."""
    de = e[1] - e[0]
    n = int(4 * width / de)
    k = np.exp(-0.5 * (np.arange(-n, n + 1) * de / width) ** 2)
    return np.convolve(y, k / k.sum(), mode='same')


def plot(tag, emin=-7.0, emax=4.5):
    d = np.load(os.path.join(RES, f'{tag}_electronic.npz'), allow_pickle=True)
    nval = 24
    doped = not tag.startswith('pristine')
    if doped:
        e0 = float(d['ef'])
        zero_label = 'E − E$_F$ (eV)'
    else:
        e0 = float(d['eig_dense'][:, nval - 1].max())
        zero_label = 'E − E$_{VBM}$ (eV)'
    x, eps = d['x'], d['eps_band'] - e0
    fig, (ax, axd) = plt.subplots(1, 2, figsize=(7.2, 4.2), sharey=True,
                                  gridspec_kw={'width_ratios': [2.3, 1], 'wspace': 0.05})
    for n in range(eps.shape[1]):
        col = BLUE if n >= nval else INK2
        ax.plot(x, eps[:, n], color=col, lw=1.3 if n >= nval else 1.0)
    for t in d['ticks']:
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)
    ax.axhline(0, color=MUTED, lw=0.8, ls='--')
    ax.set_xticks(d['ticks'])
    ax.set_xticklabels([LABEL[s] for s in d['labels']])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(emin, emax)
    ax.set_ylabel(zero_label)

    e = d['energies'] - e0
    keys = list(d['pdos_keys'])
    pd = dict(zip(keys, d['pdos']))
    tot = smooth(d['dos'], e) / 2          # per ZrNCl formula unit
    axd.fill_betweenx(e, 0, tot, color=GRID, lw=0, label='total')
    for key, col, lab in [('Zr-d', BLUE, 'Zr 4d'), ('N-p', ORANGE, 'N 2p'),
                          ('Cl-p', AQUA, 'Cl 3p')]:
        axd.plot(smooth(pd[key], e) / 2, e, color=col, lw=1.4, label=lab)
    axd.axhline(0, color=MUTED, lw=0.8, ls='--')
    axd.set_xlim(0, None)
    axd.set_xlabel('DOS (states/eV/ZrNCl)')
    axd.legend(loc='upper right')
    title = 'β-ZrNCl, PBE'
    if doped:
        title = f'β-ZrNCl + {tag.replace("x", "")} e/ZrNCl (jellium), PBE'
    fig.suptitle(title, x=0.12, ha='left', fontsize=11, color=INK)
    out = os.path.join(FIG, f'{tag}_bands_dos.png')
    fig.savefig(out, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    for t in sys.argv[1:] or ['pristine']:
        plot(t)
