"""Strict rigid-band doping of beta-ZrNCl from the pristine DOS.

The pristine conduction-band DOS (linear tetrahedron, dense k mesh) is
integrated from the CBM:  n(E) = int_CBM^E DOS dE'  (electrons per Zr2N2Cl2),
x = n/2 electrons per ZrNCl.  For each x the rigid-band Fermi level and
N(E_F) are compared with the self-consistent extra-electron (jellium) runs of
02_doped/series.py and, if present, with the explicit Co(Cp)2 intercalate.

Usage: python analysis/rigid_band.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import plt, BLUE, ORANGE, AQUA, INK, INK2, MUTED, GRID  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
HB2M = 3.80998  # hbar^2/2m_e, eV A^2


def rigid_band(tag='pristine'):
    d = np.load(os.path.join(RES, f'{tag}_electronic.npz'), allow_pickle=True)
    e, dos = d['energies'], d['dos']            # states/eV/cell (both spins)
    cbm = float(d['eig_dense'][:, 24].min())
    sel = e >= cbm
    ee, dd = e[sel], dos[sel]
    n = np.concatenate([[0], np.cumsum(0.5 * (dd[1:] + dd[:-1]) * np.diff(ee))])
    x = n / 2
    return cbm, ee, dd, x


def main():
    cbm, ee, dd, xx = rigid_band()
    js = json.load(open(os.path.join(RES, 'pristine_electronic.json')))
    area = 3.6046 ** 2 * np.sqrt(3) / 2
    mstar = js['m_eff_CB_at_K_avg']
    # 2D parabolic band, valleys K and K', both spins, per ZrNCl (2 per cell):
    # N = g_s g_v m* A_cell / (2 pi hbar^2) / 2 = (m*/m_e) A_cell / (2 pi HB2M)
    n2d = mstar * area / (2 * np.pi * HB2M)
    targets = [0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    rb = []
    for x in targets:
        ef = float(np.interp(x, xx, ee))
        nef = float(np.interp(ef, ee, dd)) / 2
        rb.append({'x': x, 'ef_minus_cbm': ef - cbm, 'N_EF_per_ZrNCl': nef})
    # DOS mass: mean DOS in [CBM+0.05, CBM+0.30] eV (single K/K' band, both spins)
    win = (ee > cbm + 0.05) & (ee < cbm + 0.30)
    ndos = float(dd[win].mean()) / 2                       # per ZrNCl
    m_dos = ndos * 2 * np.pi * HB2M / area
    out = {'cbm': cbm, 'm_eff_avg': mstar, 'N2D_parabolic_per_ZrNCl': n2d,
           'N_CB_bottom_per_ZrNCl': ndos, 'm_dos': m_dos,
           'note': 'm_eff_avg is the kz=0 in-plane curvature at K; interlayer hopping '
                   'vanishes at K but grows ~|k-K| away from it, so the kz-averaged '
                   '(DOS) mass m_dos is the relevant one for N(E_F).',
           'rigid_band': rb}

    jel = None
    p = os.path.join(RES, 'doped_series.json')
    if os.path.exists(p):
        jel = [r for r in json.load(open(p))['frozen'] if r['x'] > 0]
        out['jellium'] = [{k: r[k] for k in ('x', 'ef_minus_cbm', 'N_EF_per_ZrNCl')}
                          for r in jel]
    dense = None
    p = os.path.join(RES, 'x0.10_electronic.json')
    if os.path.exists(p):
        dense = json.load(open(p))
        out['jellium_x0.10_dense'] = {k: dense.get(k) for k in (
            'ef_tetra_minus_cbm', 'ef_dos_minus_cbm', 'N_EF_tetra_per_ZrNCl',
            'gamma_valley_minus_cbm', 'kdos')}
    inter = None
    p = os.path.join(RES, 'intercalate_electronic.json')
    if os.path.exists(p):
        inter = json.load(open(p))

    with open(os.path.join(RES, 'rigid_band.json'), 'w') as f:
        json.dump(out, f, indent=2)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.6, 3.3))
    xs = np.linspace(0.01, 0.5, 400)
    ef_abs = np.interp(xs, xx, ee)
    a1.plot(xs, 1e3 * (ef_abs - cbm), color=BLUE, label='rigid band (pristine DOS)')
    a2.plot(xs, np.interp(ef_abs, ee, dd) / 2, color=BLUE,
            label='rigid band (pristine DOS)')
    if jel:
        lab = 'extra e⁻ + jellium (SCF mesh)'
        a1.plot([r['x'] for r in jel], [1e3 * r['ef_minus_cbm'] for r in jel], 'o',
                ms=6, mfc=ORANGE, mec='white', mew=1.5, label=lab)
        a2.plot([r['x'] for r in jel], [r['N_EF_per_ZrNCl'] for r in jel], 'o',
                ms=6, mfc=ORANGE, mec='white', mew=1.5, label=lab)
    if dense and 'ef_tetra_minus_cbm' in dense:
        lab = 'extra e⁻ + jellium, dense mesh'
        a1.plot([0.10], [1e3 * dense['ef_tetra_minus_cbm']], 's', ms=7, mfc=ORANGE,
                mec=INK, mew=1.0, label=lab)
        a2.plot([0.10], [dense['N_EF_tetra_per_ZrNCl']], 's', ms=7, mfc=ORANGE,
                mec=INK, mew=1.0, label=lab)
    if inter and 'ef_minus_cbm_layer' in inter:
        a1.plot([0.10], [1e3 * inter['ef_minus_cbm_layer']], 'D', ms=7, mfc=AQUA,
                mec='white', mew=1.5, label='ZrNCl{Co(Cp)₂}₀.₁₀')
        if 'N_EF_per_ZrNCl' in inter:
            a2.plot([0.10], [inter['N_EF_per_ZrNCl']], 'D', ms=7, mfc=AQUA,
                    mec='white', mew=1.5, label='ZrNCl{Co(Cp)₂}₀.₁₀')
    for a in (a1, a2):
        a.axvline(0.10, color=GRID, lw=5, zorder=0)   # x of ZrNCl{Co(Cp)2}0.10
        a.set_xlabel('x (electrons per ZrNCl)')
        a.set_xlim(0, 0.5)
        a.grid(True, axis='y')
    a1.set_ylabel('E$_F$ − E$_{CBM}$ (meV)')
    a2.set_ylabel('N(E$_F$) (states/eV/ZrNCl)')
    a2.set_ylim(0, None)
    a1.legend(loc='upper left')
    a2.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'doping_rigid_vs_scf.png'), bbox_inches='tight')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
