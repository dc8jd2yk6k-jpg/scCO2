"""Collect all results into results/summary.md (Markdown tables) and
results/summary.json.  Missing inputs are skipped.

Usage: python analysis/summary.py
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


def load(name):
    p = os.path.join(RES, name)
    return json.load(open(p)) if os.path.exists(p) else None


def f(v, fmt='.3f'):
    return '—' if v is None else format(v, fmt)


out, js = [], {}
rel = load('pristine_relax.json')
el = load('pristine_electronic.json')
rb = load('rigid_band.json')
ser = load('doped_series.json')
x10 = load('x0.10_electronic.json')
mol = load('molecule.json')
ie = load('intercalate_electronic.json')
ig = load('intercalate_geometry.json')

if rel:
    out.append('### Pristine β-ZrNCl: structure (R-3m, 6c sites (0,0,z))\n')
    out.append('| | a (Å) | c (Å) | z(Zr) | z(N) | z(Cl) | Zr–N ×3 (Å) | Zr–N apical (Å) | Zr–Cl ×3 (Å) |')
    out.append('|---|---|---|---|---|---|---|---|---|')
    out.append('| starting model (lit.) | 3.6046 | 27.672 | 0.1196 | 0.1981 | 0.3888 | 2.126 | 2.172 | 2.735 |')
    for tag, lab in [('pbe', 'PBE, exp. cell'), ('pbed3', 'PBE+D3(BJ), exp. cell'),
                     ('pbed3_cell', 'PBE+D3(BJ), cell relaxed (800 eV)')]:
        r = rel.get('runs', {}).get(tag)
        if not r:
            continue
        g = r['geometry']
        out.append(f"| {lab} | {r['a']:.4f} | {r['c']:.3f} | {r['z']['Zr']:.4f} | "
                   f"{r['z']['N']:.4f} | {r['z']['Cl']:.4f} | {g['Zr-N'][0]:.3f} | "
                   f"{g['Zr-N'][3]:.3f} | {g['Zr-Cl'][0]:.3f} |")
    out.append('')

if el:
    out.append('### Pristine β-ZrNCl: electronic structure (PBE)\n')
    m = el['m_eff_CB_at_K']
    kz = el['kz_dispersion']
    out.append('| quantity | value |')
    out.append('|---|---|')
    out.append(f"| band gap (indirect, dense mesh) | {el['gap_dense']:.3f} eV |")
    if 'direct_gap_min' in el:
        out.append(f"| smallest direct gap | {el['direct_gap_min']:.3f} eV |")
    out.append(f"| CB curvature mass at K, k_z = 0 (K→Γ / K→M) | {m['K-G']:.3f} / {m['K-M']:.3f} m_e |")
    if rb:
        out.append(f"| CB density-of-states mass (k_z-averaged) | {rb['m_dos']:.3f} m_e |")
    out.append(f"| CB kz dispersion at K (E(K,Z) − E(K,0)) | "
               f"{1e3 * (kz['CB_K(kz=Z)'] - kz['CB_K(kz=0)']):.1f} meV |")
    out.append(f"| VB kz dispersion at Γ (E(Z) − E(Γ)) | "
               f"{1e3 * (kz['VB_Z'] - kz['VB_G(kz=0)']):.1f} meV |")
    out.append('')
    js['pristine'] = el

if rb:
    out.append('### Electron doping: strict rigid band vs extra electrons + jellium (SCF)\n')
    out.append('| x (e/ZrNCl) | E_F−E_CBM rigid (meV) | E_F−E_CBM SCF (meV) | '
               'N(E_F) rigid | N(E_F) SCF | k_F 2D (1/Å) |')
    out.append('|---|---|---|---|---|---|')
    jel = {round(r['x'], 3): r for r in (ser or {}).get('frozen', [])}
    for r in rb['rigid_band']:
        j = jel.get(round(r['x'], 3))
        out.append(f"| {r['x']:.3f} | {1e3 * r['ef_minus_cbm']:.0f} | "
                   f"{f(1e3 * j['ef_minus_cbm'], '.0f') if j else '—'} | "
                   f"{r['N_EF_per_ZrNCl']:.3f} | {f(j['N_EF_per_ZrNCl']) if j else '—'} | "
                   f"{f(j['kF_2D_parabolic']) if j else '—'} |")
    out.append('\nN(E_F) in states/eV/ZrNCl (both spins).  CB-bottom DOS '
               f"{rb['N_CB_bottom_per_ZrNCl']:.3f} states/eV/ZrNCl ⇒ m*_DOS = {rb['m_dos']:.2f} m_e "
               f"(k_z = 0 in-plane curvature: {rb['m_eff_avg']:.2f} m_e).\n")

if mol:
    out.append('### Isolated Co(Cp)₂ and Co(Cp)₂⁺ (PBE, FD, open boundaries)\n')
    out.append('| | Co–C (Å) | Co–Cp centroid (Å) | C–C (Å) |')
    out.append('|---|---|---|---|')
    for tag, lab in [('cocp2_neutral', 'Co(Cp)₂ (S=½)'), ('cocp2_cation', 'Co(Cp)₂⁺ (S=0)')]:
        g = mol[tag]['geometry']
        out.append(f"| {lab} | {g['Co-C_mean']:.3f} | {g['Co-centroid']:.3f} | {g['C-C_mean']:.3f} |")
    if ig:
        g = ig['CoCp2']
        out.append(f"| in ZrNCl{{Co(Cp)₂}}₀.₁₀ | {g['Co-C_mean']:.3f} | {g['Co-centroid']:.3f} | "
                   f"{g['C-C_mean']:.3f} |")
    out.append(f"\nIonisation energy (ΔSCF): adiabatic {mol['IE_adiabatic']:.2f} eV, "
               f"vertical {mol['IE_vertical']:.2f} eV.\n")

if ie:
    out.append('### ZrNCl{Co(Cp)₂}₀.₁₀ model: electronic structure and charge transfer\n')
    out.append('| quantity | value |')
    out.append('|---|---|')
    for key, lab, fmt, scale in [
            ('ef_minus_cbm_layer', 'E_F − E_CBM(ZrNCl layer)', '.0f', 1e3),
            ('layer_gap', 'ZrNCl layer gap (VBM→CBM)', '.3f', 1),
            ('n_el_in_layer_CB_states', 'electrons in ZrNCl conduction band (per Co(Cp)₂)', '.3f', 1),
            ('x_eff_e_per_ZrNCl', 'effective doping x (e/ZrNCl)', '.3f', 1),
            ('N_EF_per_ZrNCl', 'N(E_F) (states/eV/ZrNCl)', '.3f', 1),
            ('magmom_total', 'total moment (spin-polarised start 1 μB on Co)', '.3f', 1),
            ('E_spin_minus_nonspin', 'E(spin) − E(non-spin) (eV)', '.4f', 1)]:
        if key in ie:
            out.append(f'| {lab} | {format(ie[key] * scale, fmt)}'
                       + (' meV |' if key == 'ef_minus_cbm_layer' else ' |'))
    for key, lab in [('bader', 'Bader'), ('hirshfeld', 'Hirshfeld')]:
        if key in ie:
            out.append(f"| {lab} charge of Co(Cp)₂ | {ie[key]['molecule']:+.3f} e |")
    if 'molecular_levels_occupied_rel_EF' in ie:
        occ = ie['molecular_levels_occupied_rel_EF'][-4:]
        emp = ie['molecular_levels_empty_rel_EF'][:4]
        out.append(f"| highest occupied guest levels (E−E_F, eV) | {', '.join(f'{v:.2f}' for v in occ)} |")
        out.append(f"| lowest empty guest levels (E−E_F, eV) | {', '.join(f'{v:.2f}' for v in emp)} |")
    out.append('')

if ie:
    checks = [('PBE', ie.get('Ucheck_U0')),
              ('PBE+U (U_eff = 4 eV, Co 3d)', ie.get('Ucheck_U4')),
              ('PBE, guest levels +0.5 eV', ie.get('shiftcheck_shift0.5'))]
    checks = [(lab, c) for lab, c in checks if c]
    if checks:
        out.append('### Charge-split checks (6×6×1 k, Fermi–Dirac 0.02 eV, same geometry)\n')
        out.append('| | ' + ' | '.join(lab for lab, _ in checks) + ' |')
        out.append('|---' * (len(checks) + 1) + '|')

        def e1(c):
            g = c.get('guest_e1_band_rel_EF')
            return '—' if not g else f"{g['min']:+.2f} … {g['max']:+.2f}"

        rows = [('electrons in ZrNCl CB states, per guest',
                 lambda c: f"{c['n_el_in_layer_CB_states']:.2f}"),
                ('electrons in guest states near E_F',
                 lambda c: f"{c['n_el_in_guest_states_near_EF']:.2f}"),
                ('x_eff (e⁻/ZrNCl)', lambda c: f"{c['x_eff_e_per_ZrNCl']:.3f}"),
                ('E_F − E_CBM(layer) (meV)', lambda c: f"{1e3 * c['ef_minus_cbm_layer']:.0f}"),
                ("guest states within 1 eV of E_F (eV)", e1),
                ('top of the occupied Co 3d (a₁′/e₂′) levels (eV)',
                 lambda c: f"{max(v for v in c['molecular_levels_occupied_rel_EF'] if v < -1):+.2f}"
                 if any(v < -1 for v in c['molecular_levels_occupied_rel_EF']) else '—'),
                ('layer gap (eV)', lambda c: f"{c['layer_gap']:.2f}"),
                ('N(E_F), layer states (states/eV/ZrNCl)',
                 lambda c: f"{c['N_EF_layer_weighted_per_ZrNCl']:.3f}")]
        for lab, fn in rows:
            out.append(f'| {lab} | ' + ' | '.join(fn(c) for _, c in checks) + ' |')
        out.append('')
    sp = ie.get('spin_check')
    if sp:
        line = (f"Spin-polarised SCF (6×6×1, started from 1 μB on Co): total moment "
                f"{sp['magmom_total']:.3f} μB, Co {sp['magmom_Co']:+.3f} μB, "
                f"Σ|m| {sp['magmom_abs_sum']:.3f} μB")
        if ie.get('Ucheck_U0'):
            line += (f"; E(spin) − E(non-spin) = "
                     f"{1e3 * (sp['energy_spin'] - ie['Ucheck_U0']['energy']):.1f} meV per cell")
        out.append(line + '.\n')

if mol and 'curvature_e1_SOMO_minus_cation_LUMO' in mol:
    out.append(f"Isolated molecule: ε_SOMO(Co(Cp)₂) − ε_LUMO(Co(Cp)₂⁺) at the neutral geometry = "
               f"{mol['curvature_e1_SOMO_minus_cation_LUMO']:.2f} eV (PBE curvature of E(N) for the "
               "e₁″ level; zero for the exact functional).\n")

with open(os.path.join(RES, 'summary.md'), 'w') as fh:
    fh.write('\n'.join(out))
print('\n'.join(out))
