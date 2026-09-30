"""Electronic structure and charge transfer of the ZrNCl{Co(Cp)2}0.10 model.

Stages (run in order; each checkpoints to runs/intercalate/):
  scf        non-spin-polarised SCF, 9x9x1 k (contains the folded K/K'),
             Fermi-Dirac 0.02 eV, followed directly by the DOS analysis:
             tetrahedron DOS / PDOS, electron count in the ZrNCl conduction
             band, guest levels, N(E_F)
  spin       spin-polarised SCF (6x6x1) started from 1 mu_B on Co
  dos        (optional) the DOS analysis from a stored scf.gpw
  unfold     non-SCF along G-M-K-G of the 1x1 ZrNCl cell, spectral weights
             (band unfolding of the pseudo wave functions)
  charges    Hirshfeld and Bader (all-electron density) charges per fragment
  fragments  slab-only and molecule-only SCFs at the frozen geometry ->
             plane-averaged density difference and charge-displacement curve

Run:  mpiexec -n 4 python 03_intercalate/electronic.py --stage scf   (etc.)
"""
import argparse
import json
import os
import sys

import numpy as np
from ase.io import read
from gpaw import GPAW
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import gpaw_calc, unfold_path, unfold_weights, RESULTS, ROOT  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument('--stage', required=True,
               choices=['scf', 'spin', 'dos', 'unfold', 'charges', 'fragments', 'scfU'])
p.add_argument('--structure', default=None)
p.add_argument('--U', type=float, default=4.0, help='U_eff on Co 3d (eV) for --stage scfU')
args = p.parse_args()

RUN = os.path.join(ROOT, 'runs', 'intercalate')
NB = 190
NMOL = 21
KSCF = {'size': (9, 9, 1), 'gamma': True}   # contains the folded K/K' (0, +-1/3)
WIDTH = 0.02                                 # eV; 0.05 eV put ~0.2 e into the guest
                                             # e1'' level lying only ~0.15 eV above E_F
KAUX = {'size': (6, 6, 1), 'gamma': True}    # spin-polarised check
CONV = {'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8}
M = np.array([[2, 1, 0], [-1, 2, 0], [0, 0, 1]])   # A_SC = M a_PC


def structure():
    if args.structure:
        at = read(args.structure)
    else:
        for name in ('relaxed_pw.traj', 'relaxed_lcao.traj'):
            f = os.path.join(RUN, name)
            if os.path.exists(f):
                at = read(f)
                break
    at.set_constraint()
    return at


def update_json(new):
    if world.rank != 0:
        return
    path = os.path.join(RESULTS, 'intercalate_electronic.json')
    data = json.load(open(path)) if os.path.exists(path) else {}
    data.update(new)
    json.dump(data, open(path, 'w'), indent=2, default=float)


def layer_mol_indices(atoms):
    n = len(atoms)
    return list(range(n - NMOL)), list(range(n - NMOL, n))


# --------------------------------------------------------------------- scf
if args.stage == 'scf':
    atoms = structure()
    atoms.calc = gpaw_calc(os.path.join(RUN, 'scf.txt'), KSCF, width=WIDTH,
                           nbands=NB, convergence=dict(CONV, bands=NB - 10))
    e = atoms.get_potential_energy()
    atoms.calc.write(os.path.join(RUN, 'scf.gpw'))
    update_json({'energy_nonspin': e, 'ef_scf': atoms.calc.get_fermi_level(),
                 'kscf': KSCF['size'], 'width': WIDTH,
                 'structure_used': args.structure or 'relaxed_pw/lcao'})
    args.stage = 'dos-from-scf'
    dos_calc = atoms.calc      # analysed below by dos_analysis()

# --------------------------------------------------------------------- spin
if args.stage == 'spin':
    atoms = structure()
    mag = np.zeros(len(atoms))
    mag[len(atoms) - NMOL] = 1.0          # Co
    atoms.set_initial_magnetic_moments(mag)
    atoms.calc = gpaw_calc(os.path.join(RUN, 'scf_spin.txt'), KAUX, width=WIDTH,
                           spin=True, nbands=NB, convergence=CONV)
    e = atoms.get_potential_energy()
    mtot = atoms.calc.get_magnetic_moment()
    mloc = atoms.calc.get_magnetic_moments()
    update_json({'spin_check': {'kpts': KAUX['size'], 'width': WIDTH, 'energy_spin': e,
                                'magmom_total': mtot,
                                'magmom_Co': float(mloc[len(atoms) - NMOL]),
                                'magmom_abs_sum': float(np.abs(mloc).sum())}})

# --------------------------------------------------------------------- dos
def dos_analysis(calc, tag='', kmesh=None):
    """Tetrahedron DOS/PDOS, dense-mesh Fermi level by electron counting, and
    the split of the conduction electrons between ZrNCl-layer states and
    Co(Cp)2 guest states.  Saves results/intercalate_dos{tag}.npz, returns dict."""
    atoms = calc.get_atoms()
    ef_scf = calc.get_fermi_level()     # fixed_density keeps the SCF Fermi level
    from gpaw.dos import DOSCalculator
    dc = DOSCalculator.from_calculator(calc, shift_fermi_level=False)
    nk = len(calc.get_ibz_k_points())
    eig = np.array([calc.get_eigenvalues(kpt=k) for k in range(nk)])
    wk = np.asarray(calc.get_k_point_weights())
    nel = calc.get_number_of_electrons()
    # Fermi level of the mesh by electron counting with the tetrahedron DOS:
    # take the highest clean gap below E_F (all bands beneath it are full) and
    # integrate the DOS upward from there.
    ng = None
    for n in range(int(nel // 2), 0, -1):
        top, bot = eig[:, n - 1].max(), eig[:, n].min()
        if bot - top > 0.02 and top < ef_scf - 0.05:
            ng = n
            break
    e_gap = 0.5 * (eig[:, ng - 1].max() + eig[:, ng].min())
    egrid = np.linspace(e_gap, e_gap + 4.0, 4001)
    dgrid = dc.raw_dos(egrid, width=0.0)
    ncum = 2 * ng + np.concatenate([[0], np.cumsum(0.5 * (dgrid[1:] + dgrid[:-1]) * np.diff(egrid))])
    ef = float(np.interp(nel, ncum, egrid))
    energies = np.linspace(ef - 8, ef + 4, 2401)
    dos = dc.raw_dos(energies, width=0.0)
    sym = atoms.get_chemical_symbols()
    lay, mol = layer_mol_indices(atoms)
    comp = {}
    for key, idx, l in [('Zr-d', [i for i in lay if sym[i] == 'Zr'], 2),
                        ('N-p', [i for i in lay if sym[i] == 'N'], 1),
                        ('Cl-p', [i for i in lay if sym[i] == 'Cl'], 1),
                        ('Co-d', [i for i in mol if sym[i] == 'Co'], 2),
                        ('C-p', [i for i in mol if sym[i] == 'C'], 1),
                        ('H-s', [i for i in mol if sym[i] == 'H'], 0)]:
        comp[key] = sum(dc.raw_pdos(energies, a=i, l=l, width=0.0) for i in idx)
    # Co-d resolved: GPAW m order xy, yz, 3z2-r2, zx, x2-y2 (lab frame)
    ico = [i for i in mol if sym[i] == 'Co'][0]
    for m, name in enumerate(['dxy', 'dyz', 'dz2', 'dzx', 'dx2-y2']):
        comp[f'Co-{name}'] = dc.raw_pdos(energies, a=ico, l=2, m=m, width=0.0)

    # occupations at the tetrahedron E_F (0..2 per state; 10 meV Fermi-Dirac
    # only to split states exactly at E_F)
    occ = 2.0 / (np.exp(np.clip((eig - ef) / 0.01, -60, 60)) + 1.0)
    # state-resolved character: sum_i |<p_i|psi_nk>|^2 over all projectors
    # of the molecule's atoms vs. the layer's atoms
    from gpaw.dos import IBZWaveFunctions
    wfs = IBZWaveFunctions(calc)
    wmol = np.zeros_like(eig)
    wlay = np.zeros_like(eig)
    for a in range(len(atoms)):
        ni = calc.wfs.setups[a].ni
        w = wfs.pdos_weights(a, list(range(ni))).sum(2)      # (k, n)
        if a in mol:
            wmol += w
        else:
            wlay += w
    frac_mol = wmol / np.maximum(wmol + wlay, 1e-12)

    # classify states: molecular if >50% of the projection weight on Co(Cp)2
    ismol = frac_mol > 0.5
    # layer VBM: highest layer-type state well below E_F (the ZrNCl gap is ~1.8 eV
    # and E_F lies a few 0.1 eV inside the conduction band, so E_F - 1 eV is in the gap)
    lay_occ_max = eig[~ismol & (eig < ef - 1.0)].max()
    # layer conduction-band minimum: lowest layer-type state above the layer gap
    above = eig[~ismol & (eig > lay_occ_max + 1.0)]
    cbm_layer = float(above.min())
    emid = 0.5 * (lay_occ_max + cbm_layer)
    sel = eig > emid
    n_above_layer = float((wk[:, None] * occ * (1 - frac_mol))[sel].sum())
    n_above_mol = float((wk[:, None] * occ * frac_mol)[sel].sum())
    n_cb_layer_states = float((wk[:, None] * occ)[sel & ~ismol].sum())
    n_guest_states = float((wk[:, None] * occ)[sel & ismol].sum())
    # guest e1'' band: molecular states within +-1 eV of E_F
    g = ismol & (np.abs(eig - ef) < 1.0)
    e1 = {'min': float(eig[g].min() - ef), 'max': float(eig[g].max() - ef)} if g.any() else None
    mol_levels = sorted(set(np.round(eig[ismol & (eig > ef - 4) & (eig < ef + 3)], 2)))
    mol_occ = [float(np.round(v - ef, 3)) for v in mol_levels if v < ef]
    mol_emp = [float(np.round(v - ef, 3)) for v in mol_levels if v >= ef]
    nef_tot = float(dc.raw_dos([ef], width=0.0)[0])
    nef_layer = float(2 * dc.calculate(np.array([ef]), eig, 1 - frac_mol, width=0.0)[0])
    nef_guest = float(2 * dc.calculate(np.array([ef]), eig, frac_mol, width=0.0)[0])
    res = {'ef_dos': ef, 'ef_scf_mesh': ef_scf, 'kdos': list(kmesh or KSCF['size']),
           'clean_gap_band_index': ng, 'clean_gap_energy': float(e_gap),
           'n_el_check_sum_occ': float((wk[:, None] * occ).sum()),
           'layer_vbm': float(lay_occ_max), 'layer_cbm': cbm_layer,
           'layer_gap': float(cbm_layer - lay_occ_max),
           'ef_minus_cbm_layer': float(ef - cbm_layer),
           'n_el_above_gap_mid_layer_weighted': n_above_layer,
           'n_el_above_gap_mid_mol_weighted': n_above_mol,
           'n_el_in_layer_CB_states': n_cb_layer_states,
           'n_el_in_guest_states_near_EF': n_guest_states,
           'x_eff_e_per_ZrNCl': n_cb_layer_states / 10,
           'guest_e1_band_rel_EF': e1,
           'N_EF_per_cell': nef_tot, 'N_EF_per_ZrNCl': nef_tot / 10,
           'N_EF_layer_weighted_per_ZrNCl': nef_layer / 10,
           'N_EF_guest_weighted_per_cell': nef_guest,
           'molecular_levels_occupied_rel_EF': mol_occ[-12:],
           'molecular_levels_empty_rel_EF': mol_emp[:12]}
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, f'intercalate_dos{tag}.npz'), energies=energies,
                 dos=dos, keys=list(comp), pdos=np.array(list(comp.values())),
                 ef=ef, eig=eig, wk=wk, occ=occ, frac_mol=frac_mol,
                 wmol=wmol, wlay=wlay)
        print(json.dumps(res, indent=1), flush=True)
    return res


if args.stage in ('dos', 'dos-from-scf'):
    if args.stage == 'dos':
        dos_calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
            kpts=KSCF, nbands=NB,
            symmetry={'point_group': False, 'time_reversal': True},
            convergence={'bands': NB - 10, 'eigenstates': 1e-5},
            txt=os.path.join(RUN, 'dos.txt'))
    update_json(dos_analysis(dos_calc))

# --------------------------------------------------------------------- +U test
if args.stage == 'scfU':
    # sensitivity of the guest level / charge split to the Co 3d on-site
    # interaction (Dudarev U_eff on Co d); same geometry, 6x6x1 mesh
    atoms = structure()
    tag = f'_U{args.U:g}'
    setups = {'Co': f':d,{args.U:g}'} if args.U > 0 else 'paw'
    atoms.calc = gpaw_calc(os.path.join(RUN, f'scf{tag}.txt'), KAUX, width=WIDTH,
                           nbands=NB, convergence=dict(CONV, bands=NB - 10),
                           setups=setups)
    e = atoms.get_potential_energy()
    r = dos_analysis(atoms.calc, tag=tag, kmesh=KAUX['size'])
    r['energy'] = e
    update_json({f'Ucheck{tag}': r})

# --------------------------------------------------------------------- unfold
if args.stage == 'unfold':
    atoms = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).get_atoms()
    a_pc = np.linalg.solve(M.astype(float), atoms.cell.array)   # 1x1 cell rows
    kpc, x, ticks, Ksc, Gt = unfold_path(M, a_pc, dens=20)   # ~55 k points
    # 168 converged bands cover E_F + 2 eV at every k (the top Zr-d bands of a
    # 180-band run converge very slowly and are not plotted)
    calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
        kpts=Ksc, symmetry='off', nbands=180,
        convergence={'bands': 168, 'eigenstates': 1e-4},   # plotting accuracy
        txt=os.path.join(RUN, 'unfold.txt'))
    lay, mol = layer_mol_indices(atoms)
    e_kn, P_kn, Fmol_kn = unfold_weights(calc, M, Gt, mol_atoms=mol)
    ef = calc.get_fermi_level()
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, 'intercalate_unfold.npz'), x=x,
                 ticks=np.array(ticks), labels=np.array(['G', 'M', 'K', 'G']),
                 e_kn=e_kn, P_kn=P_kn, Fmol_kn=Fmol_kn, ef=ef, kpc=kpc)

# --------------------------------------------------------------------- charges
if args.stage == 'charges':
    calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None)
    atoms = calc.get_atoms()
    lay, mol = layer_mol_indices(atoms)
    res = {}
    try:
        from gpaw.analyse.hirshfeld import HirshfeldPartitioning
        hq = np.asarray(HirshfeldPartitioning(calc).get_charges())
        res['hirshfeld'] = {'per_atom': [float(q) for q in hq],
                            'molecule': float(np.sum(hq[mol])),
                            'layer': float(np.sum(hq[lay]))}
    except Exception as err:          # Hirshfeld helper is real-space oriented
        res['hirshfeld_error'] = repr(err)
    rho = calc.get_all_electron_density(gridrefinement=2)       # e/A^3, incl. core
    if world.rank == 0:
        np.save(os.path.join(RUN, 'rho_ae.npy'), rho)
        from pybader.interface import Bader
        from pybader.io.cube import write as cube_write
        lat = np.ascontiguousarray(atoms.cell.array, dtype=float)
        vol = abs(np.linalg.det(lat))
        info = {'filename': 'rho_ae', 'prefix': RUN, 'file_type': 'gpaw',
                'write_function': cube_write,
                'elements': atoms.get_atomic_numbers(), 'voxel_offset': np.zeros(3)}
        bader = Bader({'charge': np.ascontiguousarray(rho)}, lat,   # e/A^3
                      np.ascontiguousarray(atoms.positions, dtype=float), info,
                      threads=4, export_mode=None, output='none')
        bader()
        nb = np.asarray(bader.atoms_charge)                     # electrons per atom
        Z = atoms.get_atomic_numbers()
        qb = Z - nb
        res['bader'] = {'per_atom': [float(q) for q in qb],
                        'molecule': float(qb[mol].sum()), 'layer': float(qb[lay].sum()),
                        'total_electrons_integrated': float(nb.sum()),
                        'total_electrons_expected': float(Z.sum()),
                        'Co': float(qb[mol][0]),
                        'by_species': {s: float(np.mean([qb[i] for i in range(len(atoms))
                                                         if atoms[i].symbol == s and i in lay]))
                                       for s in ('Zr', 'N', 'Cl')}}
        print(json.dumps({k: v for k, v in res['bader'].items() if k != 'per_atom'}, indent=1))
    update_json(res)

# --------------------------------------------------------------------- fragments
if args.stage == 'fragments':
    atoms = structure()
    lay, mol = layer_mol_indices(atoms)
    ref = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None)
    dens = {'total': ref.get_all_electron_density(gridrefinement=2)}
    efrag = {}
    for name, idx in [('slab', lay), ('molecule', mol)]:
        frag = atoms[idx]
        frag.calc = gpaw_calc(os.path.join(RUN, f'fragment_{name}.txt'), KSCF,
                              width=WIDTH, convergence=CONV, gpts=ref.wfs.gd.N_c)
        efrag[name] = frag.get_potential_energy()
        dens[name] = frag.calc.get_all_electron_density(gridrefinement=2)
    drho = dens['total'] - dens['slab'] - dens['molecule']       # e/A^3
    cell = atoms.cell.array
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    nz = drho.shape[2]
    dz = cell[2, 2] / nz
    z = np.arange(nz) * dz
    prof = drho.mean(axis=(0, 1)) * area                         # e/A
    cdc = np.cumsum(prof) * dz                                   # e
    e0 = json.load(open(os.path.join(RESULTS, 'intercalate_electronic.json')))['energy_nonspin']
    from dftd3.ase import DFTD3
    ed3 = {}
    for name, sub in [('total', atoms), ('slab', atoms[lay]), ('molecule', atoms[mol])]:
        sub = sub.copy()
        sub.calc = DFTD3(method='PBE', damping='d3bj')
        ed3[name] = sub.get_potential_energy()
    e_int_pbe = e0 - efrag['slab'] - efrag['molecule']
    e_int_d3 = ed3['total'] - ed3['slab'] - ed3['molecule']
    update_json({'E_int_frozen_frags_PBE': e_int_pbe,
                 'E_int_frozen_frags_D3_part': e_int_d3,
                 'E_int_frozen_frags_PBE_D3': e_int_pbe + e_int_d3,
                 'E_fragments': efrag, 'E_D3': ed3})
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, 'intercalate_drho.npz'), z=z, drho_z=prof,
                 cdc=cdc, zatoms=atoms.positions[:, 2],
                 symbols=np.array(atoms.get_chemical_symbols()))
