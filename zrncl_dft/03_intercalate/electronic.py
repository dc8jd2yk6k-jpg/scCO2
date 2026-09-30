"""Electronic structure and charge transfer of the ZrNCl{Co(Cp)2}0.10 model.

Stages (run in order; each checkpoints to runs/intercalate/):
  scf        non-spin-polarised SCF, 6x6x1 k (contains the folded K/K')
  spin       spin-polarised SCF started from 1 mu_B on Co (moment survives?)
  dos        non-SCF 9x9x1 mesh -> tetrahedron DOS / PDOS, electron count in
             the ZrNCl conduction band, N(E_F)
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
from common import gpaw_calc, RESULTS, ROOT  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument('--stage', required=True,
               choices=['scf', 'spin', 'dos', 'unfold', 'charges', 'fragments'])
p.add_argument('--structure', default=None)
args = p.parse_args()

RUN = os.path.join(ROOT, 'runs', 'intercalate')
NB = 190
NMOL = 21
KSCF = {'size': (6, 6, 1), 'gamma': True}
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
    atoms.calc = gpaw_calc(os.path.join(RUN, 'scf.txt'), KSCF, width=0.05,
                           nbands=NB, convergence=CONV)
    e = atoms.get_potential_energy()
    atoms.calc.write(os.path.join(RUN, 'scf.gpw'))
    update_json({'energy_nonspin': e, 'ef_scf': atoms.calc.get_fermi_level(),
                 'structure_used': args.structure or 'relaxed_pw/lcao'})

# --------------------------------------------------------------------- spin
if args.stage == 'spin':
    atoms = structure()
    mag = np.zeros(len(atoms))
    mag[len(atoms) - NMOL] = 1.0          # Co
    atoms.set_initial_magnetic_moments(mag)
    atoms.calc = gpaw_calc(os.path.join(RUN, 'scf_spin.txt'), KSCF, width=0.05,
                           spin=True, nbands=NB, convergence=CONV)
    e = atoms.get_potential_energy()
    mtot = atoms.calc.get_magnetic_moment()
    mloc = atoms.calc.get_magnetic_moments()
    e0 = json.load(open(os.path.join(RESULTS, 'intercalate_electronic.json')))['energy_nonspin']
    update_json({'energy_spin': e, 'E_spin_minus_nonspin': e - e0,
                 'magmom_total': mtot, 'magmom_Co': float(mloc[len(atoms) - NMOL]),
                 'magmom_abs_sum': float(np.abs(mloc).sum())})

# --------------------------------------------------------------------- dos
if args.stage == 'dos':
    calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
        kpts={'size': (9, 9, 1), 'gamma': True}, nbands=NB,
        symmetry={'point_group': False, 'time_reversal': True},
        convergence={'bands': NB - 10, 'eigenstates': 1e-5}, txt=os.path.join(RUN, 'dos.txt'))
    atoms = calc.get_atoms()
    ef = calc.get_fermi_level()
    from gpaw.dos import DOSCalculator
    dc = DOSCalculator.from_calculator(calc, shift_fermi_level=False)
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

    # state-resolved character: sum_i |<p_i|psi_nk>|^2 over all projectors
    # of the molecule's atoms vs. the layer's atoms
    nk = len(calc.get_ibz_k_points())
    eig = np.array([calc.get_eigenvalues(kpt=k) for k in range(nk)])
    wk = np.asarray(calc.get_k_point_weights())
    f = np.array([calc.get_occupation_numbers(kpt=k) for k in range(nk)])
    nel = calc.get_number_of_electrons()
    if abs(f.sum() - nel) < 1e-3:          # GPAW's f_n already carry k weights
        occ = f / wk[:, None]              # -> 0..2 per state
    else:
        occ = f
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
    mol_levels = sorted(set(np.round(eig[ismol & (eig > ef - 4) & (eig < ef + 3)], 2)))
    mol_occ = [float(np.round(v - ef, 3)) for v in mol_levels if v < ef]
    mol_emp = [float(np.round(v - ef, 3)) for v in mol_levels if v >= ef]
    nef_tot = float(dc.raw_dos([ef], width=0.0)[0])
    res = {'ef_dos': ef, 'kdos': [9, 9, 1],
           'layer_vbm': float(lay_occ_max), 'layer_cbm': cbm_layer,
           'layer_gap': float(cbm_layer - lay_occ_max),
           'ef_minus_cbm_layer': float(ef - cbm_layer),
           'n_el_above_gap_mid_layer_weighted': n_above_layer,
           'n_el_above_gap_mid_mol_weighted': n_above_mol,
           'n_el_in_layer_CB_states': n_cb_layer_states,
           'x_eff_e_per_ZrNCl': n_cb_layer_states / 10,
           'N_EF_per_cell': nef_tot, 'N_EF_per_ZrNCl': nef_tot / 10,
           'molecular_levels_occupied_rel_EF': mol_occ[-12:],
           'molecular_levels_empty_rel_EF': mol_emp[:12]}
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, 'intercalate_dos.npz'), energies=energies,
                 dos=dos, keys=list(comp), pdos=np.array(list(comp.values())),
                 ef=ef, eig=eig, wk=wk, occ=occ, frac_mol=frac_mol,
                 wmol=wmol, wlay=wlay)
        print(json.dumps(res, indent=1))
    update_json(res)

# --------------------------------------------------------------------- unfold
if args.stage == 'unfold':
    from gpaw.unfold import find_K_from_k
    atoms = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).get_atoms()
    # primitive (1x1) in-plane path G-M-K-G, reduced PC coordinates
    pts = {'G': np.array([0, 0, 0.]), 'M': np.array([0.5, 0, 0.]),
           'K': np.array([1 / 3, 1 / 3, 0.])}
    a_pc = np.linalg.solve(M.astype(float), atoms.cell.array)   # PC cell rows
    b_pc = 2 * np.pi * np.linalg.inv(a_pc).T
    path = ['G', 'M', 'K', 'G']
    kpc, x, ticks = [], [], [0.0]
    x0 = 0.0
    for s, t in zip(path[:-1], path[1:]):
        L = np.linalg.norm((pts[t] - pts[s]) @ b_pc)
        n = max(int(round(L * 45)), 4)
        for i in range(n):
            kpc.append(pts[s] + (pts[t] - pts[s]) * i / n)
            x.append(x0 + L * i / n)
        x0 += L
        ticks.append(x0)
    kpc.append(pts['G'])
    x.append(x0)
    kpc = np.array(kpc)
    KG = [find_K_from_k(k, M) for k in kpc]
    Ksc = np.array([K for K, G in KG])
    Gt = np.array([G for K, G in KG])
    calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
        kpts=Ksc, symmetry='off', nbands=NB, convergence={'bands': NB - 10, 'eigenstates': 1e-6},
        txt=os.path.join(RUN, 'unfold.txt'))
    wfs = calc.wfs
    nk = len(Ksc)
    e_kn = np.zeros((nk, NB))
    P_kn = np.zeros((nk, NB))
    Fmol_kn = np.zeros((nk, NB))
    lay, mol = layer_mol_indices(atoms)
    Minv = np.linalg.inv(M.astype(float))
    A = atoms.cell.array
    for kpt in wfs.kpt_u:
        k = kpt.k
        G_Gv = wfs.pd.get_reciprocal_vectors(q=kpt.q, add_q=False)
        G_Gc = np.rint(G_Gv @ A.T / (2 * np.pi)).astype(int)
        n_Gc = (G_Gc - Gt[k]) @ Minv.T
        mask = np.all(np.abs(n_Gc - np.rint(n_Gc)) < 1e-6, axis=1)
        C_nG = kpt.psit_nG[:]
        norm = (np.abs(C_nG) ** 2).sum(1)
        wfs.pd.gd.comm.sum(norm)
        part = (np.abs(C_nG[:, mask]) ** 2).sum(1)
        wfs.pd.gd.comm.sum(part)
        wm = np.zeros(len(part))
        wl = np.zeros(len(part))
        for a, P_ni in kpt.projections.items():
            w = (np.abs(P_ni) ** 2).sum(1)
            if a in mol:
                wm += w
            else:
                wl += w
        if wfs.pd.gd.comm.rank == 0 and wfs.bd.comm.rank == 0:
            P_kn[k, :len(part)] = part / norm
            e_kn[k] = kpt.eps_n * 27.211386245988
            Fmol_kn[k, :len(part)] = wm / np.maximum(wm + wl, 1e-12)
    wfs.kd.comm.sum(P_kn)
    wfs.kd.comm.sum(e_kn)
    wfs.kd.comm.sum(Fmol_kn)
    ef = calc.get_fermi_level()
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, 'intercalate_unfold.npz'), x=np.array(x),
                 ticks=np.array(ticks), labels=np.array(path), e_kn=e_kn, P_kn=P_kn,
                 Fmol_kn=Fmol_kn, ef=ef, kpc=kpc)

# --------------------------------------------------------------------- charges
if args.stage == 'charges':
    calc = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None)
    atoms = calc.get_atoms()
    lay, mol = layer_mol_indices(atoms)
    from gpaw.analyse.hirshfeld import HirshfeldPartitioning
    hq = np.asarray(HirshfeldPartitioning(calc).get_charges())
    rho = calc.get_all_electron_density(gridrefinement=2)       # e/A^3, incl. core
    res = {'hirshfeld': {'per_atom': [float(q) for q in hq],
                         'molecule': float(np.sum(hq[mol])),
                         'layer': float(np.sum(hq[lay]))}}
    if world.rank == 0:
        np.save(os.path.join(RUN, 'rho_ae.npy'), rho)
        from pybader.interface import Bader
        from pybader.io.cube import write as cube_write
        lat = np.ascontiguousarray(atoms.cell.array, dtype=float)
        vol = abs(np.linalg.det(lat))
        info = {'filename': 'rho_ae', 'prefix': RUN, 'file_type': 'gpaw',
                'write_function': cube_write,
                'elements': atoms.get_atomic_numbers(), 'voxel_offset': np.zeros(3)}
        bader = Bader({'charge': np.ascontiguousarray(rho * vol)}, lat,
                      np.ascontiguousarray(atoms.positions, dtype=float), info,
                      threads=4, export_mode=None)
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
                              width=0.05, convergence=CONV, gpts=ref.wfs.gd.N_c)
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
    update_json({'E_bind_frozen_frags_PBE': e0 - efrag['slab'] - efrag['molecule'],
                 'E_fragments': efrag})
    if world.rank == 0:
        np.savez(os.path.join(RESULTS, 'intercalate_drho.npz'), z=z, drho_z=prof,
                 cdc=cdc, zatoms=atoms.positions[:, 2],
                 symbols=np.array(atoms.get_chemical_symbols()))
