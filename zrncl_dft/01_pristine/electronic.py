"""Electronic structure of pristine beta-ZrNCl (PBE, relaxed internal z).

Also serves the electron-doped jellium runs: ``--charge Q`` adds -Q electrons
per Zr2N2Cl2 cell with a compensating uniform background (GPAW `charge=Q`).

Steps: SCF -> band structure G-M-K-G-Z -> dense-mesh DOS/PDOS (tetrahedron)
-> gap, band edges and conduction-band effective masses at K.

Run:  mpiexec -n 4 python 01_pristine/electronic.py [--charge -0.2] [--tag x0.10]
"""
import argparse
import os
import sys

import numpy as np
from ase.dft.bandgap import bandgap
from gpaw import GPAW
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import (A_EXP, C_EXP, zrncl_compact, gpaw_calc, load_json,  # noqa: E402
                    hex_special_points_cart, cart_to_scaled_k, kpath_cart,
                    save_json, ROOT)

p = argparse.ArgumentParser()
p.add_argument('--charge', type=float, default=0.0)
p.add_argument('--tag', default='pristine')
p.add_argument('--kscf', type=int, nargs=3, default=[12, 12, 3])
p.add_argument('--kdos', type=int, nargs=3, default=[30, 30, 3])
p.add_argument('--zset', default='z_pbed3')
p.add_argument('--relaxed-traj', default=None,
               help='use this structure instead of the relaxed pristine z')
args = p.parse_args()

RUN = os.path.join(ROOT, 'runs', args.tag)
os.makedirs(RUN, exist_ok=True)
NB = 40            # 24 occupied + 16 empty bands (Zr-4d manifold)

if args.relaxed_traj:
    from ase.io import read
    atoms = read(args.relaxed_traj)
    atoms.set_constraint()
else:
    z = load_json('pristine_relax.json')[args.zset]
    atoms = zrncl_compact(z=z)
cell = atoms.cell.array

# ---------------------------------------------------------------- SCF
width = 0.01 if args.charge == 0 else 0.02
atoms.calc = gpaw_calc(os.path.join(RUN, 'scf.txt'),
                       kpts={'size': tuple(args.kscf), 'gamma': True},
                       width=width, charge=args.charge, nbands=NB,
                       convergence={'energy': 1e-6, 'density': 1e-5,
                                    'eigenstates': 1e-9})
etot = atoms.get_potential_energy()
ef_scf = atoms.calc.get_fermi_level()
atoms.calc.write(os.path.join(RUN, 'scf.gpw'))

# ---------------------------------------------------------------- bands
sp = hex_special_points_cart(A_EXP, C_EXP)
sp['Z'] = np.array([0, 0, 3 * np.pi / C_EXP])     # (0,0,3/2) hex = RHL Z point
labels = ['G', 'M', 'K', 'G', 'Z']
kc, x, ticks = kpath_cart([sp[s] for s in labels], npts_per_invA=70)
bs = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
    kpts=cart_to_scaled_k(kc, cell), symmetry='off', nbands=NB,
    convergence={'bands': NB - 4, 'eigenstates': 1e-6}, txt=os.path.join(RUN, 'bands.txt'))
eps_band = np.array([bs.get_eigenvalues(kpt=k) for k in range(len(kc))])

# fine lines around K for effective masses: K->G and K->M, |dk| <= 0.08 1/A
dirs = {'K-G': (sp['G'] - sp['K']) / np.linalg.norm(sp['G'] - sp['K']),
        'K-M': (sp['M'] - sp['K']) / np.linalg.norm(sp['M'] - sp['K'])}
dks = np.linspace(0, 0.08, 9)
kk = np.array([sp['K'] + d * v for v in dirs.values() for d in dks])
# kz dispersion of the band edges at K and G
kz_pts = np.array([sp['K'], sp['K'] + sp['Z'], sp['G'], sp['Z']])
kk = np.vstack([kk, kz_pts])
fm = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
    kpts=cart_to_scaled_k(kk, cell), symmetry='off', nbands=NB,
    convergence={'bands': NB - 4, 'eigenstates': 1e-7}, txt=os.path.join(RUN, 'kfine.txt'))
eps_fine = np.array([fm.get_eigenvalues(kpt=k) for k in range(len(kk))])

# ---------------------------------------------------------------- dense DOS
dos = GPAW(os.path.join(RUN, 'scf.gpw'), txt=None).fixed_density(
    kpts={'size': tuple(args.kdos), 'gamma': True}, nbands=NB,
    convergence={'bands': NB - 4, 'eigenstates': 1e-5}, txt=os.path.join(RUN, 'dos.txt'))
ef_dos = dos.get_fermi_level()
dos.write(os.path.join(RUN, 'dos.gpw'))

from gpaw.dos import DOSCalculator  # noqa: E402
dc = DOSCalculator.from_calculator(dos, shift_fermi_level=False)
energies = np.linspace(ef_dos - 12, ef_dos + 6, 3601)
d_tot = dc.raw_dos(energies, width=0.0)
sym = atoms.get_chemical_symbols()
pdos = {}
for s, l in [('Zr', 2), ('Zr', 1), ('Zr', 0), ('N', 1), ('N', 0), ('Cl', 1), ('Cl', 0)]:
    pdos[f'{s}-{"spd"[l]}'] = sum(dc.raw_pdos(energies, a=i, l=l, width=0.0)
                                  for i, t in enumerate(sym) if t == s)
# Zr-d orbital resolved (real spherical harmonics m=0..4 in GPAW's order:
# xy, yz, 3z2-r2, zx, x2-y2)
for m, name in enumerate(['dxy', 'dyz', 'dz2', 'dzx', 'dx2-y2']):
    pdos[f'Zr-{name}'] = sum(dc.raw_pdos(energies, a=i, l=2, m=m, width=0.0)
                             for i, t in enumerate(sym) if t == 'Zr')

# ---------------------------------------------------------------- analysis
nel = dos.get_number_of_electrons()
nval = 24   # occupied bands of neutral Zr2N2Cl2 (48 electrons)
eig_dense = np.array([dos.get_eigenvalues(kpt=k) for k in range(len(dos.get_ibz_k_points()))])
vbm, cbm = eig_dense[:, nval - 1].max(), eig_dense[:, nval].min()
res = {'tag': args.tag, 'charge_per_cell': args.charge,
       'x_e_per_ZrNCl': -args.charge / 2, 'energy': etot,
       'ef_scf': ef_scf, 'ef_dos': ef_dos, 'n_electrons': nel,
       'vbm_dense': float(vbm), 'cbm_dense': float(cbm), 'gap_dense': float(cbm - vbm),
       'kscf': args.kscf, 'kdos': args.kdos}
if args.charge == 0:
    g, p1, p2 = bandgap(dos, output=None)
    res['bandgap_ase'] = g
    res['vbm_k'] = dos.get_ibz_k_points()[p1[1]].tolist()
    res['cbm_k'] = dos.get_ibz_k_points()[p2[1]].tolist()
    gd, _, _ = bandgap(dos, direct=True, output=None)
    res['direct_gap_min'] = gd
# band edges on the path
ib = eps_band[:, nval]
res['cb_path_min'] = float(ib.min())
res['cb_path_min_x'] = float(x[ib.argmin()])
res['vb_path_max'] = float(eps_band[:, nval - 1].max())
res['vb_path_max_x'] = float(x[eps_band[:, nval - 1].argmax()])
# effective masses (parabola through |dk| <= 0.05 1/A)
hb2m = 3.80998  # hbar^2/2m_e in eV A^2
masses = {}
for j, name in enumerate(dirs):
    e = eps_fine[j * len(dks):(j + 1) * len(dks), nval]
    sel = dks <= 0.05 + 1e-9
    c2 = np.polyfit(dks[sel] ** 2, e[sel], 1)[0]
    masses[name] = hb2m / c2
res['m_eff_CB_at_K'] = masses
res['m_eff_CB_at_K_avg'] = float(np.mean(list(masses.values())))
ek = eps_fine[-4:]
res['kz_dispersion'] = {'CB_K(kz=0)': float(ek[0, nval]), 'CB_K(kz=Z)': float(ek[1, nval]),
                        'VB_G(kz=0)': float(ek[2, nval - 1]), 'VB_Z': float(ek[3, nval - 1])}

if world.rank == 0:
    np.savez(os.path.join(ROOT, 'results', f'{args.tag}_electronic.npz'),
             x=x, ticks=ticks, labels=labels, eps_band=eps_band, energies=energies,
             dos=d_tot, pdos_keys=list(pdos), pdos=np.array(list(pdos.values())),
             ef=ef_dos, ef_scf=ef_scf, eig_dense=eig_dense,
             ibzk=dos.get_ibz_k_points(), ibzw=dos.get_k_point_weights(),
             eps_fine=eps_fine, dks=dks, cell=cell)
    save_json(f'{args.tag}_electronic.json', res)
    print(res)
