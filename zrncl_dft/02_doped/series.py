"""Electron doping of beta-ZrNCl by adding electrons (jellium background).

For each x (electrons per ZrNCl) the Zr2N2Cl2 cell gets GPAW charge = -2x.
Structure: relaxed pristine (frozen), plus an internal relaxation at the
target x = 0.10 to quantify the structural response.  Records E_F, the
conduction-band minimum at K, N(E_F) and the Fermi wave vector, to be compared
with the strict rigid-band shift of the pristine DOS (analysis/rigid_band.py).

Run:  mpiexec -n 4 python 02_doped/series.py
"""
import os
import sys

import numpy as np
from ase.optimize import BFGS
from ase.spacegroup.symmetrize import FixSymmetry
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import (A_EXP, C_EXP, zrncl_compact, gpaw_calc, load_json,  # noqa: E402
                    hex_special_points_cart, cart_to_scaled_k, z_params,
                    save_json, ROOT)

RUN = os.path.join(ROOT, 'runs', 'doped_series')
os.makedirs(RUN, exist_ok=True)
KPTS = {'size': (18, 18, 2), 'gamma': True}
NB = 32
XS = [0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40]
CONV = {'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8}

z0 = load_json('pristine_relax.json')['z_pbed3']
sp = hex_special_points_cart(A_EXP, C_EXP)


def analyse(atoms, x):
    """E_F, band edges, N(E_F) (tetrahedron on the SCF mesh) and k_F estimate."""
    calc = atoms.calc
    ef = calc.get_fermi_level()
    kpts = calc.get_ibz_k_points()
    eig = np.array([calc.get_eigenvalues(kpt=k) for k in range(len(kpts))])
    icb = int(np.argmin(eig[:, 24]))
    kcart = kpts[icb] @ (2 * np.pi * np.linalg.inv(atoms.cell.array).T)
    from gpaw.dos import DOSCalculator
    dc = DOSCalculator.from_calculator(calc, shift_fermi_level=False)
    nef = float(dc.raw_dos([ef], width=0.0)[0])
    area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))
    kK = np.linalg.norm(sp['K'][:2])
    return {'x': x, 'charge': -2 * x, 'energy': atoms.get_potential_energy(),
            'ef': ef, 'cbm': float(eig[icb, 24]), 'cb2_at_cbm_k': float(eig[icb, 25]),
            'cbm_k_inplane_over_K': float(np.linalg.norm(kcart[:2]) / kK),
            'vbm': float(eig[:, 23].max()),
            'ef_minus_cbm': float(ef - eig[icb, 24]),
            'gap_direct_min': float((eig[:, 24] - eig[:, 23]).min()),
            'N_EF_per_cell': nef, 'N_EF_per_ZrNCl': nef / 2,
            'kF_2D_parabolic': float(np.sqrt(np.pi * 2 * x / area)),
            'z': z_params(atoms)}


res = {'frozen': [], 'relaxed_x0.10': None}
for x in XS:
    atoms = zrncl_compact(z=z0)
    atoms.calc = gpaw_calc(os.path.join(RUN, f'x{x:.3f}.txt'), KPTS,
                           width=0.02, charge=-2 * x, nbands=NB, convergence=CONV)
    atoms.get_potential_energy()
    r = analyse(atoms, x)
    res['frozen'].append(r)
    if abs(x - 0.10) < 1e-9:
        atoms.calc.write(os.path.join(RUN, 'x0.100.gpw'))
    if world.rank == 0:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        save_json('doped_series.json', res)

# k-mesh check at x = 0.10
atoms = zrncl_compact(z=z0)
atoms.calc = gpaw_calc(os.path.join(RUN, 'x0.100_k24.txt'),
                       {'size': (24, 24, 2), 'gamma': True},
                       width=0.02, charge=-0.2, nbands=NB, convergence=CONV)
atoms.get_potential_energy()
res['x0.10_k24x24x2'] = analyse(atoms, 0.10)
if world.rank == 0:
    save_json('doped_series.json', res)

# structural response at x = 0.10 (fixed experimental lattice)
atoms = zrncl_compact(z=z0)
atoms.set_constraint(FixSymmetry(atoms, symprec=1e-3))
atoms.calc = gpaw_calc(os.path.join(RUN, 'relax_x0.100.txt'), KPTS,
                       width=0.02, charge=-0.2, nbands=NB, convergence=CONV)
BFGS(atoms, logfile=os.path.join(RUN, 'relax_x0.100.log'),
     trajectory=os.path.join(RUN, 'relax_x0.100.traj')).run(fmax=0.003, steps=40)
res['relaxed_x0.10'] = analyse(atoms, 0.10)
if world.rank == 0:
    from ase.io import write
    write(os.path.join(RUN, 'relaxed_x0.100.traj'), atoms)
    save_json('doped_series.json', res)
    print('relaxed x=0.10', res['relaxed_x0.10'])
