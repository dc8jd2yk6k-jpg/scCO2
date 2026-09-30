"""Relax pristine beta-ZrNCl with R-3m symmetry enforced (FixSymmetry).

  --mode pbed3       PBE+D3(BJ), experimental a, c fixed, z relaxed  -> 'z_pbed3'
  --mode pbe         PBE,        experimental a, c fixed, z relaxed  -> 'z_pbe'
  --mode pbed3_cell  PBE+D3(BJ), a, c and z relaxed at 800 eV (lattice check)

Results are merged into results/pristine_relax.json.
Run:  mpiexec -n 4 python 01_pristine/relax.py --mode pbed3
"""
import argparse
import json
import os
import sys

import numpy as np
from ase.calculators.mixing import SumCalculator
from ase.constraints import ExpCellFilter
from ase.io import write
from ase.optimize import BFGS
from ase.spacegroup.symmetrize import FixSymmetry
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import (A_EXP, C_EXP, ECUT, zrncl_compact, gpaw_calc, z_params,  # noqa: E402
                    RESULTS, ROOT)

p = argparse.ArgumentParser()
p.add_argument('--mode', choices=['pbe', 'pbed3', 'pbed3_cell'], required=True)
args = p.parse_args()
tag = args.mode
d3 = tag.startswith('pbed3')
cellrelax = tag.endswith('cell')

RUN = os.path.join(ROOT, 'runs', 'pristine')
os.makedirs(RUN, exist_ok=True)
KPTS = {'size': (8, 8, 2), 'gamma': True}     # insulator: converged (convergence.json)
ecut = 800.0 if cellrelax else ECUT


def geometry(atoms):
    """Zr coordination shell of the relaxed structure (A)."""
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list('ijd', atoms, 3.0)
    sym = np.array(atoms.get_chemical_symbols())
    zr = int(np.where(sym == 'Zr')[0][0])
    zrn = sorted(d[(i == zr) & (sym[j] == 'N')])
    zrcl = sorted(d[(i == zr) & (sym[j] == 'Cl')])
    ii, jj, dd = neighbor_list('ijd', atoms, 4.0)
    cl = sym == 'Cl'
    clcl = sorted(dd[cl[ii] & cl[jj]])
    return {'Zr-N': [float(x) for x in zrn], 'Zr-Cl': [float(x) for x in zrcl],
            'Cl-Cl_min': float(clcl[0]) if clcl else None}


atoms = zrncl_compact()
atoms.set_constraint(FixSymmetry(atoms, symprec=1e-3))
dft = gpaw_calc(os.path.join(RUN, f'relax_{tag}.txt'), KPTS, ecut=ecut,
                convergence={'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-9})
if d3:
    from dftd3.ase import DFTD3
    atoms.calc = SumCalculator([dft, DFTD3(method='PBE', damping='d3bj')])
else:
    atoms.calc = dft
target = ExpCellFilter(atoms) if cellrelax else atoms
opt = BFGS(target, logfile=os.path.join(RUN, f'relax_{tag}.log'),
           trajectory=os.path.join(RUN, f'relax_{tag}.traj'))
opt.run(fmax=0.003 if not cellrelax else 0.005, steps=60)

a = float(np.linalg.norm(atoms.cell[0]))
c = float(3 * atoms.cell[2, 2])
res = {'energy': atoms.get_potential_energy(), 'a': a, 'c': c, 'ecut': ecut,
       'kpts': KPTS['size'], 'z': z_params(atoms, a, c), 'geometry': geometry(atoms),
       'fmax': float(np.abs(atoms.get_forces()).max()), 'steps': opt.nsteps}
if cellrelax:
    res['stress_GPa'] = (atoms.get_stress() * 160.21766).tolist()

if world.rank == 0:
    write(os.path.join(RUN, f'relaxed_{tag}.traj'), atoms)
    path = os.path.join(RESULTS, 'pristine_relax.json')
    data = json.load(open(path)) if os.path.exists(path) else {}
    data.update({'a_exp': A_EXP, 'c_exp': C_EXP})
    data.setdefault('runs', {})[tag] = res
    if tag in ('pbe', 'pbed3'):
        data[f'z_{tag}'] = res['z']
    else:
        data['a_pbed3'], data['c_pbed3'] = a, c
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    print(tag, json.dumps(res, indent=1))
