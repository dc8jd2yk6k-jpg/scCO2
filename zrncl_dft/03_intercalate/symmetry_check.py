"""Is the Cm-symmetric relaxed intercalate a true minimum?

The BFGS relaxation keeps the mirror plane (Cm) of the starting model, which
pins the Cp-Co-Cp axis parallel to the layers.  Here the guest is rotated
rigidly, by 8 deg out of the layer plane and by 6 deg about the layer normal,
which breaks the mirror, and the structure is relaxed again (PBE+D3, same
settings).  A first attempt that also gave every guest atom a random 0.03 A
kick started at 5.6 eV/A (stretched C-H bonds) and was stopped: the rigid
rotation alone tests the orientation at a fraction of the cost.
Progress: analysis/symcheck_progress.py.  Result -> results/intercalate_symcheck.json

Run:  mpiexec -n 4 python 03_intercalate/symmetry_check.py
"""
import json
import os
import sys

import numpy as np
from ase.calculators.mixing import SumCalculator
from ase.io import read, write
from ase.optimize import BFGS
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))
from common import gpaw_calc, RESULTS, ROOT  # noqa: E402
from geometry import cocp2_info  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'intercalate')
NMOL = 21
KPTS = {'size': (3, 3, 1), 'gamma': True}


def rotate(pos, center, axis, deg):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    t = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
    return (pos - center) @ R.T + center


ref = read(os.path.join(RUN, 'relaxed_pw.traj'))
e_ref = ref.get_potential_energy()
atoms = ref.copy()
mol = np.arange(len(atoms) - NMOL, len(atoms))
pos = atoms.positions.copy()
co = pos[mol[0]]
ax = pos[mol[1:6]].mean(0) - pos[mol[11:16]].mean(0)      # ring-1 minus ring-2 centroid
ax /= np.linalg.norm(ax)
tilt_axis = np.cross(ax, [0, 0, 1])                         # in-plane, perpendicular to axis
pos[mol] = rotate(pos[mol], co, tilt_axis, 8.0)
pos[mol] = rotate(pos[mol], co, [0, 0, 1], 6.0)
atoms.positions = pos

dft = gpaw_calc(os.path.join(RUN, 'symcheck.txt'), KPTS, width=0.1,
                convergence={'energy': 1e-6, 'density': 1e-4, 'eigenstates': 4e-8})
from dftd3.ase import DFTD3  # noqa: E402
atoms.calc = SumCalculator([dft, DFTD3(method='PBE', damping='d3bj')])
opt = BFGS(atoms, maxstep=0.1, logfile=os.path.join(RUN, 'symcheck.log'),
           trajectory=os.path.join(RUN, 'symcheck.traj'))
opt.run(fmax=0.05, steps=40)      # capped: ~4 min per step; the orientation is what matters
e = atoms.get_potential_energy()
res = {'E_symmetric': e_ref, 'E_broken_start_relaxed': e, 'dE_meV': 1e3 * (e - e_ref),
       'steps': opt.nsteps, 'fmax': float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
       'guest_symmetric': cocp2_info(ref.positions[-NMOL:]),
       'guest_broken_relaxed': cocp2_info(atoms.positions[-NMOL:])}
write(os.path.join(RUN, 'symcheck_relaxed.traj'), atoms)     # all ranks (ASE handles I/O)
if world.rank == 0:
    json.dump(res, open(os.path.join(RESULTS, 'intercalate_symcheck.json'), 'w'), indent=2)
    print(json.dumps(res, indent=1))
