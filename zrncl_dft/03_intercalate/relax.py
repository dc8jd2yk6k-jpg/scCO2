"""Relax the ZrNCl{Co(Cp)2}0.10 model (cell fixed: a = a_exp, d = 14.7 A).

Plane waves (production cutoff) + D3(BJ), Gamma-centred 3x3x1 k mesh, which
contains the folded K/K' points (0, +-1/3) of the supercell where the ZrNCl
conduction-band minimum sits; Fermi-Dirac 0.1 eV.  Optimiser: ASE PreconLBFGS
with the Exp preconditioner (well suited to soft molecular degrees of freedom).
(An LCAO/dzp pre-relaxation was tried and abandoned: 50 s per SCF iteration in
GPAW-LCAO for this cell, slower than plane waves with OpenBLAS.)

Restarts from runs/intercalate/relax_pw.traj if present.
Run:  mpiexec -n 4 python 03_intercalate/relax.py --stage pw
"""
import argparse
import os
import sys

import numpy as np
from ase.calculators.mixing import SumCalculator
from ase.io import read, write
from ase.optimize.precon import PreconLBFGS
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import gpaw_calc, ROOT  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument('--stage', choices=['pw'], default='pw')
p.add_argument('--fmax', type=float, default=0.03)
p.add_argument('--steps', type=int, default=250)
args = p.parse_args()

RUN = os.path.join(ROOT, 'runs', 'intercalate')
KPTS = {'size': (3, 3, 1), 'gamma': True}

src = os.path.join(RUN, 'start.traj')
traj = os.path.join(RUN, 'relax_pw.traj')
if os.path.exists(traj):
    try:
        read(traj)
        src = traj
    except Exception:
        pass
atoms = read(src)
atoms.set_constraint()
if world.rank == 0:
    print('starting from', src, flush=True)

dft = gpaw_calc(os.path.join(RUN, 'relax_pw.txt'), KPTS, width=0.1,
                convergence={'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8})
from dftd3.ase import DFTD3  # noqa: E402
atoms.calc = SumCalculator([dft, DFTD3(method='PBE', damping='d3bj')])
out_traj = traj if src != traj else os.path.join(RUN, 'relax_pw_part2.traj')
opt = PreconLBFGS(atoms, precon='Exp', use_armijo=True,
                  logfile=os.path.join(RUN, 'relax_pw.log'), trajectory=out_traj)
opt.run(fmax=args.fmax, steps=args.steps)
if world.rank == 0:
    write(os.path.join(RUN, 'relaxed_pw.traj'), atoms)
    f = atoms.get_forces()
    print('final fmax %.4f  E %.5f' % (np.linalg.norm(f, axis=1).max(),
                                        atoms.get_potential_energy()))
