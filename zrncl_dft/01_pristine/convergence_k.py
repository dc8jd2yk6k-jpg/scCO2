"""k-mesh convergence of pristine beta-ZrNCl at the production cutoff.

Complements the cutoff scan of convergence.py (whose log is parsed into
results/convergence.json); run after the cutoff scan.
Run:  mpiexec -n 4 python 01_pristine/convergence_k.py
"""
import json
import os
import sys
import time

import numpy as np
from ase.dft.bandgap import bandgap
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import ECUT, zrncl_compact, gpaw_calc, RESULTS, ROOT  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'convergence')
os.makedirs(RUN, exist_ok=True)
out = []
for k, kz in [(6, 2), (8, 2), (12, 3)]:
    atoms = zrncl_compact()
    atoms.calc = gpaw_calc(os.path.join(RUN, f'pw{int(ECUT)}_k{k}x{kz}.txt'),
                           kpts={'size': (k, k, kz), 'gamma': True},
                           convergence={'energy': 1e-6, 'density': 1e-5,
                                        'eigenstates': 1e-9})
    t0 = time.time()
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    gap, _, _ = bandgap(atoms.calc, output=None)
    r = {'ecut': ECUT, 'k': [k, k, kz], 'e_per_atom': e / len(atoms),
         'fz': [float(v) for v in f[:, 2]], 'gap': gap, 'time_s': time.time() - t0}
    out.append(r)
    if world.rank == 0:
        print('kpts', r, flush=True)
if world.rank == 0:
    path = os.path.join(RESULTS, 'convergence.json')
    data = json.load(open(path)) if os.path.exists(path) else {}
    data['kpts'] = out
    json.dump(data, open(path, 'w'), indent=2)
