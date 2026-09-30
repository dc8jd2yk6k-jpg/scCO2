"""Plane-wave cutoff and k-mesh convergence for pristine beta-ZrNCl (PBE).

Run:  mpiexec -n 4 python 01_pristine/convergence.py
"""
import os
import sys
import time

import numpy as np
from ase.dft.bandgap import bandgap
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import zrncl_compact, gpaw_calc, save_json, ROOT  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'convergence')
os.makedirs(RUN, exist_ok=True)


def run(ecut, k, kz):
    atoms = zrncl_compact()
    atoms.calc = gpaw_calc(os.path.join(RUN, f'pw{int(ecut)}_k{k}x{kz}.txt'),
                           kpts={'size': (k, k, kz), 'gamma': True}, ecut=ecut)
    t0 = time.time()
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    gap, _, _ = bandgap(atoms.calc, output=None)
    return {'ecut': ecut, 'k': [k, k, kz], 'energy': e,
            'e_per_atom': e / len(atoms), 'fz': [float(x) for x in f[:, 2]],
            'fmax_z': float(np.abs(f[:, 2]).max()), 'gap': gap,
            'time_s': time.time() - t0}


def log(tag, r):
    if world.rank == 0:
        print(tag, r['ecut'], r['k'], 'E/atom %.6f' % r['e_per_atom'],
              'Fz', np.round(r['fz'], 4), 'gap %.4f' % r['gap'],
              '%.0fs' % r['time_s'], flush=True)


res = {'ecut': [], 'kpts': []}
for ecut in [400, 500, 550, 600, 700, 800]:
    r = run(ecut, 8, 2)
    res['ecut'].append(r)
    log('ecut', r)
for k, kz in [(6, 2), (8, 2), (10, 2), (12, 3), (12, 4), (16, 4)]:
    r = run(550, k, kz)
    res['kpts'].append(r)
    log('kpts', r)
if world.rank == 0:
    save_json('convergence.json', res)
