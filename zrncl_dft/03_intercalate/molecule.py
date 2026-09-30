"""Isolated Co(Cp)2 (neutral, S=1/2) and Co(Cp)2+ (S=0) reference calculations.

Real-space (FD) PAW, open boundary conditions, so the cation needs no
charged-cell correction.  Gives the relaxed Co-C / Co-centroid distances used
as oxidation-state fingerprints, the adiabatic/vertical ionisation energies and
the frontier orbital energies.

Run:  mpiexec -n 4 python 03_intercalate/molecule.py
"""
import os
import sys

import numpy as np
from ase.io import write
from ase.optimize import BFGS
from gpaw import GPAW, FermiDirac
from gpaw.mpi import world

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import cocp2, save_json, ROOT, XC, COCP2_C, COCP2_H  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'molecule')
os.makedirs(RUN, exist_ok=True)


def calc(tag, charge, spin):
    return GPAW(mode='fd', h=0.18, xc=XC, charge=charge, spinpol=spin,
                occupations=FermiDirac(0.01, fixmagmom=spin),
                convergence={'energy': 1e-6, 'density': 1e-5},
                txt=os.path.join(RUN, f'{tag}.txt'))


def structure_info(atoms):
    co = atoms.positions[0]
    cc = atoms.positions[COCP2_C]
    hh = atoms.positions[COCP2_H]
    dco = np.linalg.norm(cc - co, axis=1)
    cent = [cc[:5].mean(0), cc[5:].mean(0)]
    ring_cc = [np.linalg.norm(r[k] - r[(k + 1) % 5]) for r in (cc[:5], cc[5:]) for k in range(5)]
    ch = [np.linalg.norm(hh[k] - cc[k]) for k in range(10)]
    return {'Co-C_mean': float(dco.mean()), 'Co-C_min': float(dco.min()),
            'Co-C_max': float(dco.max()),
            'Co-centroid': float(np.mean([np.linalg.norm(c - co) for c in cent])),
            'C-C_mean': float(np.mean(ring_cc)), 'C-H_mean': float(np.mean(ch))}


def frontier(calc, spin):
    out = {}
    ef = calc.get_fermi_level()
    for s in range(2 if spin else 1):
        e = calc.get_eigenvalues(spin=s)
        f = calc.get_occupation_numbers(spin=s)
        occ = e[f > 0.01]
        emp = e[f < 0.99]
        out[f'spin{s}'] = {'eig': e.tolist(), 'occ': f.tolist(),
                           'homo': float(occ.max()), 'lumo': float(emp.min())}
    out['fermi'] = float(ef)
    return out


res = {}
for tag, charge, spin, mag in [('cocp2_neutral', 0, True, 1.0),
                               ('cocp2_cation', 1, False, 0.0)]:
    atoms = cocp2(co_c=2.10 if charge == 0 else 2.03)
    atoms.center(vacuum=5.5)
    if spin:
        m = np.zeros(len(atoms))
        m[0] = mag
        atoms.set_initial_magnetic_moments(m)
    atoms.calc = calc(tag, charge, spin)
    BFGS(atoms, logfile=os.path.join(RUN, f'{tag}.log'),
         trajectory=os.path.join(RUN, f'{tag}.traj')).run(fmax=0.02, steps=80)
    r = {'energy': atoms.get_potential_energy(), 'geometry': structure_info(atoms),
         'levels': frontier(atoms.calc, spin)}
    if spin:
        r['magmom'] = float(atoms.calc.get_magnetic_moment())
    res[tag] = r
    if world.rank == 0:
        write(os.path.join(RUN, f'{tag}_relaxed.traj'), atoms)
        print(tag, r['energy'], r['geometry'], flush=True)

# vertical IE: cation at the neutral geometry
from ase.io import read  # noqa: E402
neu = read(os.path.join(RUN, 'cocp2_neutral_relaxed.traj'))
neu.set_initial_magnetic_moments(None)
neu.calc = calc('cocp2_cation_at_neutral', 1, False)
e_vert = neu.get_potential_energy()
res['IE_adiabatic'] = res['cocp2_cation']['energy'] - res['cocp2_neutral']['energy']
res['IE_vertical'] = e_vert - res['cocp2_neutral']['energy']
if world.rank == 0:
    print('IE adiabatic %.3f  vertical %.3f eV' % (res['IE_adiabatic'], res['IE_vertical']))
    save_json('molecule.json', res)
