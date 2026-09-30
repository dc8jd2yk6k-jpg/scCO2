"""Build the ordered ZrNCl{Co(Cp)2}0.10 model.

One ZrNCl double layer (10 ZrNCl) + one Co(Cp)2 per gallery in the 5-cell
in-plane supercell A1 = 2a1 + a2 (|A1| = sqrt3 a), A2 = -a1 + 2a2 (sqrt7 a),
i.e. x = 1/10 exactly as found by elemental analysis (Fogg et al. 1999).
Basal spacing d = 14.7 A (XRD), layers AA-stacked, Cp-Co-Cp axis parallel to
the layers (deduced in the paper from d(CoCp2) = d(CoCp'2) < d(CoCp*2)).

The axis is placed along a2 (perpendicular to A1), so neighbouring molecules
along the short 6.24 A vector meet ring-edge to ring-edge with the pentagons
interdigitated (vertex H of one ring facing an edge of the next); the top and
bottom H of each ring point into Cl hollows of the adjacent layers.
"""
import os
import sys

import numpy as np
from ase.build import make_supercell
from ase.io import write

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import (A_EXP, D_COCP2, Z_START, SITES, zrncl_slab_cell, cocp2,  # noqa: E402
                    load_json, RESULTS, ROOT)

P5 = np.array([[2, 1, 0], [-1, 2, 0], [0, 0, 1]])


def relaxed_z():
    try:
        return load_json('pristine_relax.json')['z_pbed3']
    except (FileNotFoundError, KeyError):
        print('WARNING: relaxed z not found, using starting values')
        return dict(Z_START)


def build(z=None, a=A_EXP, d=D_COCP2, co_c=2.10, theta0=0.0):
    z = relaxed_z() if z is None else z
    slab = zrncl_slab_cell(z, a=a, d=d, z0=0.0)
    sup = make_supercell(slab, P5)
    cell = sup.cell.array
    zcl = sorted(set(np.round(sup.positions[sup.numbers == 17, 2], 4)))
    z_top, z_bot_next = zcl[-1], zcl[0] + d
    zmid = 0.5 * (z_top + z_bot_next)
    a1 = np.array([a, 0, 0])
    a2 = np.array([-a / 2, a * np.sqrt(3) / 2, 0])
    u = a2 / np.linalg.norm(a2)                  # molecular axis
    e = (2 * a1 + a2) / np.linalg.norm(2 * a1 + a2)  # in-plane ring direction
    # Cl sits on site C in both galleries walls (AA stacking), so the site-A
    # hollows along a2 are 3.60 A apart, close to the 3.42 A separation of the
    # two uppermost (lowermost) H atoms of the two rings: centre those two H
    # atoms on the midpoint between neighbouring A hollows.
    target = np.array([*(SITES['A'] @ np.array([a1[:2], a2[:2]]) + a2[:2] / 2), zmid])
    mol0 = cocp2(axis=u, ring_dir=e, co_c=co_c, theta0=theta0)
    rH = np.linalg.norm(mol0.positions[11, :] - mol0.positions[0, :]
                        - np.dot(mol0.positions[11] - mol0.positions[0], u) * u)
    center = target - rH * np.cos(np.radians(72 + theta0)) * e
    mol = cocp2(axis=u, ring_dir=e, co_c=co_c, theta0=theta0, center=center)
    atoms = sup + mol
    atoms.set_cell(cell)
    return atoms   # molecule kept whole (not wrapped) for analysis


def contacts(atoms, nmol=21, rmax=4.0):
    """Shortest guest-host and guest-guest (periodic image) contacts."""
    sym = np.array(atoms.get_chemical_symbols())
    pos = atoms.positions
    cell = atoms.cell.array
    im = len(atoms) - nmol
    host, mol = np.arange(im), np.arange(im, len(atoms))
    out = {}
    shifts = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for s in shifts:
        T = np.dot(s, cell)
        for a_ in mol:
            d_h = np.linalg.norm(pos[host] + T - pos[a_], axis=1)
            for b_, dist in zip(host, d_h):
                if dist < rmax:
                    key = f'{sym[a_]}...{sym[b_]}'
                    out[key] = min(out.get(key, rmax), dist)
            if s == (0, 0, 0):
                continue
            d_m = np.linalg.norm(pos[mol] + T - pos[a_], axis=1)
            for b_, dist in zip(mol, d_m):
                if dist < rmax:
                    key = f'{sym[a_]}...{sym[b_]} (guest-guest)'
                    out[key] = min(out.get(key, rmax), dist)
    return out


if __name__ == '__main__':
    at = build()
    print(at.get_chemical_formula(), len(at), 'atoms, cell', at.cell.cellpar().round(3))
    print('area/molecule %.2f A^2' % np.linalg.norm(np.cross(at.cell[0], at.cell[1])))
    for k, v in sorted(contacts(at).items()):
        print('  %-22s %.3f' % (k, v))
    from ase.parallel import world
    if world.rank == 0:
        out = os.path.join(ROOT, 'runs', 'intercalate')
        os.makedirs(out, exist_ok=True)
        write(os.path.join(out, 'start.traj'), at)
        write(os.path.join(RESULTS, 'intercalate_start.xyz'), at, format='extxyz')
