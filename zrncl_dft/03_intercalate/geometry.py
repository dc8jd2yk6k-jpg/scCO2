"""Geometry analysis of the relaxed ZrNCl{Co(Cp)2}0.10 model -> results/intercalate_geometry.json

Guest: Co-C, Co-centroid, C-C, ring tilt, axis orientation relative to the
layers; guest-host contacts; host: Zr-N / Zr-Cl bond statistics and slab
thickness compared with the pristine double layer.

Usage: python 03_intercalate/geometry.py [structure.traj]
"""
import json
import os
import sys

import numpy as np
from ase.io import read

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common import RESULTS, ROOT, COCP2_C  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from build import contacts  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'intercalate')
NMOL = 21


def cocp2_info(pos):
    co, cc = pos[0], pos[COCP2_C]
    d = np.linalg.norm(cc - co, axis=1)
    cent = [cc[:5].mean(0), cc[5:].mean(0)]
    axis = cent[0] - cent[1]
    axis /= np.linalg.norm(axis)
    # ring normals
    normals = []
    for ring in (cc[:5], cc[5:]):
        u, s, vt = np.linalg.svd(ring - ring.mean(0))
        normals.append(vt[2])
    tilt = [np.degrees(np.arccos(min(1.0, abs(np.dot(n, axis))))) for n in normals]
    ring_cc = [np.linalg.norm(r[k] - r[(k + 1) % 5]) for r in (cc[:5], cc[5:]) for k in range(5)]
    return {'Co-C_mean': float(d.mean()), 'Co-C_min': float(d.min()), 'Co-C_max': float(d.max()),
            'Co-centroid': float(np.mean([np.linalg.norm(c - co) for c in cent])),
            'C-C_mean': float(np.mean(ring_cc)),
            'axis_angle_to_layer_plane_deg': float(np.degrees(np.arcsin(abs(axis[2])))),
            'ring_normal_vs_axis_deg': [float(t) for t in tilt]}


def host_info(at):
    from ase.neighborlist import neighbor_list
    host = at[:len(at) - NMOL]
    i, j, d = neighbor_list('ijd', host, 3.0)
    sym = np.array(host.get_chemical_symbols())
    zrn = d[(sym[i] == 'Zr') & (sym[j] == 'N')]
    zrcl = d[(sym[i] == 'Zr') & (sym[j] == 'Cl')]
    zcl = host.positions[sym == 'Cl', 2]
    zzr = host.positions[sym == 'Zr', 2]
    return {'Zr-N_mean': float(zrn.mean()), 'Zr-N_min': float(zrn.min()),
            'Zr-N_max': float(zrn.max()), 'Zr-Cl_mean': float(zrcl.mean()),
            'Zr-Cl_min': float(zrcl.min()), 'Zr-Cl_max': float(zrcl.max()),
            'Cl-Cl_slab_thickness': float(np.ptp(zcl)),
            'Zr-Zr_interplane': float(np.ptp(zzr))}


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RUN, 'relaxed_pw.traj')
    at = read(src)
    mol_pos = at.positions[-NMOL:]
    zcl = at.positions[:len(at) - NMOL][np.array(at.get_chemical_symbols()[:len(at) - NMOL]) == 'Cl', 2]
    zco = mol_pos[0, 2]
    res = {'structure': os.path.relpath(src, ROOT), 'CoCp2': cocp2_info(mol_pos),
           'host': host_info(at),
           'contacts': {k: float(v) for k, v in contacts(at).items()},
           'Co_height_above_top_Cl_plane': float(zco - zcl.max()) if zco > zcl.max()
           else float(zco + at.cell[2, 2] - zcl.max())}
    start = read(os.path.join(RUN, 'start.traj'))
    res['start'] = {'CoCp2': cocp2_info(start.positions[-NMOL:]), 'host': host_info(start)}
    json.dump(res, open(os.path.join(RESULTS, 'intercalate_geometry.json'), 'w'), indent=2)
    print(json.dumps(res, indent=1))
