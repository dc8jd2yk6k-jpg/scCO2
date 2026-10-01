"""Progress of the symmetry-broken relaxation (03_intercalate/symmetry_check.py):
per BFGS step the energy relative to the Cm-symmetric minimum, the residual
force, and the orientation of the Cp-Co-Cp axis (tilt out of the layer plane,
in-plane rotation relative to the symmetric structure).

Usage: python analysis/symcheck_progress.py  -> results/intercalate_symcheck_progress.json
"""
import json
import os
import sys

import numpy as np
from ase.io import read

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import RESULTS, ROOT, COCP2_C  # noqa: E402

RUN = os.path.join(ROOT, 'runs', 'intercalate')
NMOL = 21


def axis(pos):
    cc = pos[-NMOL:][COCP2_C]
    a = cc[:5].mean(0) - cc[5:].mean(0)
    return a / np.linalg.norm(a)


ref = read(os.path.join(RUN, 'relaxed_pw.traj'))
e_ref = ref.get_potential_energy()
a_ref = axis(ref.positions)
a_ref_xy = a_ref[:2] / np.linalg.norm(a_ref[:2])
rows = []
for k, at in enumerate(read(os.path.join(RUN, 'symcheck.traj'), ':')):
    a = axis(at.positions)
    axy = a[:2] / np.linalg.norm(a[:2])
    rot = np.degrees(np.arctan2(a_ref_xy[0] * axy[1] - a_ref_xy[1] * axy[0], np.dot(a_ref_xy, axy)))
    rows.append({'step': k, 'dE_meV': 1e3 * (at.get_potential_energy() - e_ref),
                 'fmax': float(np.linalg.norm(at.get_forces(), axis=1).max()),
                 'tilt_deg': float(np.degrees(np.arcsin(abs(a[2])))),
                 'inplane_rotation_deg': float(rot)})
    r = rows[-1]
    print(f"{k:3d}  dE {r['dE_meV']:8.1f} meV  fmax {r['fmax']:.3f}  tilt {r['tilt_deg']:5.2f}  "
          f"rot {r['inplane_rotation_deg']:+6.2f}")
json.dump({'E_symmetric': e_ref, 'steps': rows},
          open(os.path.join(RESULTS, 'intercalate_symcheck_progress.json'), 'w'), indent=1)
