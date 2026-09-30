"""Collect the cutoff scan (runs/convergence.log) into results/convergence.json
and plot energy/force/gap convergence."""
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from style import plt, BLUE, ORANGE, AQUA, MUTED  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
rows = []
for line in open(os.path.join(ROOT, 'runs', 'convergence.log')):
    m = re.match(r'ecut (\d+) \[(\d+), (\d+), (\d+)\] E/atom (\S+) Fz \[(.*)\] gap (\S+) (\d+)s', line)
    if m:
        fz = [float(v) for v in m.group(6).split()]
        rows.append({'ecut': float(m.group(1)), 'k': [int(m.group(i)) for i in (2, 3, 4)],
                     'e_per_atom': float(m.group(5)), 'fz': fz, 'gap': float(m.group(7)),
                     'time_s': float(m.group(8))})
path = os.path.join(ROOT, 'results', 'convergence.json')
data = json.load(open(path)) if os.path.exists(path) else {}
data['ecut'] = rows
json.dump(data, open(path, 'w'), indent=2)

ec = np.array([r['ecut'] for r in rows])
e = np.array([r['e_per_atom'] for r in rows])
f = np.array([r['fz'] for r in rows])
g = np.array([r['gap'] for r in rows])
fig, ax = plt.subplots(1, 3, figsize=(9, 2.8))
ax[0].plot(ec, 1e3 * (e - e[-1]), 'o-', color=BLUE, ms=5)
ax[0].set_ylabel('E − E(800 eV) (meV/atom)')
for i, (lab, col) in enumerate([('Zr', BLUE), ('N', ORANGE), ('Cl', AQUA)]):
    ax[1].plot(ec, f[:, 2 * i], 'o-', color=col, ms=5, label=lab)
ax[1].set_ylabel('F$_z$ at unrelaxed z (eV/Å)')
ax[1].legend()
ax[2].plot(ec, g, 'o-', color=BLUE, ms=5)
ax[2].set_ylabel('PBE gap (eV)')
for a in ax:
    a.set_xlabel('PW cutoff (eV)')
    a.axvline(600, color=MUTED, lw=0.8, ls='--')
    a.grid(True, axis='y')
fig.tight_layout()
os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
fig.savefig(os.path.join(ROOT, 'figures', 'convergence_cutoff.png'), bbox_inches='tight')
print(json.dumps(rows, indent=1))
