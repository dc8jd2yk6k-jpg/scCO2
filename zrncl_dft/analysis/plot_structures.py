"""Structure panels: pristine beta-ZrNCl (hexagonal cell) and the
ZrNCl{Co(Cp)2}0.10 model (side view along the Cp-Co-Cp axis, top view).

Usage: python analysis/plot_structures.py
"""
import json
import os
import sys

import numpy as np
from ase.io import read, write

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from style import plt, INK2  # noqa: E402
from common import zrncl_hexagonal  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
FIG = os.path.join(ROOT, 'figures')
RUN = os.path.join(ROOT, 'runs', 'intercalate')
TMP = os.path.join(ROOT, 'runs', 'img')
os.makedirs(TMP, exist_ok=True)

z = json.load(open(os.path.join(ROOT, 'results', 'pristine_relax.json')))['z_pbed3']
hexcell = zrncl_hexagonal(z=z).repeat((3, 1, 1))
write(os.path.join(TMP, 'pristine.png'), hexcell, rotation='-90x', radii=0.5, scale=22)

for name in ('relaxed_pw.traj', 'relaxed_lcao.traj', 'start.traj'):
    f = os.path.join(RUN, name)
    if os.path.exists(f):
        inter = read(f)
        src = name
        break
# side view along the molecular axis (a2 direction, 120 deg from x)
big = inter.repeat((2, 1, 2))
write(os.path.join(TMP, 'inter_side.png'), big, rotation='-30z,-90x', radii=0.5, scale=22)
# top view: guest + the Cl sheet underneath it
zc = inter.cell[2, 2]
sel = [i for i, a in enumerate(inter)
       if a.symbol in ('Co', 'C', 'H') or (a.symbol == 'Cl' and a.position[2] > 0.3 * zc)]
top = inter[sel].repeat((3, 2, 1))
write(os.path.join(TMP, 'inter_top.png'), top, rotation='0x', radii=0.5, scale=22)

fig, ax = plt.subplots(1, 3, figsize=(10, 4.6),
                       gridspec_kw={'width_ratios': [1.0, 1.1, 1.3]})
for a, fn, title in [(ax[0], 'pristine.png', 'β-ZrNCl (R-3m), 3 layers'),
                     (ax[1], 'inter_side.png', 'ZrNCl{Co(Cp)₂}₀.₁₀, d = 14.7 Å'),
                     (ax[2], 'inter_top.png', 'guest layer on Cl sheet (top)')]:
    a.imshow(plt.imread(os.path.join(TMP, fn)))
    a.set_title(title, fontsize=9.5, color=INK2)
    a.axis('off')
fig.text(0.01, 0.02, 'Zr cyan · N blue · Cl green · Co pink · C grey · H white   '
         f'(model geometry: {src})', fontsize=8, color=INK2)
fig.savefig(os.path.join(FIG, 'structures.png'), bbox_inches='tight')
print('wrote figures/structures.png from', src)
