"""Shared settings and structure builders for the beta-ZrNCl DFT workflow.

Code: GPAW 24.1 (PAW, plane waves), ASE 3.22, PBE; PBE+D3(BJ) for relaxations
with a molecular guest.  All energies in eV, lengths in Angstrom.

Structure of beta-ZrNCl: SmSI type, R-3m (No. 166), all atoms on 6c (0,0,z).
Lattice constants are the experimental room-temperature values used throughout
the beta-ZrNCl literature (a = 3.6046 A, c = 27.672 A, hexagonal setting); the
internal z parameters are only starting guesses and are relaxed with DFT.
"""
import json
import os

import numpy as np
from ase import Atoms
from ase.spacegroup import crystal

ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, 'results')
os.makedirs(RESULTS, exist_ok=True)

# --- experimental cell of pristine beta-ZrNCl (hexagonal setting of R-3m) ---
A_EXP = 3.6046
C_EXP = 27.672
Z_START = {'Zr': 0.1196, 'N': 0.1981, 'Cl': 0.3888}   # starting guesses

# basal (interlayer) spacing of ZrNCl{Co(Cp)2}0.10, Fogg et al. Chem. Mater. 1999
D_COCP2 = 14.7

# --- numerical settings (see results/convergence.json) ---
ECUT = 600.0            # plane-wave cutoff, eV (see results/convergence.json)
XC = 'PBE'


def hex_cell(a=A_EXP, c=C_EXP):
    """Hexagonal lattice vectors a1=(a,0,0), a2=(-a/2, a*sqrt3/2, 0), a3=(0,0,c)."""
    return np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])


def zrncl_primitive(a=A_EXP, c=C_EXP, z=None):
    """Rhombohedral primitive cell of beta-ZrNCl: Zr2N2Cl2 (one ZrN double layer)."""
    z = dict(Z_START if z is None else z)
    at = crystal(['Zr', 'N', 'Cl'],
                 basis=[(0, 0, z['Zr']), (0, 0, z['N']), (0, 0, z['Cl'])],
                 spacegroup=166, cellpar=[a, a, c, 90, 90, 120],
                 primitive_cell=True)
    return Atoms(at.get_chemical_symbols(), positions=at.positions,
                 cell=at.cell, pbc=True)


def compact_cell(a=A_EXP, c=C_EXP):
    """Primitive cell of the R-3m lattice spanned by the in-plane hexagonal
    vectors a1, a2 and the centring vector a3 = (2/3, 1/3, 1/3)_hex.

    Same lattice as the symmetric rhombohedral cell (alpha ~ 22 deg), but the
    plane-wave FFT grid is ~7x smaller because a1, a2 are short.  Used for all
    production calculations on the pristine/doped host (without point-group
    symmetry, see gpaw_calc(); GPAW would otherwise symmetrise the grid back
    to N1 = N2 = N3).
    """
    return np.array([[a, 0, 0],
                     [-a / 2, a * np.sqrt(3) / 2, 0],
                     [a / 2, a * np.sqrt(3) / 6, c / 3]])


def zrncl_compact(a=A_EXP, c=C_EXP, z=None):
    """beta-ZrNCl (Zr2N2Cl2) in the compact primitive cell (see compact_cell)."""
    rh = zrncl_primitive(a, c, z)
    at = Atoms(rh.get_chemical_symbols(), positions=rh.positions,
               cell=compact_cell(a, c), pbc=True)
    at.wrap()
    return at


def zrncl_hexagonal(a=A_EXP, c=C_EXP, z=None):
    """Conventional hexagonal cell (18 atoms, 3 double layers)."""
    z = dict(Z_START if z is None else z)
    at = crystal(['Zr', 'N', 'Cl'],
                 basis=[(0, 0, z['Zr']), (0, 0, z['N']), (0, 0, z['Cl'])],
                 spacegroup=166, cellpar=[a, a, c, 90, 90, 120])
    return Atoms(at.get_chemical_symbols(), positions=at.positions,
                 cell=at.cell, pbc=True)


def z_params(atoms, a=A_EXP, c=C_EXP):
    """Recover the 6c z parameters (hexagonal setting) of Zr, N, Cl.

    Each atom is mapped back onto the (0,0,z) line with the R-centring
    translations (2/3,1/3,1/3); of the +-z pair the value in (0, 1/2) is kept.
    Only meaningful for a structure that still has R-3m symmetry.
    """
    s_hex = np.dot(atoms.positions, np.linalg.inv(hex_cell(a, c)))
    out = {}
    for s in ('Zr', 'N', 'Cl'):
        vals = []
        for sh, sym in zip(s_hex, atoms.get_chemical_symbols()):
            if sym != s:
                continue
            for m in range(3):
                d = sh[:2] - m * np.array([2 / 3, 1 / 3])
                d -= np.round(d)
                if np.allclose(d, 0, atol=1e-3):
                    vals.append((sh[2] - m / 3) % 1.0)
        vals = [v if v < 0.5 else 1 - v for v in vals]
        out[s] = float(np.mean(vals))
    return out


# --- ZrNCl double layer and the Co(Cp)2 guest ---------------------------------
SITES = {'A': np.array([0.0, 0.0]), 'B': np.array([2 / 3, 1 / 3]),
         'C': np.array([1 / 3, 2 / 3])}


def slab_layers(z, c=C_EXP):
    """Atomic planes of one Cl-Zr-N-N-Zr-Cl double layer of the R-3m host:
    list of (symbol, in-plane site, height in A), bottom Cl at the lowest z.
    Derived from the 6c positions (0,0,+-z) and the R-centring (2/3,1/3,1/3).
    """
    planes = [('Cl', 'C', (z['Cl'] - 1 / 3) * c),
              ('Zr', 'A', z['Zr'] * c),
              ('N', 'B', (1 / 3 - z['N']) * c),
              ('N', 'A', z['N'] * c),
              ('Zr', 'B', (1 / 3 - z['Zr']) * c),
              ('Cl', 'C', (2 / 3 - z['Cl']) * c)]
    h0 = planes[0][2]
    return [(s, site, h - h0) for s, site, h in planes]


def zrncl_slab_cell(z, a=A_EXP, c_host=C_EXP, d=D_COCP2, z0=None):
    """Single ZrNCl double layer in a hexagonal cell of height d (AA-stacked
    layers separated by galleries).  Bottom Cl plane at z0 (default: centred)."""
    planes = slab_layers(z, c_host)
    thick = planes[-1][2]
    if z0 is None:
        z0 = 0.5 * (d - thick)
    cell = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, d]])
    syms, pos = [], []
    for s, site, h in planes:
        xy = SITES[site] @ cell[:2, :2]
        syms.append(s)
        pos.append([xy[0], xy[1], z0 + h])
    return Atoms(syms, positions=pos, cell=cell, pbc=True)


# atom order produced by cocp2(): Co, ring-1 C x5, ring-1 H x5, ring-2 C x5, ring-2 H x5
COCP2_C = [1, 2, 3, 4, 5, 11, 12, 13, 14, 15]
COCP2_H = [6, 7, 8, 9, 10, 16, 17, 18, 19, 20]


def cocp2(axis=(1, 0, 0), ring_dir=(0, 1, 0), co_c=2.10, cc=1.43, ch=1.09,
          theta0=0.0, staggered=False, center=(0, 0, 0)):
    """Co(Cp)2 with its Cp-Co-Cp axis along `axis`.

    Ring atoms sit at centroid + r(cos t * ring_dir + sin t * n) with
    n = axis x ring_dir; t = theta0 + 72k (deg).  Eclipsed (D5h) by default.
    """
    u = np.asarray(axis, float)
    u /= np.linalg.norm(u)
    e = np.asarray(ring_dir, float)
    e -= u * np.dot(e, u)
    e /= np.linalg.norm(e)
    n = np.cross(u, e)
    r = cc / (2 * np.sin(np.pi / 5))
    h = np.sqrt(co_c ** 2 - r ** 2)
    syms, pos = ['Co'], [np.zeros(3)]
    for side in (+1, -1):
        off = 36.0 if (staggered and side < 0) else 0.0
        for k in range(5):
            t = np.radians(theta0 + off + 72 * k)
            d = np.cos(t) * e + np.sin(t) * n
            syms.append('C')
            pos.append(side * h * u + r * d)
        for k in range(5):
            t = np.radians(theta0 + off + 72 * k)
            d = np.cos(t) * e + np.sin(t) * n
            syms.append('H')
            pos.append(side * h * u + (r + ch) * d)
    return Atoms(syms, positions=np.array(pos) + np.asarray(center))


def gpaw_calc(txt, kpts, ecut=ECUT, width=0.01, charge=0.0, spin=False,
              symmetry='time-reversal-only', **kw):
    """GPAW plane-wave calculator with the project-wide settings.

    symmetry='time-reversal-only' keeps k -> -k but drops point-group
    operations, which lets GPAW use the small FFT grid of compact_cell().
    """
    from gpaw import GPAW, PW, FermiDirac
    if symmetry == 'time-reversal-only':
        symmetry = {'point_group': False, 'time_reversal': True}
    params = dict(mode=PW(ecut), xc=XC, kpts=kpts,
                  occupations=FermiDirac(width), charge=charge,
                  spinpol=spin, symmetry=symmetry,
                  convergence={'energy': 1e-6, 'density': 1e-6,
                               'eigenstates': 1e-9},
                  txt=txt)
    params.update(kw)
    return GPAW(**params)


# --- band unfolding -----------------------------------------------------------
def unfold_path(M, a_pc, labels=('G', 'M', 'K', 'G'), dens=45):
    """Primitive-cell k path (reduced PC coords), x axis, ticks, and the
    corresponding supercell k points K (reduced SC coords) with the integer
    shifts G_t such that M k = K + G_t (gpaw.unfold.find_K_from_k)."""
    from gpaw.unfold import find_K_from_k
    pts = {'G': np.array([0, 0, 0.]), 'M': np.array([0.5, 0, 0.]),
           'K': np.array([1 / 3, 1 / 3, 0.])}
    b_pc = 2 * np.pi * np.linalg.inv(a_pc).T
    kpc, x, ticks, x0 = [], [], [0.0], 0.0
    for s, t in zip(labels[:-1], labels[1:]):
        L = np.linalg.norm((pts[t] - pts[s]) @ b_pc)
        n = max(int(round(L * dens)), 4)
        for i in range(n):
            kpc.append(pts[s] + (pts[t] - pts[s]) * i / n)
            x.append(x0 + L * i / n)
        x0 += L
        ticks.append(x0)
    kpc.append(pts[labels[-1]])
    x.append(x0)
    kpc = np.array(kpc)
    KG = [find_K_from_k(k, M) for k in kpc]
    return kpc, np.array(x), ticks, np.array([K for K, G in KG]), \
        np.array([G for K, G in KG])


def unfold_weights(calc, M, Gt, mol_atoms=()):
    """Spectral weights P_n(k) = sum_{g in PC lattice} |C_n(k - K + g)|^2 of the
    supercell pseudo wave functions (PW coefficients), for a calculator whose
    k-point list is the unfolding list K (one per PC k, symmetry off).

    Also returns eigenvalues (eV) and the fraction of PAW-projector weight on
    `mol_atoms`.  Collective over the calculator's communicators."""
    wfs = calc.wfs
    nk, nb = len(Gt), wfs.bd.nbands
    e_kn, P_kn, F_kn = np.zeros((nk, nb)), np.zeros((nk, nb)), np.zeros((nk, nb))
    Minv = np.linalg.inv(np.asarray(M, float))
    A = calc.wfs.gd.cell_cv          # bohr, the units of GPAW's G vectors
    mol_atoms = set(mol_atoms)
    for kpt in wfs.kpt_u:
        k = kpt.k
        G_Gv = wfs.pd.get_reciprocal_vectors(q=kpt.q, add_q=False)
        G_Gc = np.rint(G_Gv @ A.T / (2 * np.pi)).astype(int)
        n_Gc = (G_Gc - Gt[k]) @ Minv.T
        mask = np.all(np.abs(n_Gc - np.rint(n_Gc)) < 1e-6, axis=1)
        C_nG = kpt.psit_nG[:]
        norm = (np.abs(C_nG) ** 2).sum(1)
        part = (np.abs(C_nG[:, mask]) ** 2).sum(1)
        wfs.pd.gd.comm.sum(norm)
        wfs.pd.gd.comm.sum(part)
        wm, wl = np.zeros(len(part)), np.zeros(len(part))
        for a, P_ni in kpt.projections.items():
            w = (np.abs(P_ni) ** 2).sum(1)
            if a in mol_atoms:
                wm += w
            else:
                wl += w
        if wfs.pd.gd.comm.rank == 0 and wfs.bd.comm.rank == 0:
            P_kn[k, :len(part)] = part / norm
            e_kn[k] = kpt.eps_n * 27.211386245988
            F_kn[k, :len(part)] = wm / np.maximum(wm + wl, 1e-12)
    for arr in (P_kn, e_kn, F_kn):
        wfs.kd.comm.sum(arr)
    return e_kn, P_kn, F_kn


# --- k-space helpers --------------------------------------------------------
def hex_special_points_cart(a=A_EXP, c=C_EXP):
    """Cartesian special points of the 2D hexagonal BZ (kz=0) plus A, in 1/A."""
    b = 2 * np.pi * np.linalg.inv(hex_cell(a, c)).T
    return {'G': np.zeros(3),
            'M': 0.5 * b[0],
            'K': (b[0] + b[1]) / 3.0,
            'A': 0.5 * b[2]}


def cart_to_scaled_k(kcart, cell):
    """Convert Cartesian k (1/A, including 2pi) to fractional coordinates of the
    reciprocal basis of `cell` (GPAW's convention for explicit k-point lists)."""
    return np.dot(kcart, np.asarray(cell).T) / (2 * np.pi)


def kpath_cart(points, npts_per_invA=60):
    """Piecewise-linear path through Cartesian points; returns (kcart, x, xticks)."""
    ks, xs, ticks = [], [], [0.0]
    x0 = 0.0
    for p, q in zip(points[:-1], points[1:]):
        L = np.linalg.norm(q - p)
        n = max(int(L * npts_per_invA), 2)
        for i in range(n):
            t = i / n
            ks.append(p + t * (q - p))
            xs.append(x0 + t * L)
        x0 += L
        ticks.append(x0)
    ks.append(points[-1])
    xs.append(x0)
    return np.array(ks), np.array(xs), ticks


def save_json(name, data):
    path = os.path.join(RESULTS, name)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=float)
    return path


def load_json(name):
    with open(os.path.join(RESULTS, name)) as f:
        return json.load(f)
