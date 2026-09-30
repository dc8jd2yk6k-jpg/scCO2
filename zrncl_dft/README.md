# β-ZrNCl → electron-doped β-ZrNCl → ZrNCl{Co(Cp)₂}₀.₁₀: a DFT study

Companion calculations to **A. M. Fogg, V. M. Green, D. O'Hare, "Superconducting Metallocene Intercalation
Compounds of β-ZrNCl", *Chem. Mater.* 11, 216 (1999)**. That paper reports:
- Co(Cp)₂, Co(Cp′)₂ and Co(Cp*)₂ intercalate into β-ZrNCl at x ≈ 0.09–0.15 guests per ZrNCl.
- The basal spacing rises from 9.2 Å to 14.7 Å (Cp, Cp′) and 16.5 Å (Cp*), with the metallocene axis parallel to the
  layers.
- The room-temperature conductivity rises from 4×10⁻⁷ to 5.6 S cm⁻¹.
- All three compounds are type-II bulk superconductors at T_c = 14 K, the same T_c as the alkali-metal intercalates.
- The authors close by hoping that "future band structure calculations on β-ZrNCl will provide further insights."

This directory does that in three steps: (1) the pristine host, (2) the host doped with extra electrons, and (3) the full
guest–host compound.

> Method in one line: GPAW 24.1 plane-wave PAW, PBE (+D3(BJ) for geometries), 600 eV, experimental lattice.
> [METHODS.md](METHODS.md) has the details.

![structures](figures/structures.png)

## 1. Pristine β-ZrNCl

**Structure.** The internal coordinates barely move under PBE+D3 relaxation at the experimental cell.
The largest shift is 0.022 Å, on Cl.

PRISTINE_STRUCTURE_TABLE

**Electronic structure** ([figure](figures/pristine_bands_dos.png)):

![pristine](figures/pristine_bands_dos.png)

- **Band insulator.** The Zr⁴⁺ d⁰ ions give a PBE gap of **1.77 eV**, indirect from Γ (VBM) to K (CBM). The smallest
  direct gap is 2.11 eV, at K.
  - Semilocal DFT underestimates the gap. The experimental optical gap is larger, of order 3 eV.
  - The conductivity of 4×10⁻⁷ S cm⁻¹ quoted in the paper is that of an insulator.
- **Valence-band top.** Mainly N 2p (with Cl 3p); it disperses by 0.26 eV along k_z.
- **Conduction-band bottom.** Zr 4d at K: 83 % of the PAW-projected weight is Zr 4d, almost all of it in-plane
  d_xy/d_x²−y², with 16 % N 2p.
  - The next conduction band at K lies 1.10 eV higher.
  - A second valley at Γ lies 0.94 eV above the CBM.
  - Light electron doping therefore fills only the two K/K′ valleys.
- **The conduction band is two-dimensional exactly at K.**
  - E(K, k_z) varies by only 0.1 meV: in the rhombohedral (ABC) stacking the interlayer hopping cancels at K by the
    three-fold phase factors.
  - Away from K the interlayer terms grow roughly linearly. At |k−K| = 0.067 Å⁻¹, E varies by 16 meV along k_z; at
    0.134 Å⁻¹, by 31 meV.
  - The in-plane curvature at k_z = 0 gives m* = 0.85 (K→Γ) / 0.68 (K→M) mₑ. That point is the bonding bottom of the
    k_z band.
  - The k_z-averaged density of states corresponds to a **DOS mass m*_DOS = 0.58 mₑ**, i.e. N_CB = 0.27–0.29
    states eV⁻¹ ZrNCl⁻¹ just above the CBM.

## 2. Electron-doped β-ZrNCl: rigid band vs. extra electrons

Doping is modelled in two ways, x electrons per ZrNCl in both cases:

- **Strict rigid band.** The pristine bands are frozen and E_F is moved until the conduction band holds x electrons.
  E_F comes from electron counting with the tetrahedron DOS on the dense 30×30×3 mesh.
- **Extra electrons with a jellium background.** The calculation is self-consistent with −2x electrons per Zr₂N₂Cl₂.
  The bands relax around the added charge, which sits on a uniform positive background.
  - Run on a series of x values with the 18×18×2 SCF mesh.
  - At x = 0.10 it is repeated with the full band/DOS treatment (dense 30×30×3 mesh), because that is the doping
    level of ZrNCl{Co(Cp)₂}₀.₁₀.

DOPED_RESULTS

## 3. The full intercalate ZrNCl{Co(Cp)₂}₀.₁₀

INTERCALATE_RESULTS

## 4. Three levels of description side by side

COMPARISON

## 5. What this says about the paper's observations

DISCUSSION

## Caveats

- **Functional.** Semilocal PBE underestimates the host gap, and so misplaces how deep the guest levels sit relative to
  the ZrNCl bands.
  - The qualitative conclusion (complete ionisation of Co(Cp)₂) rests on an energy separation. The report quotes that
    separation below so it can be judged against a typical PBE error.
- **Ordered model.** The intercalate is an ordered, AA-stacked model at exactly x = 1/10.
  - The real compound is powder-crystalline, probably with disordered guests and its own stacking.
  - The resistivity upturn below ~70 K reported in the paper (weak localisation) is a disorder effect that this model
    cannot capture.
- **Geometry constraints.** The cell is fixed at the measured basal spacing and the experimental in-plane lattice
  constant.
  - The relaxation keeps the mirror symmetry (space group Cm) of the starting model, so the "axis parallel to the
    layers" orientation is imposed, not predicted.
- **Jellium doping.** It spreads the compensating charge uniformly, including through the van der Waals gap. It is a
  model of doping, not of any particular donor.
- **Scope.** No electron–phonon coupling was computed, so nothing here predicts T_c.

## Reproducing

All jobs were run through `runs/queue.sh` (4 MPI ranks each), in this order:

```bash
mpiexec -n 4 python 01_pristine/convergence.py            # cutoff scan -> analysis/parse_convergence.py
mpiexec -n 4 python 01_pristine/relax.py --mode pbed3     # also --mode pbe, --mode pbed3_cell
python 03_intercalate/build.py                             # ordered x = 1/10 model
mpiexec -n 4 python 01_pristine/electronic.py              # bands, DOS/PDOS, masses
mpiexec -n 4 python 03_intercalate/relax.py --stage pw     # PBE+D3 relaxation of the intercalate
mpiexec -n 4 python 02_doped/series.py                     # jellium doping series + relaxation at x=0.1
mpiexec -n 4 python 01_pristine/electronic.py --charge -0.2 --tag x0.10 --kscf 18 18 3
mpiexec -n 4 python 03_intercalate/electronic.py --stage scf   # then dos, unfold, charges, spin, fragments
mpiexec -n 4 python 03_intercalate/molecule.py             # isolated Co(Cp)2 / Co(Cp)2+
mpiexec -n 4 python 01_pristine/convergence_k.py
python 03_intercalate/geometry.py
python analysis/rigid_band.py && python analysis/plot_*.py && python analysis/summary.py
```

Software: GPAW 24.1.0, ASE 3.22.1, libxc 5.2.3, OpenBLAS 0.3.26, simple-dftd3 1.6.0, pybader, spglib 2.3.1
(Ubuntu 24.04 packages plus pip). Large restart files (`*.gpw`) and densities (`*.npy`) are not committed.
The GPAW text outputs of every run are in `runs/`.
