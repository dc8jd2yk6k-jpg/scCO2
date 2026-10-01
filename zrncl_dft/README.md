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

![doping](figures/doping_rigid_vs_scf.png)

| x (e/ZrNCl) | E_F−E_CBM rigid (meV) | E_F−E_CBM SCF (meV) | N(E_F) rigid | N(E_F) SCF | k_F 2D (1/Å) |
|---|---|---|---|---|---|
| 0.025 | 100 | 112 | 0.283 | 0.239 | 0.118 |
| 0.050 | 192 | 209 | 0.282 | 0.278 | 0.167 |
| 0.100 | 371 | 390 | 0.289 | 0.276 | 0.236 |
| 0.150 | 542 | 567 | 0.302 | 0.302 | 0.289 |
| 0.200 | 702 | 731 | 0.324 | 0.330 | 0.334 |
| 0.300 | 954 | 935 | 1.127 | 1.251 | 0.409 |
| 0.400 | 1032 | 996 | 1.520 | 1.568 | 0.473 |

N(E_F) in states eV⁻¹ ZrNCl⁻¹ (both spins); self-consistent values from tetrahedron electron counting on the 18×18×2 SCF mesh (at x = 0.10 a 24×24×2 SCF gives 382 meV / 0.281, the dense-mesh result below 378 meV / 0.284).

**Findings.**

1. **The rigid-band picture holds well for the filling.**
   - For x ≤ 0.2 the self-consistent E_F − E_CBM exceeds the rigid-band value by only 4–8 % (at x = 0.10: 390 meV vs
     371 meV).
   - N(E_F) is the same within the mesh accuracy.
   - At x = 0.10 the carriers fill two small, nearly circular pockets at K and K′ with k_F ≈ 0.24 Å⁻¹. Each pocket
     covers 5 % of the Brillouin zone.
   - N(E_F) ≈ 0.28 states eV⁻¹ ZrNCl⁻¹ and is flat in x up to x ≈ 0.25, the two-dimensional constant-DOS signature.
   - Near x ≈ 0.28 E_F reaches the Γ valley (0.9 eV above the K minimum) and N(E_F) jumps three- to fourfold.
2. **Where the rigid band fails.**
   - The self-consistent charge separates the host bands: the K-conduction-band bottom rises relative to the Γ-valence
     top by about +65 meV per 0.1 e⁻, from 1.77 eV at x = 0 to 2.02 eV at x = 0.4.
   - The Γ valley comes down relative to K by up to 50 meV.
   - Both effects are electrostatic: the added electrons sit in the Zr planes and the compensating charge is spread
     through the van der Waals gap. In the real intercalate that charge sits on the guests in the gallery.
3. **The lattice responds.** Relaxing the internal coordinates at x = 0.10 changes the bonds as follows:
   - Zr–Cl lengthens by 0.020 Å and Zr–N shortens by 0.006 Å.
   - The energy drops by 18 meV per cell.
   - The Γ valley moves up to 1.06 eV above the K minimum, while E_F − E_CBM barely changes (398 meV).

**x = 0.10 on the dense mesh.** The SCF used 18×18×3 k-points; the bands and DOS were evaluated on 30×30×3 with the
tetrahedron method.

| | strict rigid band | extra e⁻ + jellium (frozen lattice) |
|---|---|---|
| E_F − E_CBM(K) | 371 meV | 378 meV |
| N(E_F) (states eV⁻¹ ZrNCl⁻¹, both spins) | 0.289 | 0.284 |
| gap VBM(Γ) → CBM(K) | 1.772 eV | 1.838 eV |
| Γ valley above CBM | 0.935 eV | 0.918 eV |
| in-plane CB mass at K, k_z = 0 (K→Γ / K→M) | 0.85 / 0.68 | 0.83 / 0.66 |
| k_F (Luttinger, 2 valleys) | 0.236 Å⁻¹ | 0.236 Å⁻¹ |

![x010](figures/x0.10_bands_dos.png)
![fermi](figures/fermi_surface_x0.10.png)

The Fermi contours of the two descriptions lie on top of each other: two slightly trigonal pockets of radius
≈ 0.24 Å⁻¹ centred on K and K′. For the occupied states, the rigid band is essentially exact at this doping. What it
misses is the ≈ 65 meV relative shift of the valence band.

## 3. The full intercalate ZrNCl{Co(Cp)₂}₀.₁₀

**Model.**
- Composition: one ZrNCl double layer (10 formula units) plus one Co(Cp)₂ per gallery (51 atoms), i.e. exactly the
  analysed composition ZrNCl{Co(Cp)₂}₀.₁₀.
- Cell: basal spacing fixed at the measured 14.7 Å; in-plane √3a × √7a supercell giving 56.3 Å² per guest.
- Guest orientation: Cp–Co–Cp axis parallel to the layers, as the paper infers from the Cp/Cp′/Cp* series. Along the
  6.24 Å cell vector the guests interdigitate ring-edge to ring-edge.
- Relaxation: PBE+D3(BJ), all atoms, residual force 0.021 eV/Å.

**Relaxed geometry** ([figure](figures/structures.png)):

| ZrNCl layer | pristine (PBE+D3) | jellium x = 0.10 (PBE) | in ZrNCl{Co(Cp)₂}₀.₁₀ (PBE+D3) |
|---|---|---|---|
| Zr–N (mean of 4) | 2.136 Å | 2.130 Å | 2.131 Å (2.120–2.157) |
| Zr–Cl (×3) | 2.749 Å | 2.770 Å | 2.780 Å (2.769–2.790) |
| Cl–Cl slab thickness | 6.199 Å | 6.227 Å | 6.264 Å |

GUEST_TABLE

- The guest shrinks from the neutral-cobaltocene starting geometry (Co–C 2.10 Å) to Co–C = 2.058 Å and
  Co–(Cp centroid) = 1.659 Å, which is cobaltocenium-like. See the isolated-molecule references in the table.
- The uppermost and lowermost Cp hydrogens sit in the Cl hollows of the two walls: H···Cl 2.70 Å, Co 4.2 Å above the
  Cl plane.
- The host responds as it does to jellium doping, only more strongly:
  - Zr–Cl lengthens from 2.749 to 2.780 Å (jellium x = 0.10: 2.770 Å).
  - The Cl–Cl slab thickness grows by 0.065 Å.
  - Zr–N is unchanged, 2.131 Å on average.

**Electronic structure and charge transfer** (PBE; SCF on 9×9×1 k with 0.02 eV Fermi–Dirac smearing; tetrahedron
DOS; E_F from electron counting):

![pdos](figures/intercalate_pdos.png)

| quantity (per Co(Cp)₂ unless stated) | PBE, 9×9×1, σ = 0.02 eV | PBE, 6×6×1, σ = 0.05 eV |
|---|---|---|
| electrons in ZrNCl conduction-band states | **0.76** | 0.73 |
| electrons left in the guest e₁″ band | 0.24–0.28 | 0.27 |
| effective doping x_eff (e⁻/ZrNCl) | **0.076** | 0.073 |
| E_F − E_CBM(ZrNCl layer) | 339 meV | 280 meV (smearing-limited) |
| guest e₁″ band relative to E_F | −0.00 … +0.11 eV (pinned at E_F) | +0.05 … +0.18 eV |
| occupied Co 3d (a₁′, e₂′) | −2.9 … −3.2 eV (below the layer VBM at −2.24 eV) | – |
| ZrNCl layer gap (VBM → CBM) | 1.90 eV | – |
| layer-type DOS at the CB bottom | 0.26 states eV⁻¹ ZrNCl⁻¹ | – |
| Bader charge of Co(Cp)₂ (all-electron density) | **+0.64 e** (Co +0.59) | – |
| Bader: Zr / N / Cl (pristine: +2.27 / −1.60 / −0.67) | +2.24 / −1.61 / −0.69 | – |

MAGNETISM_AND_U

**Unfolded band structure.** The supercell states are projected back onto the 1×1 ZrNCl zone. The spectral weights
are clean: 93 % of the layer states have weight > 0.8 or < 0.2.

![unfold](figures/intercalate_unfolded_bands.png)

- **The host bands are almost those of the jellium-doped host** (grey lines). This holds for the valence band, the
  K-valley conduction band (bottom at E_F − 0.34 eV, k_F = 0.21–0.23 Å⁻¹) and the higher Zr-4d bands.
  - The main visible difference is that the Γ valley sits higher: 1.08 eV above the K minimum, as in the jellium host
    after internal relaxation (1.06 eV). The Zr–Cl expansion causes this.
  - Neither the guest's electrostatic potential nor its hybridisation reshapes the Zr-4d band.
- **The guest contributes flat levels.** The occupied Co 3d (a₁′/e₂′) levels lie at −2.9 and −3.2 eV, inside the
  ZrNCl valence-band energy range. The e₁″ band (bandwidth ≈ 0.1 eV, from guest–guest overlap along the 6.24 Å rows)
  lies at E_F.
- **Luttinger count.** The Fermi wave vector corresponds to 0.08–0.09 e⁻/ZrNCl, consistent with the 0.076 e⁻/ZrNCl
  from state counting. Both lie below the 0.10 that complete ionisation of every guest would give.

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
