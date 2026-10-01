# Computational details

## Code and physics
- **GPAW 24.1** (projector-augmented waves, plane-wave basis) driven by **ASE 3.22**, standard GPAW PAW setups 0.9.20000
  (Zr: 4s²4p⁶4d²5s², 12 e⁻; N 5; Cl 7; Co 3d⁷4s², 9; C 4; H 1).
- Exchange–correlation: **PBE**. Dispersion: **Grimme D3(BJ)** (`simple-dftd3`) added for every geometry relaxation.
  D3 changes only the geometry, not the Kohn–Sham Hamiltonian.
- Plane-wave cutoff **600 eV**. Relative to 800 eV, forces agree within 0.01 eV/Å and the gap within 2 meV
  (`figures/convergence_cutoff.png`, `results/convergence.json`).
- Occupations: Fermi–Dirac smearing, 0.01 eV (insulator), 0.02 eV (doped host and intercalate electronic structure),
  0.1 eV (intercalate relaxation only).
  - DOS and N(E_F) use the linear tetrahedron method, either on non-self-consistent dense meshes or on the SCF mesh.
  - E_F is always obtained by counting electrons with the tetrahedron DOS, integrated up from a clean gap below E_F.
  - For the intercalate, 0.05 eV smearing is too coarse: it put ≈ 0.2 e too many into the guest level, which lies
    within 0.1–0.2 eV of E_F (`results/intercalate_electronic_k6_w0.05.json`).
- BLAS: OpenBLAS. The Ubuntu image linked GPAW to reference BLAS at first; everything after the cutoff scan used OpenBLAS.

## Structures
- **β-ZrNCl**: SmSI type, R-3m, all atoms on 6c (0,0,z). The experimental lattice a = 3.6046 Å, c = 27.672 Å is kept fixed.
  Only the internal z parameters are relaxed, with FixSymmetry.
  - Cell choice: compact primitive cell a1, a2 (in-plane) and a3 = (2/3,1/3,1/3)ₕₑₓ. This is the same lattice as the
    α = 22° rhombohedral cell, but the FFT grid is 7× smaller (21×21×54 instead of 54×54×54 at 600 eV).
  - Point-group symmetry is switched off and only time-reversal is kept; GPAW would otherwise symmetrise the grid back to
    N₁ = N₂ = N₃.
  - As a lattice check, a PBE+D3 variable-cell relaxation is done at 800 eV to suppress Pulay stress.
- **Electron-doped host (jellium doping)**: the same cell with GPAW `charge = −2x` per Zr₂N₂Cl₂.
  The extra x electrons per ZrNCl are compensated by a uniform positive background (the G = 0 term of the Hartree
  potential is dropped).
  - "Strict" rigid band: the pristine tetrahedron DOS is integrated upward from the CBM.
- **ZrNCl{Co(Cp)₂}₀.₁₀**: an ordered model of the experimental stoichiometry.
  - Cell: one ZrNCl double layer (Cl–Zr–N–N–Zr–Cl, 10 ZrNCl) plus one Co(Cp)₂ per gallery, 51 atoms.
  - In-plane supercell A₁ = 2a₁ + a₂ (|A₁| = √3a = 6.24 Å), A₂ = −a₁ + 2a₂ (√7a = 9.54 Å): 56.3 Å² per guest, i.e.
    exactly x = 1/10.
  - Basal spacing fixed at the XRD value d = 14.7 Å, layers AA-stacked.
  - The Cp–Co–Cp axis lies parallel to the layers, as deduced in Fogg et al. from d(CoCp₂) = d(CoCp′₂) < d(CoCp*₂).
  - Guests along the 6.24 Å vector meet ring-edge to ring-edge with interdigitated pentagons (shortest guest–guest
    H···H 2.48 Å). The uppermost and lowermost H of each ring point into Cl hollows (H···Cl 2.77 Å at the start).
  - Relaxation: all atoms relaxed with the cell fixed, using plane waves (600 eV) + D3(BJ) and ASE BFGS (maximum
    step 0.1 Å) to a residual force of 0.021 eV/Å.
    - The k mesh is a Γ-centred 3×3×1, which contains the folded K/K′ points (0, ±1/3) where the host conduction band
      has its minimum.
    - The relaxation keeps the Cm mirror plane of the starting model.
    - Abandoned approaches: an LCAO/dzp pre-relaxation (slower per SCF than plane waves) and PreconLBFGS (its line
      search cost 2–3 SCFs per step).
  - `03_intercalate/symmetry_check.py` tests whether the axis-parallel orientation is a minimum.
    - It restarts from the relaxed structure with the guest rotated rigidly, by 8° out of the layer plane and 6° about
      the layer normal, which breaks the mirror.
    - It then relaxes with the same settings (BFGS, fmax 0.05 eV/Å, at most 40 steps).
    - `analysis/symcheck_progress.py` follows the energy and the axis orientation step by step.
- **Isolated Co(Cp)₂ (S = ½) and Co(Cp)₂⁺ (S = 0)**: real-space PAW (h = 0.18 Å, 5.5 Å vacuum) with open boundaries,
  so the cation needs no charged-cell correction.
  - These calculations give the Co–C fingerprints of the oxidation state, the ΔSCF ionisation energies and PBE's
    curvature for the e₁″ level, ε_SOMO(neutral) − ε_LUMO(cation) at one geometry.
  - Fermi–Dirac smearing of 0.1 eV keeps the neutral's singly occupied e₁″ pair half/half filled, i.e. the
    D₅h-averaged molecule.
  - The neutral relaxation was stopped at fmax 0.05 eV/Å. Beyond that point the geometry starts to follow the
    Jahn–Teller distortion and the SCF jumps between orbitally polarised solutions.
  - The cation was relaxed to 0.02 eV/Å.

## k-point meshes (Γ-centred)
| system | SCF | DOS / N(E_F) | bands |
|---|---|---|---|
| β-ZrNCl, relaxation | 8×8×2 | – | – |
| β-ZrNCl | 12×12×3 | 30×30×3 (tetrahedron) | Γ–M–K–Γ–Z (kz = 0 in-plane path) |
| doped series | 18×18×2 (check 24×24×2) | tetrahedron on SCF mesh | – |
| doped x = 0.10 | 18×18×3 | 30×30×3 | Γ–M–K–Γ–Z |
| intercalate, relaxation | 3×3×1 | – | – |
| intercalate | 9×9×1 (σ = 0.02 eV) | tetrahedron on the 9×9×1 SCF mesh | unfolded onto Γ–M–K–Γ of the 1×1 cell (56 k) |
| intercalate checks (+U, guest shift, spin, fragments) | 6×6×1 (σ = 0.02 eV) | tetrahedron on the SCF mesh | – |

All electronic-structure meshes contain the in-plane K point, the conduction-band minimum.

k-point convergence (pristine, 600 eV, `01_pristine/convergence_k.py`, `results/convergence.json`):
- The total energy differs by ≤ 0.1 meV/atom between the 6×6×2, 8×8×2 and 12×12×3 meshes.
- The forces differ by < 1 meV/Å.
- The 8×8×2 relaxation mesh does not contain K. That matters for the gap (1.90 instead of 1.78 eV), not for the energy
  or forces of the insulator.

## Analysis
- **Charge transfer in the intercalate**, four independent measures:
  1. The number of electrons in host-type conduction-band states: Kohn–Sham states whose PAW-projection weight lies
     > 50 % on the layer, above the middle of the layer gap.
  2. The positions of the guest levels relative to E_F: the occupied Co 3d "a₁′/e₂′" set and the empty e₁″* level,
     which is the SOMO of neutral Co(Cp)₂.
  3. Bader charges from the all-electron density (`gridrefinement=2`, `pybader`). Hirshfeld partitioning is not
     available for GPAW plane-wave mode.
  4. The plane-averaged density difference Δρ(z) = ρ[intercalate] − ρ[slab] − ρ[guest] at frozen geometry, and its
     integral (the charge-displacement curve).
- **Band unfolding**: spectral weights P_Km(k) = Σ_g |C_Km(k−K+g)|² from the pseudo-wavefunction plane-wave
  coefficients, using the 1×1 ZrNCl cell of the same layer.
  - The run is non-self-consistent from the 9×9×1 density, along Γ–M–K–Γ with 20 points per Å⁻¹, using 180 bands of
    which the lowest 168 are converged to 10⁻⁴ eV².
  - The supercell matrix is M = [[2,1,0],[−1,2,0],[0,0,1]]. GPAW's G vectors are in bohr⁻¹, so the cell enters in bohr.
  - Each state is labelled layer or guest by its PAW-projection weight. Fermi wave vectors are taken where the
    unfolded layer band crosses E_F.
- **Sensitivity checks of the charge split** (6×6×1, σ = 0.02 eV, same geometry), compared with plain PBE on the
  same mesh:
  - PBE+U: Dudarev U_eff = 4 eV on Co 3d.
  - Guest-level shift: a smooth potential step of +0.5 eV for electrons on the guest. It is applied on the union of
    spheres R = 1.8 Å (erfc edge 0.25 Å) around Co and the ten C atoms; the Cl nuclei lie ≥ 3.5 Å away. This is a
    scissor-like stand-in for a correction of the guest's frontier level, not a functional.
  - Spin polarisation: the SCF is started from 1 μ_B on Co.
- Conduction-band effective masses: parabolic fits within |k − K| ≤ 0.05 Å⁻¹ along K→Γ and K→M.
