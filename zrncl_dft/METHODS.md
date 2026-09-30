# Computational details

## Code and physics
- **GPAW 24.1** (projector-augmented waves, plane-wave basis) driven by **ASE 3.22**, standard GPAW PAW setups 0.9.20000
  (Zr: 4s²4p⁶4d²5s², 12 e⁻; N 5; Cl 7; Co 3d⁷4s², 9; C 4; H 1).
- Exchange–correlation: **PBE**. Dispersion: **Grimme D3(BJ)** (`simple-dftd3`) added for every geometry relaxation.
  D3 changes only the geometry, not the Kohn–Sham Hamiltonian.
- Plane-wave cutoff **600 eV**. Relative to 800 eV, forces agree within 0.01 eV/Å and the gap within 2 meV
  (`figures/convergence_cutoff.png`, `results/convergence.json`).
- Occupations: Fermi–Dirac smearing, 0.01 eV (insulator), 0.02 eV (doped host) and 0.05–0.1 eV (intercalate).
  DOS and N(E_F) use the linear tetrahedron method on non-self-consistent dense meshes.
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
  - Relaxation: all atoms relaxed with the cell fixed. First LCAO/dzp+D3 (fmax 0.05 eV/Å), then plane waves+D3
    (fmax 0.02–0.03 eV/Å), both on a Γ-centred 3×3×1 mesh. That mesh contains the folded K/K′ points (0, ±1/3) where
    the host conduction band has its minimum.
- **Isolated Co(Cp)₂ (S = ½) and Co(Cp)₂⁺ (S = 0)**: real-space PAW (h = 0.18 Å) with open boundaries, so the cation
  needs no charged-cell correction. Gives Co–C fingerprints of the oxidation state and the ΔSCF ionisation energy.

## k-point meshes (Γ-centred)
| system | SCF | DOS / N(E_F) | bands |
|---|---|---|---|
| β-ZrNCl, relaxation | 8×8×2 | – | – |
| β-ZrNCl | 12×12×3 | 30×30×3 (tetrahedron) | Γ–M–K–Γ–Z (kz = 0 in-plane path) |
| doped series | 18×18×2 (check 24×24×2) | tetrahedron on SCF mesh | – |
| doped x = 0.10 | 18×18×3 | 30×30×3 | Γ–M–K–Γ–Z |
| intercalate | 6×6×1 | 9×9×1 (tetrahedron) | unfolded onto Γ–M–K–Γ of the 1×1 cell |

All meshes contain the in-plane K point, the conduction-band minimum.

## Analysis
- **Charge transfer in the intercalate**, four independent measures:
  1. The number of electrons in host-type conduction-band states: Kohn–Sham states whose PAW-projection weight lies
     > 50 % on the layer, above the middle of the layer gap.
  2. The positions of the guest levels relative to E_F: the occupied Co 3d "a₁′/e₂′" set and the empty e₁″* level,
     which is the SOMO of neutral Co(Cp)₂.
  3. Bader (all-electron density, `pybader`) and Hirshfeld charges.
  4. The plane-averaged density difference Δρ(z) = ρ[intercalate] − ρ[slab] − ρ[guest] at frozen geometry, and its
     integral (the charge-displacement curve).
- **Band unfolding**: spectral weights P_Km(k) = Σ_g |C_Km(k−K+g)|² from the pseudo-wavefunction plane-wave
  coefficients, using the 1×1 ZrNCl cell of the same layer.
- Conduction-band effective masses: parabolic fits within |k − K| ≤ 0.05 Å⁻¹ along K→Γ and K→M.
