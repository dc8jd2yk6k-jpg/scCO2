# scCO2

## `zrncl_dft/` — DFT study of Co(Cp)₂-intercalated β-ZrNCl

A plane-wave PAW-DFT (GPAW, PBE / PBE+D3) study that follows Fogg, Green & O'Hare, *Chem. Mater.* **11**, 216 (1999)
("Superconducting metallocene intercalation compounds of β-ZrNCl") in three steps:

1. **pristine β-ZrNCl** (R-3m, SmSI type): relaxation, bands, DOS/PDOS, gap, conduction-band masses;
2. **electron-doped β-ZrNCl**: extra electrons on a jellium background, compared against a strict rigid-band shift
   of the pristine DOS for x = 0.025–0.4 e⁻/ZrNCl;
3. **explicit ZrNCl{Co(Cp)₂}₀.₁₀**: an ordered 51-atom model at the experimental stoichiometry and basal spacing
   (14.7 Å). It is used for charge transfer, guest levels, magnetism and unfolded bands.

**Key results** (PBE; details and caveats in the directory README):
- **Pristine β-ZrNCl** is a band insulator with a 1.77 eV indirect gap (Γ → K).
  - The conduction-band bottom is in-plane Zr 4d at K.
  - At K it is strictly two-dimensional: 0.1 meV of k_z dispersion.
- **Electron doping.** The rigid-band picture gets the filling right: E_F agrees with the self-consistent jellium
  calculation to within 4–8 % up to x = 0.2.
  - N(E_F) ≈ 0.28 states eV⁻¹ ZrNCl⁻¹ and nearly flat in x (+15 % up to x = 0.2), until E_F reaches the Γ valley at
    x ≈ 0.28.
  - The self-consistent charge mainly shifts the valence band, by +65 meV of gap per 0.1 e⁻.
- **ZrNCl{Co(Cp)₂}₀.₁₀.** The unfolded host bands are those of the doped host.
  - Each Co(Cp)₂ gives ≈ 0.8 e to the ZrNCl conduction band. Three independent measures give 0.76–0.87 e: state
    counting, the Fermi wave vector and the guest's Co–C bond length. Bader puts +0.64 e on the guest.
  - The rest stays in the guest's e₁″ level, which the electrostatics of the charged gallery pin at E_F.
  - The guest is non-magnetic in PBE, and the result survives +U on Co and a +0.5 eV shift of the guest levels.
- **What this explains.** A two-dimensional conduction band with a flat N(E_F) explains why the paper's T_c (14 K)
  does not depend on the guest, the gallery height or the doping level. It does not explain the value of T_c.

Start with [`zrncl_dft/README.md`](zrncl_dft/README.md) for results and discussion.
[`zrncl_dft/METHODS.md`](zrncl_dft/METHODS.md) has the computational details.
