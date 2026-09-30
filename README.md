# scCO2

## `zrncl_dft/` — DFT study of Co(Cp)₂-intercalated β-ZrNCl

A plane-wave PAW-DFT (GPAW, PBE / PBE+D3) study that follows Fogg, Green & O'Hare, *Chem. Mater.* **11**, 216 (1999)
("Superconducting metallocene intercalation compounds of β-ZrNCl") in three steps:

1. **pristine β-ZrNCl** (R-3m, SmSI type): relaxation, bands, DOS/PDOS, gap, conduction-band masses;
2. **electron-doped β-ZrNCl**: extra electrons on a jellium background, compared against a strict rigid-band shift
   of the pristine DOS for x = 0.025–0.4 e⁻/ZrNCl;
3. **explicit ZrNCl{Co(Cp)₂}₀.₁₀**: an ordered 51-atom model at the experimental stoichiometry and basal spacing
   (14.7 Å). It is used for charge transfer, guest levels, magnetism and unfolded bands.

Start with [`zrncl_dft/README.md`](zrncl_dft/README.md) for results and discussion.
[`zrncl_dft/METHODS.md`](zrncl_dft/METHODS.md) has the computational details.
