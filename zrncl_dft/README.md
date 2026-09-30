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

RESULTS_PLACEHOLDER
