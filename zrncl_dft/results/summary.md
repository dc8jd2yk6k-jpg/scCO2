### Pristine β-ZrNCl: structure (R-3m, 6c sites (0,0,z))

| | a (Å) | c (Å) | z(Zr) | z(N) | z(Cl) | Zr–N ×3 (Å) | Zr–N apical (Å) | Zr–Cl ×3 (Å) |
|---|---|---|---|---|---|---|---|---|
| starting model (exp. cell) | 3.6046 | 27.672 | 0.1196 | 0.1981 | 0.3888 | 2.126 | 2.172 | 2.735 |
| PBE, exp. cell | 3.6046 | 27.672 | 0.1195 | 0.1976 | 0.3878 | 2.129 | 2.164 | 2.751 |
| PBE+D3(BJ), exp. cell | 3.6046 | 27.672 | 0.1196 | 0.1976 | 0.3880 | 2.129 | 2.159 | 2.749 |
| PBE+D3(BJ), cell relaxed (800 eV) | 3.5919 | 27.405 | 0.1190 | 0.1979 | 0.3868 | 2.122 | 2.161 | 2.744 |

### Pristine β-ZrNCl: electronic structure (PBE)

| quantity | value |
|---|---|
| band gap (indirect, dense mesh) | 1.772 eV |
| smallest direct gap | 2.109 eV |
| CB curvature mass at K, k_z = 0 (K→Γ / K→M) | 0.854 / 0.682 m_e |
| CB density-of-states mass (k_z-averaged) | 0.583 m_e |
| CB kz dispersion at K (E(K,Z) − E(K,0)) | -0.1 meV |
| VB kz dispersion at Γ (E(Z) − E(Γ)) | -256.8 meV |

### Electron doping: strict rigid band vs extra electrons + jellium (SCF)

| x (e/ZrNCl) | E_F−E_CBM rigid (meV) | E_F−E_CBM SCF (meV) | N(E_F) rigid | N(E_F) SCF | k_F 2D (1/Å) |
|---|---|---|---|---|---|
| 0.025 | 100 | 112 | 0.283 | 0.239 | 0.118 |
| 0.050 | 192 | 209 | 0.282 | 0.278 | 0.167 |
| 0.075 | 283 | — | 0.275 | — | — |
| 0.100 | 371 | 390 | 0.289 | 0.276 | 0.236 |
| 0.125 | 457 | — | 0.291 | — | — |
| 0.150 | 542 | 567 | 0.302 | 0.302 | 0.289 |
| 0.200 | 702 | 731 | 0.324 | 0.330 | 0.334 |
| 0.250 | 849 | — | 0.361 | — | — |
| 0.300 | 954 | 935 | 1.127 | 1.251 | 0.409 |
| 0.400 | 1032 | 996 | 1.520 | 1.568 | 0.473 |
| 0.500 | 1093 | — | 1.676 | — | — |

N(E_F) in states/eV/ZrNCl (both spins).  CB-bottom DOS 0.274 states/eV/ZrNCl ⇒ m*_DOS = 0.58 m_e (k_z = 0 in-plane curvature: 0.77 m_e).

### Isolated Co(Cp)₂ and Co(Cp)₂⁺ (PBE, FD, open boundaries)

| | Co–C (Å) | Co–Cp centroid (Å) | C–C (Å) |
|---|---|---|---|
| Co(Cp)₂ (S=½) | 2.102 | 1.713 | 1.432 |
| Co(Cp)₂⁺ (S=0) | 2.045 | 1.642 | 1.433 |
| in ZrNCl{Co(Cp)₂}₀.₁₀ | 2.058 | 1.659 | 1.431 |

Ionisation energy (ΔSCF): adiabatic 5.01 eV, vertical 5.13 eV.

### ZrNCl{Co(Cp)₂}₀.₁₀ model: electronic structure and charge transfer

| quantity | value |
|---|---|
| E_F − E_CBM(ZrNCl layer) | 339 meV |
| ZrNCl layer gap (VBM→CBM) | 1.905 |
| electrons in ZrNCl conduction band (per Co(Cp)₂) | 0.761 |
| effective doping x (e/ZrNCl) | 0.076 |
| N(E_F) (states/eV/ZrNCl) | 13.825 |
| Bader charge of Co(Cp)₂ | +0.640 e |
| highest occupied guest levels (E−E_F, eV) | -2.87, -2.86, -2.85, -0.01 |
| lowest empty guest levels (E−E_F, eV) | 0.00, 0.01, 0.02, 0.04 |

### Charge-split checks (6×6×1 k, Fermi–Dirac 0.02 eV, same geometry)

| | PBE | PBE+U (U_eff = 4 eV, Co 3d) | PBE, guest levels +0.5 eV |
|---|---|---|---|
| electrons in ZrNCl CB states, per guest | 0.78 | 0.80 | 0.82 |
| electrons in guest states near E_F | 0.25 | 0.16 | 0.09 |
| x_eff (e⁻/ZrNCl) | 0.078 | 0.080 | 0.082 |
| E_F − E_CBM(layer) (meV) | 292 | 321 | 333 |
| guest states within 1 eV of E_F (eV) | -0.00 … +0.12 | +0.01 … +0.21 | +0.01 … +0.14 |
| top of the occupied Co 3d (a₁′/e₂′) levels (eV) | -2.84 | -3.97 | -2.82 |
| layer gap (eV) | 1.91 | 1.91 | 1.91 |
| N(E_F), layer states (states/eV/ZrNCl) | 0.266 | 0.261 | 0.268 |

Spin-polarised SCF (6×6×1, started from 1 μB on Co): total moment -0.000 μB, Co -0.000 μB, Σ|m| 0.001 μB; E(spin) − E(non-spin) = -0.0 meV per cell.

Isolated molecule: ε_SOMO(Co(Cp)₂) − ε_LUMO(Co(Cp)₂⁺) at the neutral geometry = 4.91 eV (PBE curvature of E(N) for the e₁″ level; zero for the exact functional).

Symmetry check (guest rotated 8° out of plane + 6° in plane, relaxed 28 BFGS steps, fmax 0.046 eV/Å): E − E(Cm) = -4.7 meV, final axis tilt 0.95°, Co–C 2.0586 Å.
