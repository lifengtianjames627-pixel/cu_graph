# Landmark structure sources (starting geometries)

Starting `.xyz` files in `ligands/` and `complexes/` are **computational starting guesses** (JT-distorted octahedra / chelate sketches). Use these links as crystallographic / database landmarks when checking bond lengths after ORCA opt.

## Free ligands / ions (PubChem)

| Species | Link |
|---------|------|
| H₂O | https://pubchem.ncbi.nlm.nih.gov/compound/962 |
| Cl⁻ | https://pubchem.ncbi.nlm.nih.gov/compound/312 |
| Br⁻ | https://pubchem.ncbi.nlm.nih.gov/compound/259 |
| OH⁻ | https://pubchem.ncbi.nlm.nih.gov/compound/961 |
| NO₃⁻ | https://pubchem.ncbi.nlm.nih.gov/compound/943 |
| SCN⁻ | https://pubchem.ncbi.nlm.nih.gov/compound/9322 |
| EDTA (H₄EDTA) | https://pubchem.ncbi.nlm.nih.gov/compound/6049 |

## Copper complexes (COD / PubChem)

| Species | Landmark | Link |
|---------|----------|------|
| Hexaaqua Cu(II) salt with nitrate | COD **8102245** `[Cu(H₂O)₆](NO₃)₂·…` | https://www.crystallography.net/cod/8102245.html |
| CuSO₄·5H₂O (chalcanthite; lab salt) | COD **9011339** | https://www.crystallography.net/cod/9011339.html |
| Hexaaqua Cu(II) sulfonate salt | COD **2214985** | https://www.crystallography.net/cod/2214985.html |
| Cu–EDTA | PubChem **23673837** | https://pubchem.ncbi.nlm.nih.gov/compound/23673837 |
| COD search: copper + edta | — | https://www.crystallography.net/cod/result.php?text=copper+edta |

## Notes for this project

- Lab form of EDTA is usually **Na₂H₂EDTA**; DFT uses **EDTA⁴⁻** + `[Cu(EDTA)]²⁻` as the fully deprotonated computational limit (common for thermochem comparisons). Mention this on the poster.
- Thiocyanate on Cu(II) is modeled **N-bound** (`NCS⁻` / isothiocyanate), which is the usual linkage isomer for Cu(II).
- Download CIF from COD if you want to rebuild a crystal fragment later; gas-phase/CPCM cluster models do **not** need the full unit cell.
