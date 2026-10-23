# CSP Validation Data

Predicted crystal structure polymorphs for 66 organic molecules (1,273 structures) generated using Schrödinger Crystal Structure Prediction (CSP). Structures are rank-ordered by r2SCAN-D3 energy within a 3 kcal/mol energy window.

## Citation

D. Zhou, I. Bier, B. Santra, L. D. Jacobson, C. Wu, A. Garaizar Suarez, B. Ramirez Almaguer, H. Yu, R. Abel, R. A. Friesner, L. Wang, "A robust crystal structure prediction method to support small molecule drug development with large scale validation and blind study", *Nature Communications*, **16**, 2210 (2025). DOI: [10.1038/s41467-025-57479-1](https://doi.org/10.1038/s41467-025-57479-1)

## Repository Structure

```
csp-validation-data/
├── README.md
├── LICENSE
├── Sysinfo.csv                                  # CSD refcodes and common names for all 66 systems
├── best_match_str_manuscript_all_systems.csv     # Best-match structures to 143 experimental forms
├── bin/                                          # Utility scripts
│   ├── generate_csv.py                           # Generate CSV summaries from .maegz files
│   ├── generate_cif.py                           # Generate crystallographic CIF files from .maegz
│   └── primitive_to_conventional.py              # Convert primitive cells to conventional cells
├── data/                                         # Crystal structure data (66 folders)
│   └── <SystemName>_<CSD>/                       # e.g., ROY_QAXMEH, Asprin_ACSALA
│       ├── r2scand3_<CSD>.maegz                  # All polymorphs in Maestro format
│       ├── r2scand3_<CSD>.csv                    # Summary table
│       └── r2scand3_<CSD>_XX.cif                 # Individual CIF files (XX = rank)
└── docs/                                         # Publication and supplementary materials
    ├── s41467-025-57479-1.pdf
    ├── 41467_2025_57479_MOESM1_ESM.pdf
    └── 41467_2025_57479_MOESM4_ESM.xlsx
```

## Data Description

### Maestro files (.maegz)

Each `.maegz` file contains all predicted polymorphs for a given molecule, sorted by `r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)`. Key structure-level properties include:

| Property | Description |
|---|---|
| `s_m_title` | Structure identifier with space group, CSD refcode, and form name |
| `r_matsci_Etot_(Ry)` | Total DFT energy in Rydberg; may refer to the original calculation cell (see note below) |
| `r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)` | Relative energy per molecule in kcal/mol, independent of cell representation |
| `r_m_Unit_Cell_Density/g/cm^3` | Unit cell density |
| `s_pdb_PDB_CRYST1_Space_Group` | Space group |
| `s_m_Unit_Cell_Formula` | Unit cell formula with Z value |
| `s_csp_NAME_EXP_FORM` | Matched experimental form name (e.g., form-I, form-alpha) or `none` |
| `s_csp_CSD_REFCODE` | Matched CSD refcode or `none` |

**Energy and cell representation:** Many DFT calculations were performed using primitive cells, while the distributed structures may use conventional cells. Consequently, `r_matsci_Etot_(Ry)` does not always correspond to the conventional cell currently stored in a structure; it may refer to the original DFT calculation cell. Do not compare these total energies directly across structures without checking the calculation cell and the number of molecules it contains. The normalized `r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)` remains valid for either primitive or conventional cell representations and should be used to compare polymorph energies within each molecular system. This is also the quantity exported as `r2scand3_relative_energy` in the CSV summaries.

### CIF files

Crystallographic CIF files with unit cell parameters, space group, Z value, and fractional coordinates. Each file includes the publication citation.

### CSV files

Summary tables with columns: `pbed3_density`, `r2scand3_relative_energy`, `space_group`, `title`, `NAME_EXP_FORM`, `CSD_REFCODE`.

### Experimental form matching

The file `best_match_str_manuscript_all_systems.csv` contains the best-match predicted structures for 143 experimentally known polymorphic forms across 66 molecules, including RMSD values and energy rankings.

## Utility Scripts

All scripts require the Schrödinger Suite Python environment:

```bash
$SCHRODINGER/run python3 bin/<script>.py [options]
```

| Script | Description |
|---|---|
| `generate_csv.py` | Generate CSV summaries from `.maegz` files. Use `--overwrite` to regenerate. |
| `generate_cif.py` | Generate crystallographic CIF files from `.maegz` files. Use `--overwrite` to regenerate. Use `--folder <name>` for a specific system. |
| `primitive_to_conventional.py` | Convert primitive cell structures to conventional cells using `schrodinger.comparison.struc.get_conventional_unit_cell`. Use `--symprec` to adjust symmetry tolerance. |

## DFT Method

All structures were optimized using r2SCAN-D3 density functional with tight convergence settings via Quantum ESPRESSO through the Schrödinger periodic DFT interface.
