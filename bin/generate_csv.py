#!/usr/bin/env python3
"""
Generate r2scand3_<CSD>.csv files from r2scand3_<CSD>.maegz files
for all molecule folders in the csp-validation-data repository.

CSV format:
    pbed3_density, r2scand3_relative_energy, space_group, title

Usage:
    $SCHRODINGER/run python3 generate_csv.py [--overwrite]

Options:
    --overwrite  Overwrite existing CSV files (default: skip)
"""

import os
import sys
from schrodinger import structure


def extract_csd_code(folder_name):
    """Extract CSD code from folder name (e.g., 'Asprin_ACSALA' -> 'ACSALA')."""
    parts = folder_name.rsplit('_', 1)
    return parts[-1] if len(parts) > 1 else folder_name


def generate_csv(maegz_path, csv_path):
    """Read a maegz file and write the corresponding CSV."""
    reader = structure.StructureReader(maegz_path)
    structs = list(reader)

    lines = ['pbed3_density,r2scand3_relative_energy,space_group,title,NAME_EXP_FORM,CSD_REFCODE']
    for st in structs:
        density = st.property.get('r_m_Unit_Cell_Density/g/cm^3', 0.0)
        rel_energy = st.property.get(
            'r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)', 0.0)
        space_group = st.property.get(
            's_pdb_PDB_CRYST1_Space_Group', '').strip()
        title = st.property.get('s_m_title', '').strip()
        form = st.property.get('s_csp_NAME_EXP_FORM', 'none')
        ref = st.property.get('s_csp_CSD_REFCODE', 'none')
        lines.append(
            f'{density:>10.5f},{rel_energy:>11.5f},{space_group:>11},{title:>30},{form},{ref}')

    with open(csv_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    return len(structs)


def main():
    overwrite = '--overwrite' in sys.argv
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    base = os.path.normpath(base)

    folders = sorted(os.listdir(base))
    folders = [f for f in folders if os.path.isdir(os.path.join(base, f))]

    created = 0
    skipped = 0
    errors = 0

    for folder in folders:
        csd = extract_csd_code(folder)
        maegz = os.path.join(base, folder, f'r2scand3_{csd}.maegz')
        csv_path = os.path.join(base, folder, f'r2scand3_{csd}.csv')

        if not os.path.exists(maegz):
            print(f'MISSING maegz: {folder}')
            errors += 1
            continue

        if os.path.exists(csv_path) and not overwrite:
            print(f'SKIP (CSV exists): {folder}')
            skipped += 1
            continue

        n = generate_csv(maegz, csv_path)
        action = 'OVERWRITTEN' if os.path.exists(csv_path) else 'CREATED'
        print(f'{action}: r2scand3_{csd}.csv ({n} structures) in {folder}')
        created += 1

    print(f'\nSummary: {created} created, {skipped} skipped, {errors} errors')


if __name__ == '__main__':
    main()
