#!/usr/bin/env python3
"""
Generate crystallographic CIF files from r2scand3_<CSD>.maegz files
for all molecule folders in the csp-validation-data repository.

CIF format includes:
    - Space group (_symmetry_space_group_name_H-M)
    - Unit cell parameters (_cell_length_a/b/c, _cell_angle_alpha/beta/gamma)
    - Z value (_cell_formula_units_Z)
    - Fractional coordinates (_atom_site_fract_x/y/z)
    - CSP annotations (_csp_NAME_EXP_FORM, _csp_CSD_REFCODE, _csp_title)

Usage:
    $SCHRODINGER/run python3 generate_cif.py [--overwrite] [--folder FOLDER]

Options:
    --overwrite      Overwrite existing CIF files (default: skip)
    --folder FOLDER  Only process a specific folder (default: all)
"""

import os
import sys
import re
import glob
import numpy as np
from schrodinger import structure


def get_cell_params(st):
    """Extract cell parameters from lattice vectors."""
    a_vec = np.array([
        st.property['r_chorus_box_ax'],
        st.property['r_chorus_box_ay'],
        st.property['r_chorus_box_az']
    ])
    b_vec = np.array([
        st.property['r_chorus_box_bx'],
        st.property['r_chorus_box_by'],
        st.property['r_chorus_box_bz']
    ])
    c_vec = np.array([
        st.property['r_chorus_box_cx'],
        st.property['r_chorus_box_cy'],
        st.property['r_chorus_box_cz']
    ])

    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)
    alpha = np.degrees(
        np.arccos(np.clip(np.dot(b_vec, c_vec) / (b * c), -1, 1)))
    beta = np.degrees(
        np.arccos(np.clip(np.dot(a_vec, c_vec) / (a * c), -1, 1)))
    gamma = np.degrees(
        np.arccos(np.clip(np.dot(a_vec, b_vec) / (a * b), -1, 1)))

    lat = np.array([a_vec, b_vec, c_vec])
    return a, b, c, alpha, beta, gamma, lat


def cart_to_frac(st, lat):
    """Convert Cartesian coordinates to fractional."""
    coords = np.array([[atom.x, atom.y, atom.z] for atom in st.atom])
    return coords @ np.linalg.inv(lat)


def extract_z(formula):
    """Extract Z from unit cell formula string, e.g. '(C9H8O4)4' -> 4."""
    m = re.search(r'\)(\d+)', formula)
    return m.group(1) if m else '?'


def write_cif(st, cif_path):
    """Write a proper crystallographic CIF file."""
    a, b, c, alpha, beta, gamma, lat = get_cell_params(st)
    frac = cart_to_frac(st, lat)

    sg = st.property.get('s_pdb_PDB_CRYST1_Space_Group', 'P 1').strip()
    title = st.property.get('s_m_title', '').strip()
    formula = st.property.get('s_m_Unit_Cell_Formula', '').strip()
    form = st.property.get('s_csp_NAME_EXP_FORM', 'none')
    ref = st.property.get('s_csp_CSD_REFCODE', 'none')
    z_val = extract_z(formula)

    data_name = os.path.basename(cif_path).replace('.cif', '')

    lines = []
    lines.append('data_' + data_name)
    lines.append('')
    lines.append('_symmetry_space_group_name_H-M   "' + sg + '"')
    lines.append('_cell_length_a   ' + f'{a:.6f}')
    lines.append('_cell_length_b   ' + f'{b:.6f}')
    lines.append('_cell_length_c   ' + f'{c:.6f}')
    lines.append('_cell_angle_alpha   ' + f'{alpha:.4f}')
    lines.append('_cell_angle_beta   ' + f'{beta:.4f}')
    lines.append('_cell_angle_gamma   ' + f'{gamma:.4f}')
    lines.append('_cell_formula_units_Z   ' + z_val)
    lines.append('')
    lines.append('_csp_NAME_EXP_FORM   ' + form)
    lines.append('_csp_CSD_REFCODE   ' + ref)
    lines.append('_csp_title   "' + title + '"')
    lines.append('')
    lines.append('loop_')
    lines.append('_atom_site_label')
    lines.append('_atom_site_type_symbol')
    lines.append('_atom_site_fract_x')
    lines.append('_atom_site_fract_y')
    lines.append('_atom_site_fract_z')

    elem_count = {}
    for i, atom in enumerate(st.atom):
        elem = atom.element
        elem_count[elem] = elem_count.get(elem, 0) + 1
        label = elem + str(elem_count[elem])
        fx, fy, fz = frac[i]
        lines.append(
            label.rjust(6) + ' ' + elem.rjust(2) + '   '
            + f'{fx:12.8f}' + '  '
            + f'{fy:12.8f}' + '  '
            + f'{fz:12.8f}')

    lines.append('')

    with open(cif_path, 'w') as f:
        f.write('\n'.join(lines))


def extract_csd_code(folder_name):
    """Extract CSD code from folder name (e.g., 'Asprin_ACSALA' -> 'ACSALA')."""
    parts = folder_name.rsplit('_', 1)
    return parts[-1] if len(parts) > 1 else folder_name


def process_folder(folder_path, csd, overwrite=False):
    """Process a single folder: generate CIFs from maegz."""
    maegz_path = os.path.join(folder_path, f'r2scand3_{csd}.maegz')
    if not os.path.exists(maegz_path):
        print(f'MISSING maegz: {folder_path}')
        return 0

    structs = list(structure.StructureReader(maegz_path))

    if overwrite:
        old_cifs = glob.glob(
            os.path.join(folder_path, f'r2scand3_{csd}_*.cif'))
        for f in old_cifs:
            os.remove(f)

    created = 0
    for i, st in enumerate(structs):
        rank = i + 1
        cif_path = os.path.join(
            folder_path, f'r2scand3_{csd}_{rank:02d}.cif')

        if os.path.exists(cif_path) and not overwrite:
            continue

        write_cif(st, cif_path)
        created += 1

    return created


def main():
    overwrite = '--overwrite' in sys.argv
    target_folder = None
    if '--folder' in sys.argv:
        idx = sys.argv.index('--folder')
        if idx + 1 < len(sys.argv):
            target_folder = sys.argv[idx + 1]

    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    base = os.path.normpath(base)

    folders = sorted(os.listdir(base))
    folders = [f for f in folders
               if os.path.isdir(os.path.join(base, f))]

    if target_folder:
        folders = [f for f in folders if f == target_folder]
        if not folders:
            print(f'Folder not found: {target_folder}')
            return

    total_created = 0
    total_folders = 0

    for folder in folders:
        csd = extract_csd_code(folder)
        folder_path = os.path.join(base, folder)
        created = process_folder(folder_path, csd, overwrite)
        if created > 0:
            print(f'CREATED: {created} CIF files in {folder}')
        else:
            print(f'SKIP: {folder}')
        total_created += created
        total_folders += 1

    print(f'\nSummary: {total_created} CIF files created '
          f'across {total_folders} folders')


if __name__ == '__main__':
    main()
