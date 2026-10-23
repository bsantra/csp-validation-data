#!/usr/bin/env python3
"""
Convert primitive cell structures to conventional cells in r2scand3_<CSD>.maegz
files. Structures with '_primitive' in the title are converted using
schrodinger.comparison.struc.get_conventional_unit_cell, which preserves
bonding, atom ordering, atom properties, and uses the correct Maestro lattice
convention.

Usage:
    $SCHRODINGER/run python3 primitive_to_conventional.py [--folder FOLDER] [--symprec SYMPREC]

Options:
    --folder FOLDER   Only process a specific folder (default: all)
    --symprec SYMPREC Symmetry precision for spglib (default: 0.01)
"""

import os
import sys
import re
import numpy as np
import spglib
from schrodinger import structure
from schrodinger.comparison.struc import get_conventional_unit_cell


def extract_csd_code(folder_name):
    """Extract CSD code from folder name (e.g., 'Asprin_ACSALA' -> 'ACSALA')."""
    parts = folder_name.rsplit('_', 1)
    return parts[-1] if len(parts) > 1 else folder_name


def update_title(st):
    """Replace _primitive with _conventional in title."""
    title = st.property.get('s_m_title', '')
    if '_primitive' in title:
        st.property['s_m_title'] = title.replace('_primitive', '_conventional')


def update_cell_properties(conv_st, original_st):
    """Refresh cell volume and scale formula Z and total cell energy."""
    lattice = np.array([
        [conv_st.property['r_chorus_box_' + axis + component]
         for component in 'xyz']
        for axis in 'abc'
    ])
    conv_st.property['r_m_Unit_Cell_Volume/Ang.^3'] = float(
        abs(np.linalg.det(lattice)))

    old_natoms = original_st.atom_total
    new_natoms = conv_st.atom_total
    multiplier = new_natoms // old_natoms

    if multiplier <= 1:
        return

    # Update formula Z
    old_formula = original_st.property.get('s_m_Unit_Cell_Formula', '')
    z_match = re.search(r'\)(\d+)', old_formula)
    if z_match:
        old_z = int(z_match.group(1))
        new_z = old_z * multiplier
        new_formula = old_formula[:old_formula.rfind(')') + 1] + str(new_z)
        conv_st.property['s_m_Unit_Cell_Formula'] = new_formula

    # Scale total energy by the multiplier
    if 'r_matsci_Etot_(Ry)' in conv_st.property:
        old_etot = conv_st.property['r_matsci_Etot_(Ry)']
        conv_st.property['r_matsci_Etot_(Ry)'] = old_etot * multiplier


def is_primitive(st, symprec=0.01):
    """Return True if the conventional standard cell is larger than the input
    cell, i.e. the input is a primitive cell of a centered lattice.

    Detection is done with spglib rather than by matching the space-group
    prefix in the title. Prefix matching misses centerings that the title does
    not advertise (C2, R-3, Cmce, ...) and trusts '_primitive'/'_conventional'
    tags that may be wrong in the source data.
    """
    a = [st.property['r_chorus_box_ax'], st.property['r_chorus_box_ay'],
         st.property['r_chorus_box_az']]
    b = [st.property['r_chorus_box_bx'], st.property['r_chorus_box_by'],
         st.property['r_chorus_box_bz']]
    c = [st.property['r_chorus_box_cx'], st.property['r_chorus_box_cy'],
         st.property['r_chorus_box_cz']]
    lat = np.array([a, b, c])
    cart = np.array([[at.x, at.y, at.z] for at in st.atom])
    frac = cart @ np.linalg.inv(lat)
    nums = [at.atomic_number for at in st.atom]
    try:
        ds = spglib.get_symmetry_dataset((lat, frac, nums), symprec=symprec)
    except Exception:
        return False
    if ds is None:
        return False
    std = ds.std_positions if hasattr(ds, 'std_positions') else ds['std_positions']
    return len(std) > st.atom_total


def process_folder(folder_path, csd, symprec=0.01):
    """Process a single folder: convert primitive cells to conventional."""
    maegz_path = os.path.join(folder_path, 'r2scand3_' + csd + '.maegz')
    if not os.path.exists(maegz_path):
        return 0

    structs = list(structure.StructureReader(maegz_path))
    converted = 0
    modified = False

    for i, st in enumerate(structs):
        title = st.property.get('s_m_title', '')
        if not is_primitive(st, symprec=symprec):
            continue

        try:
            conv_st = get_conventional_unit_cell(
                st, symprec=symprec, reassign_bonds=True)

            # Update title
            update_title(conv_st)

            # Update cell volume, formula Z and total energy
            update_cell_properties(conv_st, st)

            conv_sg = conv_st.property.get(
                's_pdb_PDB_CRYST1_Space_Group', '').strip()

            structs[i] = conv_st
            converted += 1
            modified = True
            print('  Rank ' + str(i + 1) + ': ' + title
                  + ' -> ' + conv_sg
                  + ' (' + str(conv_st.atom_total) + ' atoms)')
        except Exception as e:
            print('  ERROR Rank ' + str(i + 1) + ': ' + title
                  + ' -> ' + str(e))

    if modified:
        writer = structure.StructureWriter(maegz_path)
        for st in structs:
            writer.append(st)
        writer.close()

    return converted


def main():
    symprec = 0.01
    target_folder = None

    if '--folder' in sys.argv:
        idx = sys.argv.index('--folder')
        if idx + 1 < len(sys.argv):
            target_folder = sys.argv[idx + 1]

    if '--symprec' in sys.argv:
        idx = sys.argv.index('--symprec')
        if idx + 1 < len(sys.argv):
            symprec = float(sys.argv[idx + 1])

    base = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '..', 'data')
    base = os.path.normpath(base)

    folders = sorted(os.listdir(base))
    folders = [f for f in folders if os.path.isdir(os.path.join(base, f))]

    if target_folder:
        folders = [f for f in folders if f == target_folder]
        if not folders:
            print('Folder not found: ' + target_folder)
            return

    total_converted = 0

    for folder in folders:
        csd = extract_csd_code(folder)
        folder_path = os.path.join(base, folder)
        converted = process_folder(folder_path, csd, symprec)
        if converted > 0:
            print(folder + ': converted ' + str(converted) + ' structures')
        total_converted += converted

    print('\nTotal converted: ' + str(total_converted) + ' structures')


if __name__ == '__main__':
    main()
