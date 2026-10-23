#!/usr/bin/env python3
"""
Detect and repair space groups in every r2scand3_<CSD>.maegz structure,
including cells that do not need expansion. Preserve the existing cell setting
where possible; otherwise adjust its origin or basis to match the assigned
symmetry operations. Primitive centered cells are expanded using Schrodinger's
conventional-cell converter. Stored DFT energies are preserved.

Usage:
    $SCHRODINGER/run python3 primitive_to_conventional.py [--folder FOLDER] [--symprec SYMPREC]

Options:
    --folder FOLDER   Only process a specific folder (default: all)
    --symprec SYMPREC Symmetry precision for spglib (default: 0.01)
"""

import os
import sys
import re
import tempfile
import numpy as np
import spglib
from schrodinger import structure
from schrodinger.application.matsci.nano import space_groups


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
    """Refresh volume and formula Z, preserving original DFT energy metadata."""
    volume = float(abs(np.linalg.det(get_lattice(conv_st))))
    old_volume = conv_st.property.get('r_m_Unit_Cell_Volume/Ang.^3')
    if old_volume is None or not np.isclose(old_volume, volume, rtol=1e-7,
                                           atol=1e-6):
        conv_st.property['r_m_Unit_Cell_Volume/Ang.^3'] = volume

    # Total DFT energy can refer to an earlier calculation cell. Atom count
    # alone cannot establish its normalization, so never rescale it here.
    for key in ('r_matsci_Etot_(Ry)',
                'r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)'):
        if key in original_st.property:
            conv_st.property[key] = original_st.property[key]

    old_natoms = original_st.atom_total
    new_natoms = conv_st.atom_total
    if new_natoms % old_natoms:
        raise ValueError('Conventional-cell atom count is not an integer multiple')
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


def get_lattice(st):
    """Return lattice vectors as rows, in Maestro's Cartesian coordinates."""
    return np.array([[st.property['r_chorus_box_' + axis + component]
                      for component in 'xyz'] for axis in 'abc'])


def get_symmetry_dataset(st, symprec=0.01):
    """Detect symmetry; failures must not silently become P1 assignments."""
    lattice = get_lattice(st)
    fractions = st.getXYZ() @ np.linalg.inv(lattice)
    cell = (lattice, fractions, [atom.atomic_number for atom in st.atom])
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
    if dataset is None:
        raise ValueError('Symmetry detection failed for ' + st.title)
    return dataset


def dataset_value(dataset, key):
    """Support spglib releases returning dictionaries or dataset objects."""
    return dataset[key] if isinstance(dataset, dict) else getattr(dataset, key)


def matching_space_group(dataset, preferred=''):
    """Match complete operations, including translations, in the current basis.

    The dataset's international symbol names the standard setting, which may
    differ from the current cell axes or origin. Do not copy it blindly.
    """
    operations = space_groups.get_symmops_from_spglib({
        key: dataset_value(dataset, key) for key in ('rotations', 'translations')
    })
    candidates = [obj.space_group_short_name for obj in
                  space_groups.get_spacegroups().spg_objs[
                      int(dataset_value(dataset, 'number'))]
                  if space_groups.equal_rotations(obj.symmetry_opers, operations)]
    normalize = lambda name: name.replace(' ', '').replace('_', '')
    if any(normalize(preferred) == normalize(name) for name in candidates):
        return preferred
    return candidates[0] if candidates else None


def set_lattice(st, lattice):
    """Synchronize lattice-vector and CRYST1 properties without idealization."""
    for i, axis in enumerate('abc'):
        for j, component in enumerate('xyz'):
            st.property['r_chorus_box_' + axis + component] = float(lattice[i, j])
    lengths = np.linalg.norm(lattice, axis=1)
    for axis, length in zip('abc', lengths):
        st.property['r_pdb_PDB_CRYST1_' + axis] = float(length)
    for angle, i, j in (('alpha', 1, 2), ('beta', 0, 2), ('gamma', 0, 1)):
        cosine = np.dot(lattice[i], lattice[j]) / (lengths[i] * lengths[j])
        st.property['r_pdb_PDB_CRYST1_' + angle] = float(
            np.degrees(np.arccos(np.clip(cosine, -1, 1))))


def correct_space_group(st, symprec=0.01):
    """Repair labels for every cell, retaining atomic geometry and topology."""
    output = st.copy()
    dataset = get_symmetry_dataset(output, symprec)
    if len(dataset_value(dataset, 'std_positions')) > output.atom_total:
        # Keep the heavier dependency optional for metadata-only repairs.
        from schrodinger.comparison.struc import get_conventional_unit_cell
        output = get_conventional_unit_cell(output, symprec=symprec,
                                             reassign_bonds=True)
        update_title(output)
        dataset = get_symmetry_dataset(output, symprec)

    preferred = st.property.get('s_pdb_PDB_CRYST1_Space_Group', '').strip()
    symbol = matching_space_group(dataset, preferred)
    if symbol is None:
        # A change of origin is a uniform Cartesian translation and leaves
        # all molecular geometry, atom ordering, and bonding intact.
        lattice = get_lattice(output)
        transform = dataset_value(dataset, 'transformation_matrix')
        origin = dataset_value(dataset, 'origin_shift')
        shifted = output.copy()
        shifted.setXYZ(output.getXYZ() + np.linalg.solve(transform, origin) @ lattice)
        shifted_dataset = get_symmetry_dataset(shifted, symprec)
        symbol = matching_space_group(shifted_dataset, preferred)
        if symbol is not None:
            output, dataset = shifted, shifted_dataset
        else:
            # Change basis, not physical geometry. Do not use std_lattice,
            # which idealizes the cell metric at the chosen tolerance.
            if not np.isclose(abs(np.linalg.det(transform)), 1., atol=1e-7):
                raise ValueError('Unsupported cell-size reduction for ' + st.title)
            basis = np.linalg.inv(transform).T @ lattice
            standardized = output.copy()
            standardized.setXYZ(output.getXYZ() + origin @ basis)
            set_lattice(standardized, basis)
            dataset = get_symmetry_dataset(standardized, symprec)
            symbol = matching_space_group(dataset, preferred)
            if symbol is None:
                raise ValueError('No matching space-group setting for ' + st.title)
            output = standardized

    output.property['s_pdb_PDB_CRYST1_Space_Group'] = symbol
    output.property['i_matsci_Space_Group_ID'] = int(dataset_value(dataset, 'number'))
    update_cell_properties(output, st)
    return output


def is_primitive(st, symprec=0.01):
    """Return whether the standard conventional cell contains more atoms."""
    ds = get_symmetry_dataset(st, symprec)
    std = dataset_value(ds, 'std_positions')
    return len(std) > st.atom_total


def process_folder(folder_path, csd, symprec=0.01):
    """Repair all structures and replace the input only after full success."""
    maegz_path = os.path.join(folder_path, 'r2scand3_' + csd + '.maegz')
    if not os.path.exists(maegz_path):
        return 0

    structs = list(structure.StructureReader(maegz_path))
    converted = 0
    modified = False

    for i, st in enumerate(structs):
        title = st.property.get('s_m_title', '')
        try:
            conv_st = correct_space_group(st, symprec=symprec)
            if (dict(conv_st.property) == dict(st.property)
                    and np.array_equal(conv_st.getXYZ(), st.getXYZ())
                    and conv_st.atom_total == st.atom_total):
                continue

            conv_sg = conv_st.property.get(
                's_pdb_PDB_CRYST1_Space_Group', '').strip()

            structs[i] = conv_st
            converted += 1
            modified = True
            print('  Rank ' + str(i + 1) + ': ' + title
                  + ' -> ' + conv_sg
                  + ' (' + str(conv_st.atom_total) + ' atoms)')
        except Exception as e:
            raise RuntimeError('Rank ' + str(i + 1) + ': ' + title) from e

    if modified:
        fd, temporary = tempfile.mkstemp(suffix='.maegz', dir=folder_path)
        os.close(fd)
        try:
            with structure.StructureWriter(temporary) as writer:
                for st in structs:
                    writer.append(st)
            os.replace(temporary, maegz_path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

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
            print(folder + ': repaired ' + str(converted) + ' structures')
        total_converted += converted

    print('\nTotal repaired: ' + str(total_converted) + ' structures')


if __name__ == '__main__':
    main()
