"""Run with the Schrödinger Python environment: python3 -m unittest discover -s tests."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from schrodinger import structure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import primitive_to_conventional as repair


def read(folder, code, rank):
    path = ROOT / 'data' / folder / f'r2scand3_{code}.maegz'
    return list(structure.StructureReader(str(path)))[rank - 1]


class SpaceGroupTests(unittest.TestCase):
    def check_group(self, original, number):
        output = repair.correct_space_group(original)
        self.assertEqual(output.atom_total, original.atom_total)
        self.assertEqual(output.property['i_matsci_Space_Group_ID'], number)
        dataset = repair.get_symmetry_dataset(output)
        self.assertEqual(repair.dataset_value(dataset, 'number'), number)
        symbol = output.property['s_pdb_PDB_CRYST1_Space_Group']
        self.assertEqual(repair.matching_space_group(dataset, symbol), symbol)
        # Coordinate changes must be a uniform translation, not a distortion.
        delta = output.getXYZ() - original.getXYZ()
        np.testing.assert_allclose(delta, np.tile(delta[0], (len(delta), 1)), atol=1e-9)
        self.assertAlmostEqual(abs(np.linalg.det(repair.get_lattice(original))),
                               abs(np.linalg.det(repair.get_lattice(output))), places=7)
        for key in ('r_matsci_Etot_(Ry)',
                    'r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)'):
            self.assertEqual(original.property.get(key), output.property.get(key))
        again = repair.correct_space_group(output)
        self.assertEqual(dict(again.property), dict(output.property))
        np.testing.assert_array_equal(again.getXYZ(), output.getXYZ())
        return output

    def test_wrong_p1_without_cell_expansion(self):
        st = read('Cocaine_COCAIN', 'COCAIN', 1)
        st.property['s_pdb_PDB_CRYST1_Space_Group'] = 'P 1'
        st.property['i_matsci_Space_Group_ID'] = 1
        self.check_group(st, 4)

    def test_origin_and_basis_settings(self):
        for folder, code, rank, number in (
                ('Cocaine_COCAIN', 'COCAIN', 27, 144),
                ('Glycine_GLYCIN', 'GLYCIN', 1, 144),
                ('PF-998425_VIPRES', 'VIPRES', 26, 78),
                ('PF-998425_VIPRES', 'VIPRES', 40, 78),
                ('Rotigotine_RALMOG', 'RALMOG', 3, 145)):
            with self.subTest(folder=folder, rank=rank):
                st = read(folder, code, rank)
                # Force a nonstandard basis and origin even after the checked-in
                # regression structures themselves have been repaired.
                shear = np.array([[1., 0, 0], [1., 1., 0], [0, 0, 1.]])
                repair.set_lattice(st, shear @ repair.get_lattice(st))
                st.setXYZ(st.getXYZ() + np.array([.137, .219, .311]))
                st.property['s_pdb_PDB_CRYST1_Space_Group'] = 'P 1'
                self.check_group(st, number)

    def test_keep_valid_nonstandard_aspirin_setting(self):
        st = read('Asprin_ACSALA', 'ACSALA', 1)
        output = self.check_group(st, 14)
        self.assertEqual(st.property['s_pdb_PDB_CRYST1_Space_Group'],
                         output.property['s_pdb_PDB_CRYST1_Space_Group'])
        np.testing.assert_array_equal(repair.get_lattice(st), repair.get_lattice(output))
        np.testing.assert_array_equal(st.getXYZ(), output.getXYZ())

    def test_failed_detection_is_not_p1(self):
        st = read('Cocaine_COCAIN', 'COCAIN', 1)
        with patch.object(repair.spglib, 'get_symmetry_dataset', return_value=None):
            with self.assertRaises(ValueError):
                repair.correct_space_group(st)

    def test_expanded_volume_and_z_do_not_rescale_dft_energy(self):
        original = SimpleNamespace(atom_total=10, property={
            's_m_Unit_Cell_Formula': '(C5H5)2', 'r_matsci_Etot_(Ry)': -123.45,
            'r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)': 1.23})
        lattice = np.array([[2., 0, 0], [1., 3, 0], [.2, .1, 4]])
        props = {f'r_chorus_box_{a}{c}': float(lattice[i, j])
                 for i, a in enumerate('abc') for j, c in enumerate('xyz')}
        props['r_m_Unit_Cell_Volume/Ang.^3'] = 12.
        converted = SimpleNamespace(atom_total=20, property=props)
        repair.update_cell_properties(converted, original)
        self.assertAlmostEqual(converted.property['r_m_Unit_Cell_Volume/Ang.^3'], 24.)
        self.assertEqual(converted.property['s_m_Unit_Cell_Formula'], '(C5H5)4')
        self.assertEqual(converted.property['r_matsci_Etot_(Ry)'], -123.45)
        self.assertEqual(converted.property['r_matsci_Relative_DFT_energy_per_mol_(kcal/mol)'], 1.23)


if __name__ == '__main__':
    unittest.main()
