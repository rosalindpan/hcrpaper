import unittest
import numpy as np
import parseaf as pa

# python3 -W ignore parseaf_test.py

path_to_af_data = '/mnt/d/research/drummond-lab/data/yeast-alphafold-output-unzipped/'

class TestModel(unittest.TestCase):

	def test_read_af_output1(self):
		uniprot_id = 'A2P2R3'
		af_pdb = pa.read_af_output(path_to_af_data, uniprot_id)
		self.assertTrue(af_pdb is not None)

	def test_read_af_output2(self):
		uniprot_id = 'P0CX25'
		af_pdb = pa.read_af_output(path_to_af_data, uniprot_id)
		self.assertTrue(af_pdb is not None)

	def test_read_af_output3(self):
		uniprot_id = 'TRASH'
		af_pdb = pa.read_af_output(path_to_af_data, uniprot_id)
		self.assertTrue(af_pdb is None)

	def test_get_percent_helix1(self):
		res = pa.get_percent_helix(path_to_af_data, 'P31376', 0, 20)
		expected = 0
		self.assertEqual(res, expected)

	def test_get_percent_helix2(self):
		res = pa.get_percent_helix(path_to_af_data, 'P31376', 50, 70)
		expected = 11 / 21
		self.assertEqual(res, expected)

	def test_get_percent_helix3(self):
		res = pa.get_percent_helix(path_to_af_data, 'P39730', 70, 90)
		expected = 14 / 21
		self.assertEqual(res, expected)

	def test_get_freq_values_in_range1(self):
		arr = np.array([1, 2, 3, 4, 5])
		res = pa.get_freq_values_in_range(arr, 3, 6)
		exp = 0.6
		self.assertEqual(res, exp)

	def test_get_freq_values_in_range2(self):
		arr = np.array([1, 2, 3, 4, 5])
		res = pa.get_freq_values_in_range(arr, 7, 8)
		exp = 0
		self.assertEqual(res, exp)

	def test_get_freq_values_in_range3(self):
		arr = np.array([1, 2, 3, 4, 5])
		res = pa.get_freq_values_in_range(arr, 1, 5)
		exp = 1
		self.assertEqual(res, exp)

	def test_get_structure_label1(self):
		res = pa.get_structure_label(path_to_af_data, 'P31376', 0, 20)
		expected = 'disordered'
		self.assertEqual(res, expected)

	def test_get_structure_label2(self):
		res = pa.get_structure_label(path_to_af_data, 'P31376', 50, 70)
		expected = 'unclassified'
		self.assertEqual(res, expected)

	def test_get_structure_label3(self):
		res = pa.get_structure_label(path_to_af_data, 'P39730', 100, 120)
		expected = 'helix'
		self.assertEqual(res, expected)

	def test_get_structure_label4(self):
		res = pa.get_structure_label(path_to_af_data, 'P39730', 100, 120,
									bfactor_cutoff=0.95)
		expected = 'unclassified'
		self.assertEqual(res, expected)

	def test_get_structure_label5(self):
		res = pa.get_structure_label(path_to_af_data, 'P39730', 100, 150,
									helical_cutoff=0.9)
		expected = 'unclassified'
		self.assertEqual(res, expected)





if __name__ == '__main__':
    unittest.main(verbosity=3)
    print("If you're not Rosalind, please change path_to_af_data before running the tests")