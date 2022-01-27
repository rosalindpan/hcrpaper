import unittest
import crutil

class TestModel(unittest.TestCase):

	def test_remove_gaps1(self):
		seq = 'A-P-A-PAAAP'
		res = crutil.remove_gaps(seq)
		expected = 'APAPAAAP'
		self.assertTrue(res == expected)

	def test_remove_gaps2(self):
		seq = 'A-P-A-PAAAP-'
		res = crutil.remove_gaps(seq)
		expected = 'APAPAAAP'
		self.assertTrue(res == expected)

	def test_remove_gaps3(self):
		seq = '-A-P-A-PAAAP'
		res = crutil.remove_gaps(seq)
		expected = 'APAPAAAP'
		self.assertTrue(res == expected)

	def test_remove_gaps4(self):
		seq = 'A-P-A---PAAAP'
		res = crutil.remove_gaps(seq)
		expected = 'APAPAAAP'
		self.assertTrue(res == expected)

	def test_get_mapping1(self):
		seq = 'A-P-A--P'
		expected = {0: 0, 1:2, 2:4, 3:7}
		res = crutil.get_mapping(seq)
		self.assertEqual(len(expected), len(res))
		for k,v in res.items():
			self.assertEqual(v, expected[k])

	def test_get_mapping2(self):
		seq = '-A-P-A--P'
		expected = {0: 1, 1:3, 2:5, 3:8}
		res = crutil.get_mapping(seq)
		self.assertEqual(len(expected), len(res))
		for k,v in res.items():
			self.assertEqual(v, expected[k])

	def test_get_mapping3(self):
		seq = '-A-P-A--P-'
		expected = {0: 1, 1:3, 2:5, 3:8}
		res = crutil.get_mapping(seq)
		self.assertEqual(len(expected), len(res))
		for k,v in res.items():
			self.assertEqual(v, expected[k])

	def test_get_mapping4(self):
		testseq = "EKE-EKEK--EK"
		res = crutil.get_mapping(testseq)
		expected = {0:0, 1:1, 2:2, 3:4, 4:5, 5:6, 6:7, 7:10, 8:11}
		for k,v in expected.items():
			self.assertEqual(res[k], v)

	def test_get_mapping5(self):
		testseq = "-KE-EKEK--E-"
		res = crutil.get_mapping(testseq)
		expected = {0:1, 1:2, 2:4, 3:5, 4:6, 5:7, 6:10}
		for k,v in expected.items():
			self.assertEqual(res[k], v)

	def test_combine_regions1(self):
		regions = [(1, 11), (11, 20)]
		res = crutil.combine_regions(regions)
		expected = [(1, 20)]
		self.assertListEqual(res, expected)

	def test_combine_regions2(self):
		regions = [(1, 12), (11, 20)]
		res = crutil.combine_regions(regions)
		expected = [(1, 20)]
		self.assertListEqual(res, expected)

	def test_combine_regions3(self):
		regions = [(11, 20), (1, 12)]
		res = crutil.combine_regions(regions)
		expected = [(1, 20)]
		self.assertListEqual(res, expected)

	def test_combine_regions4(self):
		regions = [(1, 12), (11, 20), (30, 40)]
		res = crutil.combine_regions(regions)
		expected = [(1, 20), (30, 40)]
		self.assertListEqual(res, expected)

	def test_combine_regions5(self):
		regions = [(1, 12), (11, 20), (19, 40)]
		res = crutil.combine_regions(regions)
		expected = [(1, 40)]
		self.assertListEqual(res, expected)

	def test_combine_regions6(self):
		regions = [(1, 12), (19, 40), (11, 20)]
		res = crutil.combine_regions(regions)
		expected = [(1, 40)]
		self.assertListEqual(res, expected)

	def test_combine_regions7(self):
		regions = [(1, 10), (12, 20)]
		expected = regions
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions8(self):
		regions = [(1, 10), (9, 20)]
		expected = [(1, 20)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions9(self):
		regions = [(1, 10), (9, 20), (18, 100)]
		expected = [(1, 100)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions10(self):
		regions = [(1, 10), (9, 20), (50, 100)]
		expected = [(1, 20), (50, 100)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions11(self):
		regions = [(1, 10), (12, 20), (18, 100)]
		expected = [(1, 10), (12, 100)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions12(self):
		regions = [(1, 10), (8, 20), (25, 50), (30, 70)]
		expected = [(1, 20), (25, 70)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions13(self):
		regions = [(1, 10), (8, 20), (21, 22), (25, 50), (30, 70)]
		expected = [(1, 20), (21, 22), (25, 70)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

	def test_combine_regions14(self):
		regions = [(1, 10), (10, 20)]
		expected = [(1, 20)]
		combined = crutil.combine_regions(regions)
		self.assertListEqual(expected, combined)

if __name__ == '__main__':
    unittest.main(verbosity=3)