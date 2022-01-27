from math import exp
import unittest
from unittest.case import expectedFailure
import fractional_charge as fc

class TestModel(unittest.TestCase):

	def test_get_fractional_charge(self):
		testseq = "AAAAEKDR"
		res = fc.get_fractional_charge(testseq)
		expected = 0.5
		self.assertAlmostEqual(res, expected)

	def test_make_weights_odd(self):
		w = fc.make_weights(5)
		expected = [0, 0.5, 1, 0.5, 0]
		for i in range(5):
			self.assertEqual(w[i], expected[i])

	def test_make_weights_even(self):
		w = fc.make_weights(6)
		expected = [0, 0.5, 1, 1, 0.5, 0]
		for i in range(6):
			self.assertEqual(w[i], expected[i])

	def test_make_weights_larger(self):
		w = fc.make_weights(7)
		expected = [0, 1/3, 2/3, 1, 2/3, 1/3, 0]
		for i in range(7):
			self.assertEqual(w[i], expected[i])

	def test_get_fractional_charge_weighted_fully_charged(self):
		testseq = "EKEKEK"
		res = fc.get_fractional_charge_weighted(testseq)
		self.assertEqual(res, 1)

	def test_get_fractional_charge_weighted_uncharged(self):
		testseq = "PAPAPAPA"
		res = fc.get_fractional_charge_weighted(testseq)
		self.assertEqual(res, 0)

	def test_get_fractional_charge_weighted_partially_charged_1(self):
		testseq = "AEKEKEK"
		res = fc.get_fractional_charge_weighted(testseq)
		self.assertEqual(res, 1)

	def test_get_fractional_charge_weighted_partially_charged_2(self):
		testseq = "EKEAKEK"
		res = fc.get_fractional_charge_weighted(testseq)
		self.assertEqual(res, 2/3)

	def test_rolling_fractional_charges_weighted_fully_charged(self):
		# expected result = [1, 1, 1, 1]
		testseq = "EKEKEK"
		res = fc.rolling_fractional_charges_weighted(testseq, 3)
		for i in res:
			self.assertEqual(i, 1)
		self.assertEqual(len(res), 4)

	def test_rolling_fractional_charges_weighted_uncharged(self):
		# expected result = [0, 0, 0, 0, 0]
		testseq = "AAAAAAAA"
		res = fc.rolling_fractional_charges_weighted(testseq, 4)
		for i in res:
			self.assertEqual(i, 0)
		self.assertEqual(len(res), 5)

	def test_rolling_fractional_charges_partially_charged(self):
		testseq = "AAEKEAEKEAA"
		res = fc.rolling_fractional_charges_weighted(testseq, 5)
		self.assertEqual(len(res), 7)
		expected = [0.75, 1, 0.75, 0.5, 0.75, 1, 0.75]
		for i in range(7):
			self.assertEqual(res[i], expected[i])

	def test_is_above_threshold_all(self):
		f_charges = [0.5, 0.6, 0.7, 0.8, 0.9, 1]
		res = fc.is_above_threshold(f_charges, 0.2)
		expected = [(0, 0.5), (1, 0.6), (2, 0.7),
					(3, 0.8), (4, 0.9), (5, 1)]
		for i in range(6):
			self.assertEqual(res[i], expected[i])

	def test_is_above_threshold_all(self):
		f_charges = [0.5, 0.6, 0.7, 0.8, 0.9, 1]
		res = fc.is_above_threshold(f_charges, 0.2)
		expected = [(0, 0.5), (1, 0.6), (2, 0.7),
					(3, 0.8), (4, 0.9), (5, 1)]
		for i in range(6):
			self.assertEqual(res[i], expected[i])

	def test_is_above_threshold_none(self):
		f_charges = [0.1, 0.05, 0.1, 0, 0, 0]
		res = fc.is_above_threshold(f_charges, 0.2)
		self.assertTrue(len(res) == 0)

	def test_is_above_threshold_some_1(self):
		f_charges = [0.2, 0.05, 0.3, 0.4, 0.2, 0.0, 0.1]
		res = fc.is_above_threshold(f_charges, 0.2)
		expected = [(0, 0.2), (2, 0.3), (3, 0.4), (4, 0.2)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_is_above_threshold_some_2(self):
		f_charges = [0.0, 0.05, 0.6, 0.4, 0.2, 0.0, 0.3]
		res = fc.is_above_threshold(f_charges, 0.3)
		expected = [(2, 0.6), (3, 0.4), (6, 0.3)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_identify_charged_regions_all(self):
		f_charges = [0.2, 0.2, 0.3, 0.2, 0.3, 0.4, 0.2]
		res = fc.identify_charged_regions(f_charges, 0.2)
		expected = [[(0, 0.2), (1, 0.2), (2, 0.3),
			 		(3, 0.2), (4, 0.3), (5, 0.4), (6, 0.2)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_none(self):
		f_charges = [0.1, 0.05, 0.1, 0, 0, 0]
		res = fc.identify_charged_regions(f_charges, 0.2)
		self.assertEqual(len(res), 0)

	def test_identify_charged_regions_some_1(self):
		f_charges = [0.2, 0.2, 0.3, 0.2, 0.3, 0.4, 0.2, 0,
					 0.2, 0.2, 0.3, 0.2, 0.3, 0.4, 0.2]
		res = fc.identify_charged_regions(f_charges, 0.2)
		expected = [[(0, 0.2), (1, 0.2), (2, 0.3),
			 		(3, 0.2), (4, 0.3), (5, 0.4), (6, 0.2)],
			 		[(8, 0.2), (9, 0.2), (10, 0.3),
			 		(11, 0.2), (12, 0.3), (13, 0.4), (14, 0.2)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_some_2(self):
		f_charges = [0, 0.2, 0.3, 0.2, 0.3, 0.4, 0, 0.2, 0.4]
		res = fc.identify_charged_regions(f_charges, 0.2)
		expected = [[(1, 0.2), (2, 0.3), (3, 0.2), (4, 0.3), (5, 0.4)],
					[(7, 0.2), (8, 0.4)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_some_3(self):
		f_charges = [0, 0.5, 0.7, 0, 0.1, 0, 0.2, 0.3,
		             0.2, 0.3, 0.4, 0, 0.2, 0.4, 0.35]
		res = fc.identify_charged_regions(f_charges, 0.2)
		expected = [[(1, 0.5), (2, 0.7)],
					[(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
					[(12, 0.2), (13, 0.4), (14, 0.35)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_with_tolerance_1(self):
		f_charges = [0, 0.5, 0.7, 0, 0.1, 0, 0.2, 0.3,
		             0.2, 0.3, 0.4, 0, 0, 0.2, 0.4, 0.35]
		res = fc.identify_charged_regions(f_charges, 0.2, tolerance=2)
		expected = [[(1, 0.5), (2, 0.7)],
					[(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4),
					 (13, 0.2), (14, 0.4), (15, 0.35)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_with_tolerance_2(self):
		f_charges = [0, 0.5, 0.7, 0, 0.1, 0, 0.2, 0.3,
		             0.2, 0.3, 0.4, 0, 0, 0.2, 0.4, 0.35, 0.0]
		res = fc.identify_charged_regions(f_charges, 0.2, tolerance=5)
		expected = [[(1, 0.5), (2, 0.7), (6, 0.2), (7, 0.3), (8, 0.2),
		             (9, 0.3), (10, 0.4), (13, 0.2), (14, 0.4), (15, 0.35)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_identify_charged_regions_with_tolerance_3(self):
		f_charges = [0, 0.5, 0.7, 0, 0.2, 0.3,
		             0.2, 0.3, 0.4, 0, 0.2, 0.4, 0.35, 0.0]
		res = fc.identify_charged_regions(f_charges, 0.2, tolerance=3)
		expected = [[(1, 0.5), (2, 0.7), (4, 0.2), (5, 0.3), (6, 0.2),
		             (7, 0.3), (8, 0.4), (10, 0.2), (11, 0.4), (12, 0.35)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_select_min_length_1(self):
		regions = [[(1, 0.5), (2, 0.7)],
				   [(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
				   [(12, 0.2), (13, 0.4), (14, 0.35)]]
		res = fc.select_min_length(regions, 4)
		expected = [[(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])
	
	def test_select_min_length_2(self):
		regions = [[(1, 0.5), (2, 0.7)],
				   [(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
				   [(12, 0.2), (13, 0.4), (14, 0.35)]]
		res = fc.select_min_length(regions, 5)
		expected = [[(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_select_min_length_3(self):
		regions = [[(1, 0.5), (2, 0.7)],
				   [(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
				   [(12, 0.2), (13, 0.4), (14, 0.35)]]
		res = fc.select_min_length(regions, 6)
		self.assertEqual(len(res), 0)

	def test_select_min_length_4(self):
		regions = [[(1, 0.5), (2, 0.7)],
				   [(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
				   [(12, 0.2), (13, 0.4), (14, 0.35)]]
		res = fc.select_min_length(regions, 2)
		expected = [[(1, 0.5), (2, 0.7)],
				   [(6, 0.2), (7, 0.3), (8, 0.2), (9, 0.3), (10, 0.4)],
				   [(12, 0.2), (13, 0.4), (14, 0.35)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_charged_regions_1(self):
		testseq = "EKEKEKEKEKEK"
		res = fc.find_charged_regions(testseq, 5, 0.2, 5, 0)
		expected = [[(0, 1), (1, 1), (2, 1), (3, 1),
		            (4, 1), (5, 1), (6, 1), (7, 1)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_charged_regions_2(self):
		testseq = "AAAAAAA"
		res = fc.find_charged_regions(testseq, 5, 0.2, 5, 0)
		self.assertEqual(len(res), 0)

	def test_find_charged_regions_3(self):
		testseq = "AAEKEAEKEAA"
		# f_charges = [0.75, 1, 0.75, 0.5, 0.75, 1, 0.75]
		res = fc.find_charged_regions(testseq, 5, 0.7, 2, 0)
		expected = [[(0, 0.75), (1, 1), (2, 0.75)],
		            [(4, 0.75), (5, 1), (6, 0.75)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_range_to_seq_1(self):
		testseq = "MVL--TIYPD--EL-VQIVSDKIASNKGKI"
		region = (0, 6)
		res = fc.range_to_seq(testseq, region)
		expected = "MVLTI"
		self.assertEqual(res, expected)

	def test_range_to_seq_2(self):
		testseq = "MVL--TIYPD--EL-VQIVSDKIASNKGKI"
		region = (25, 29)
		res = fc.range_to_seq(testseq, region)
		expected = "NKGKI"
		self.assertEqual(res, expected)

	def test_filter_regions(self):
		testseq = "AAAEKEAEKEA"
		charged_regions = [(0, 5), (6, 10)]
		charge_threshold = 0.6
		res = fc.filter_regions(testseq, charged_regions, charge_threshold)
		expected = [(6,10)]
		self.assertEqual(res, expected)

	def test_find_start_end(self):
		testseq = "--AAE---KEAE--KEA-A"
		res = fc.find_start_end(testseq, 5, 0.5, 2, 0)
		expected = [(2, 18)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_start_end_long(self):
		testseq = "EK" * 100
		res = fc.find_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_start_end_long_2(self):
		testseq = "A" * 100 + "EK" * 100 + "A" * 100
		res = fc.find_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(98, 301)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_start_end_unbalanced(self):
		testseq = "E" * 100 + "K" * 100
		res = fc.find_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_start_end_separated(self):
		testseq = "EK" * 10 + "A" * 20 + "EK" * 10
		res = fc.find_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 21), (38, 59)]
		# not [(0, 19), (40, 59)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_charged_seqs(self):
		testseq = "--AAE---KEAE--KEA-A"
		res = fc.find_charged_seqs(testseq, 5, 0.5, 2, 0)
		expected = ["AAEKEAEKEAA"]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])


# CT
	def test_get_alphabet_fraction(self):
		testseq="AAAAZZZZ"
		res = fc.get_alphabet_fraction(testseq)
		expected = 0.5
		self.assertEqual(res,expected)

	def test_get_alphabet_fraction_lower(self):
		testseq="aaaazzzz"
		res = fc.get_alphabet_fraction(testseq)
		expected = 0.5
		self.assertEqual(res,expected)

	def test_get_alphabet_fraction_weighted_one(self):
		testseq = "AAABBBCDEFGHIJD"
		res = fc.get_alphabet_fraction_weighted(testseq)
		expected = 1
		self.assertEqual(res, expected)

	def test_get_alphabet_fraction_weighted_half(self):
		testseq = "AAAAAZZZZZ"
		res = fc.get_alphabet_fraction_weighted(testseq)
		expected = 0.5
		self.assertEqual(res, expected)


	def test_rolling_alphabet_fraction_weighted_all_first(self):
		testseq = "DDEEEFFGGHHH"
		res = fc.rolling_alphabet_fraction_weighted(testseq, 3)
		#print(len(res))
		for i in res:
			self.assertEqual(i, 1)
		self.assertEqual(len(res), 10)

	def test_rolling_alphabet_fraction_weighted_all_last(self):
		# expected result = [0, 0, 0, 0, 0]
		testseq = "ZZZZZZZZ"
		res = fc.rolling_alphabet_fraction_weighted(testseq, 4)
		for i in res:
			self.assertEqual(i, 0)
		self.assertEqual(len(res), 5)

	def test_rolling_alphabet_fraction_weighted_partial(self):
		# Expected result [1.0, 0.75, 0.25, 0.0]
		testseq = "AAAAZZZZ"
		res = fc.rolling_alphabet_fraction_weighted(testseq, 5)
		self.assertEqual(len(res), 4)
		expected = [1.0, 0.75, 0.25, 0.0]
		for i in range(len(expected)):
			self.assertEqual(res[i], expected[i])
	
	def test_get_alphabet_fraction_half(self):
		testseq = "AZAAAZZZAZ"
		res = fc.get_alphabet_fraction(testseq)
		expected = 0.5
		self.assertAlmostEqual(res, expected)

	def test_filter_alphabet_regions(self):
		testseq = "ZZZAAAZAAAZ"
		regions = [(0, 5), (6, 10)]
		threshold = 0.6
		res = fc.filter_alphabet_regions(testseq, regions, threshold)
		expected = [(6,10)]
		self.assertEqual(res, expected)

	def test_find_alphabet_regions(self):
		testseq = "AABBAABBAABB"
		res = fc.find_alphabet_regions(testseq, 5, 0.2, 5, 0)
		expected = [[(0, 1), (1, 1), (2, 1), (3, 1),
		            (4, 1), (5, 1), (6, 1), (7, 1)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_alphabet_regions_2(self):
		testseq = "ZZZZZZZ"
		res = fc.find_alphabet_regions(testseq, 5, 0.2, 5, 0)
		self.assertEqual(len(res), 0)

	def test_find_alphabet_regions_3(self):
		testseq = "ZZAAAZAAAZZ"
		# f_charges = [0.75, 1, 0.75, 0.5, 0.75, 1, 0.75]
		res = fc.find_alphabet_regions(testseq, 5, 0.7, 2, 0)
		expected = [[(0, 0.75), (1, 1), (2, 0.75)],
		            [(4, 0.75), (5, 1), (6, 0.75)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_alphabet_start_end(self):
		#testseq = "--AAE---KEAE--KEA-A"
		testseq = "--ZZA---AAZA--AAZ-Z"
		res = fc.find_alphabet_start_end(testseq, 5, 0.5, 2, 0)
		expected = [(2, 18)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_alphabet_start_end_long(self):
		testseq = "AA" * 100
		res = fc.find_alphabet_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_alphabet_start_end_long_2(self):
		testseq = "Z" * 100 + "AA" * 100 + "Z" * 100
		res = fc.find_alphabet_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(98, 301)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_alphabet_start_end_unbalanced(self):
		testseq = "A" * 100 + "B" * 100
		res = fc.find_alphabet_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_alphabet_start_end_separated(self):
		#testseq = "EK" * 10 + "A" * 20 + "EK" * 10
		testseq = "AB" * 10 + "X" * 20 + "CD" * 10
		res = fc.find_alphabet_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 21), (38, 59)]
		# not [(0, 19), (40, 59)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_alphabet_start_end_separated_2(self):
		testseq = "LM" * 10 + "N" * 20 + "IJ" * 10
		res = fc.find_alphabet_start_end(testseq, 5, 0.7, 2, 0)
		expected = [(0, 21), (38, 59)]
		# not [(0, 19), (40, 59)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_enriched_fraction(self):
		testseq="AAAAAAAA"
		targets=["L", "V", "I"]
		res = fc.get_enriched_fraction(testseq, targets)
		expected = 0
		self.assertEqual(res,expected)
	
	def test_enriched_fraction_lower(self):
		testseq="aaaallll"
		targets=["L", "V", "I"]
		res = fc.get_enriched_fraction(testseq, targets)
		expected = 0.5
		self.assertEqual(res,expected)
	
	def test_enriched_fraction_all(self):
		testseq="vvvvviiiillllll"
		targets=["L", "V", "I"]
		res = fc.get_enriched_fraction(testseq, targets)
		expected = 1
		self.assertEqual(res,expected)

	def test_get_enriched_fraction_weighted_half(self):
		testseq = "AAAAAZZZZZ"
		targets = ["A", "B", "C"]
		res = fc.get_enriched_fraction_weighted(testseq, targets)
		expected = 0.5
		self.assertEqual(res, expected)
	
	def test_rolling_enriched_fraction_weighted_all_in(self):
		testseq = "AABBBCCCCCCC"
		targets = ["A", "B", "C"]
		res = fc.rolling_enriched_fraction_weighted(testseq, targets, 3)
		#print(len(res))
		for i in res:
			self.assertEqual(i, 1)
		self.assertEqual(len(res), 10)
	
	def test_rolling_enriched_fraction_weighted_partial(self):
		# Expected result [1.0, 0.75, 0.25, 0.0]
		testseq = "AAAAZZZZ"
		targets = ["A", "B", "C"]
		res = fc.rolling_enriched_fraction_weighted(testseq, targets, 5)
		self.assertEqual(len(res), 4)
		expected = [1.0, 0.75, 0.25, 0.0]
		for i in range(len(expected)):
			self.assertEqual(res[i], expected[i])



	def test_filter_enriched_regions(self):
		testseq = "ZZZAAAZAAAZ"
		targets = ["A", "B", "C"]
		regions = [(0, 5), (6, 10)]
		threshold = 0.6
		res = fc.filter_enriched_regions(testseq, targets, regions, threshold)
		expected = [(6,10)]
		self.assertEqual(res, expected)

	def test_find_enriched_regions_allin(self):
		testseq = "AABBAABBAABB"
		targets = ["A", "B", "C"]
		res = fc.find_enriched_regions(testseq, targets, 5, 0.2, 5, 0)
		expected = [[(0, 1), (1, 1), (2, 1), (3, 1),
		            (4, 1), (5, 1), (6, 1), (7, 1)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_enriched_regions_2(self):
		testseq = "ZZZZZZZ"
		targets = ["A", "B", "C"]
		res = fc.find_enriched_regions(testseq, targets, 5, 0.2, 5, 0)
		self.assertEqual(len(res), 0)

	def test_find_enriched_regions_3(self):
		testseq = "ZZAAAZAAAZZ"
		targets = ["A", "B", "C"]
		# f_charges = [0.75, 1, 0.75, 0.5, 0.75, 1, 0.75]
		res = fc.find_enriched_regions(testseq, targets, 5, 0.7, 2, 0)
		expected = [[(0, 0.75), (1, 1), (2, 0.75)],
		            [(4, 0.75), (5, 1), (6, 0.75)]]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertListEqual(res[i], expected[i])

	def test_find_alphabet_start_end(self):
		testseq = "--ZZA---AAZA--AAZ-Z"
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.5, 2, 0)
		expected = [(2, 18)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_enriched_start_end_long(self):
		testseq = "AA" * 100
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_enriched_start_end_long_2(self):
		testseq = "Z" * 100 + "AA" * 100 + "Z" * 100
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.7, 2, 0)
		expected = [(98, 301)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_enriched_start_end_unbalanced(self):
		testseq = "A" * 100 + "B" * 100
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.7, 2, 0)
		expected = [(0, 199)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_enriched_start_end_separated(self):
		#testseq = "EK" * 10 + "A" * 20 + "EK" * 10
		testseq = "AB" * 10 + "X" * 20 + "CC" * 10
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.7, 2, 0)
		expected = [(0, 21), (38, 59)]
		self.assertEqual(len(res), len(expected))
		for i in range(len(res)):
			self.assertEqual(res[i], expected[i])

	def test_find_enriched_start_end_none(self):
		testseq = "LM" * 10 + "N" * 20 + "IJ" * 10
		targets = ["A", "B", "C"]
		res = fc.find_enriched_start_end(testseq, targets, 5, 0.7, 2, 0)
		expected = []
		self.assertEqual(len(res), len(expected))

if __name__ == '__main__':
    unittest.main(verbosity=3)