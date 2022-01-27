import unittest
import charge_distribution as cd

class TestModel(unittest.TestCase):

	def test_repeats_2(self):
		r2 = 'EK' * 25
		kappa = cd.get_avg_kappa(r2, 5, 6)
		self.assertAlmostEqual(kappa, 0.0009, delta=1e-4)

	def test_repeats_3(self):
		r3 = 'EEEKKK' * 8 + 'EK'
		kappa = cd.get_avg_kappa(r3, 5, 6)
		self.assertAlmostEqual(kappa, 0.0025, delta=1e-4)

	def test_repeats_4(self):
		# sv7 in Das paper
		r4 = 'EEEEKKKK' * 6 + 'EK'
		kappa = cd.get_avg_kappa(r4, 5, 6)
		self.assertAlmostEqual(kappa, 0.0450, delta=1e-3)

	def test_repeats_5(self):
		# sv19 in Das paper
		r5 = 'EEEEEKKKKK' * 5
		kappa = cd.get_avg_kappa(r5, 5, 6)
		print(kappa)
		self.assertAlmostEqual(kappa, 0.1941, delta=1e-2)

	def test_repeats_25(self):
		# sv30 in Das paper
		r25 = 'E' * 25 + 'K' * 25
		kappa = cd.get_avg_kappa(r25, 5, 6)
		self.assertAlmostEqual(kappa, 1.0000, delta=1e-3)

	def test_random_1(self):
		# sv6 in Das paper
		ran1 = 'EEEKKEKKEEKEEKKEKKEKEEEKKKEKEEKKEEEKKKEKEEEEKKKKEK'
		kappa = cd.get_avg_kappa(ran1, 5, 6)
		self.assertAlmostEqual(kappa, 0.0273, delta=1e-2)

	def test_longer_length(self):
		r2 = 'EK' * 50
		kappa = cd.get_avg_kappa(r2, 5, 6)
		self.assertAlmostEqual(kappa, 0.0009, delta=1e-4)

	def test_longer_length_2(self):
		r5 = 'EEEEEKKKKK' * 50
		kappa = cd.get_avg_kappa(r5, 5, 6)
		self.assertAlmostEqual(kappa, 0.1941, delta=1e-2)

if __name__ == '__main__':
    unittest.main(verbosity=3)