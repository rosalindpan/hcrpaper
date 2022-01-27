from doctest import testsource
import re
import unittest

import composition as comp
import numpy as np
import string


class TestModel(unittest.TestCase):

    def test_get_character_freqs(self):
        testseq = list(comp.aas[0])
        res = list(comp.get_character_freqs(comp.aas)[0])
        expected = list(np.ones(20)*0.05)
        self.assertEqual(res, expected)

    def test_get_aa_freqs(self):
        # should work on lowercase things with incorrect characters
        testseq = string.ascii_lowercase + "*"
        res = comp.get_aa_freqs(testseq)
        expected = np.ones(20)*0.05
        for i, freq in enumerate(res):
            self.assertEqual(freq, expected[i])
    
    def test_get_aa_freqs_sparse(self):
        testseq = "A"
        res = comp.get_aa_freqs(testseq)
        expected = [1] + list(np.zeros(19))
        for i, freq in enumerate(res):
            self.assertEqual(freq, expected[i])

    def test_get_aa_counts(self):
        # should work on lowercase things with incorrect characters
        testseq = string.ascii_lowercase + "*"
        res = comp.get_aa_counts(testseq)
        expected = np.ones(20)
        for i, freq in enumerate(res):
            self.assertEqual(freq, expected[i])

    def test_get_aa_counts_otherchars(self):
        # should work on lowercase things with incorrect characters
        testseq = string.ascii_lowercase + "ZZZZZZZZZZZZZZ"
        res = comp.get_aa_counts(testseq)
        expected = np.ones(20)
        print(expected)
        for i, freq in enumerate(res):
            self.assertEqual(freq, expected[i])

    def test_seq_complexity_k1(self):
        # complexity should be the same if only the aa part of the sequence is considered
        res1 = comp.calculate_seq_complexity("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        res2 = comp.calculate_seq_complexity("ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZZZZZZZZZZZZZZZ")
        self.assertEqual(res1, res2)

    def test_seq_complexity_zero(self):
        testseq = "AAAAAAAAZZZZZZZ"
        res = comp.calculate_seq_complexity(testseq)
        self.assertEqual(res, 0)

    def test_seq_complexity_k2(self):
        testseq = "AAAAAAA"
        res = comp.calculate_seq_complexity(testseq, typ="k2", norm="uniform")
        self.assertEqual(res, 0)

if __name__ == '__main__':
    unittest.main(verbosity=3)