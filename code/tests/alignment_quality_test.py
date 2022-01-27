import unittest
from Bio import AlignIO
from scipy import stats
import crutil
import alignment_quality as aq

all_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

class TestModel(unittest.TestCase):

    def test_get_region_range1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        res = aq.get_region_range(msa, 418, 501)
        expected = (533, 791)
        self.assertEqual(res[0], expected[0])
        self.assertEqual(res[1], expected[1])

    def test_get_region_range2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        for record in msa:
            if record.id == "Saccharomyces":
                refseq = crutil.remove_gaps(record.seq)
        res = aq.get_region_range(msa, 418, 501, refseq=refseq)
        expected = (533, 791)
        self.assertEqual(res[0], expected[0])
        self.assertEqual(res[1], expected[1])

    def test_get_region_range3(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        refseq = "YQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ"
        res = aq.get_region_range(msa, 418, 501, refseq=refseq)
        expected = (533, 791)
        self.assertEqual(res[0], expected[0])
        self.assertEqual(res[1], expected[1])

    def test_get_region_range4(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        res = aq.get_region_range(msa, 418, 3000)
        expected = (533, 875)
        self.assertEqual(res[0], expected[0])
        self.assertEqual(res[1], expected[1])

    def test_get_region_range5(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        res = aq.get_region_range(msa, -10, 501)
        expected = (13, 791)
        self.assertEqual(res[0], expected[0])
        self.assertEqual(res[1], expected[1])

    def test_extract_region_msa(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 418, 501)
        for record in region_msa:
            if record.id == "Saccharomyces":
                result = crutil.remove_gaps(record.seq)
        expected = "YQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ"
        self.assertEqual(result, expected)

    def test_extract_random_region_from_msa1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        random_msa = aq.extract_random_region_from_msa(msa, 50)
        self.assertTrue(len(random_msa[0].seq), 50)

    def test_extract_random_region_from_msa2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[:, :50]
        random_msa = aq.extract_random_region_from_msa(msa, 100)
        self.assertTrue(len(random_msa[0].seq), 50)

    def test_get_seq_len(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        result = aq.get_seq_len(msa)
        self.assertEqual(len(result), 34)
        self.assertEqual(result['Saccharomyces cerevisiae sce sce|PABP_YEAST'], 577)

    def test_filter_msa1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        filtered_msa = aq.filter_msa(msa)
        removed = ['Neurospora crassa ncr ncr|PABP_NEUCR',
                   'Pichia kudriavzevii pku pku|A0A099NS00_PICKU']
        descs = []
        for record in filtered_msa:
            descs.append(record.description)
        self.assertEqual(len(descs), 32)
        for r in removed:
            self.assertTrue(r not in descs)

    def test_filter_msa2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 418, 501)
        filtered_msa = aq.filter_msa(region_msa)
        removed = ['Lipomyces starkeyi lst lst|Lipst1_1_2201',
                   'Rhodotorula graminis rgm rgm|Rhoba1_1_66719',
                   'Pichia kudriavzevii pku pku|A0A099NS00_PICKU']
        descs = []
        for record in filtered_msa:
            descs.append(record.description)
        self.assertEqual(len(descs), 31)
        for r in removed:
            self.assertTrue(r not in descs)

    def test_get_gap_only_indices1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 418, 438)
        filtered_msa = aq.filter_msa(region_msa)
        result = aq.get_gap_only_indices(filtered_msa)
        expected = [16, 17, 18, 19, 20]
        self.assertListEqual(result, expected)

    def test_get_gap_only_indices2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 465, 480)
        filtered_msa = aq.filter_msa(region_msa)
        result = aq.get_gap_only_indices(filtered_msa)
        self.assertTrue(len(result) != 0)
        expected = [20, 39, 40]
        self.assertListEqual(result, expected)

    def test_remove_gap_only_indices1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 418, 438)
        filtered_msa = aq.filter_msa(region_msa)
        initial_len = len(filtered_msa[0].seq)
        aq.remove_gap_only_indices(filtered_msa)
        final_len = len(filtered_msa[0].seq)
        self.assertEqual(initial_len, final_len + 5)

    def test_remove_gap_only_indices2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 465, 480)
        filtered_msa = aq.filter_msa(region_msa)
        initial_len = len(filtered_msa[0].seq)
        aq.remove_gap_only_indices(filtered_msa)
        final_len = len(filtered_msa[0].seq)
        self.assertEqual(initial_len, final_len + 3)

    def test_get_gap_frequency_per_seq1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")
        region_msa = aq.extract_region_msa(msa, 420, 425)[:5, :]
        result = aq.get_gap_frequency_per_seq(region_msa)
        self.assertEqual(result, 2/7)

    def test_get_gap_frequency_per_seq2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[29:, :]
        region_msa = aq.extract_region_msa(msa, 420, 425)
        result = aq.get_gap_frequency_per_seq(region_msa)
        self.assertEqual(result, 1/7)

    def test_compute_alignment_quality(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, :]
        region_msa = aq.extract_region_msa(msa, 420, 425)
        result = aq.compute_alignment_quality(region_msa)
        self.assertEqual(result, 1/12)

    def test_count_aas_in_column1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:50]
        expected = {'D': 6, 'E': 1}
        result = aq.count_aas_in_column(msa, 0, ['D', 'E'])
        self.assertEqual(len(expected), len(result))
        for k,v in expected.items():
            self.assertEqual(v, result[k])

    def test_count_aas_in_column2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:50]
        expected = {'A': 1, 'I': 4, 'Q': 1, 'V': 1, 'P': 0}
        result = aq.count_aas_in_column(msa, 1, ['A', 'I', 'Q', 'V', 'P'])
        self.assertEqual(len(expected), len(result))
        for k,v in expected.items():
            self.assertEqual(v, result[k])

    def test_get_column_aa_distribution1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:50]
        expected = [6/7, 1/7]
        result = aq.get_column_aa_distribution(msa, 0, ['D', 'E'])
        self.assertEqual(len(expected), len(result))
        for i in range(2):
            self.assertAlmostEqual(expected[i], result[i], delta=1e-6)

    def test_get_column_aa_distribution2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:50]
        expected = [1/7, 4/7, 1/7, 1/7, 0]
        result = aq.get_column_aa_distribution(msa, 1, ['A', 'I', 'Q', 'V', 'P'])
        self.assertEqual(len(expected), len(result))
        for i in range(5):
            self.assertAlmostEqual(expected[i], result[i], delta=1e-6)

    def test_get_seq_divergence1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:20]
        prob1 = [6/7, 1/7]
        prob2 = [1/7, 4/7, 1/7, 1/7]
        prob3 = [1/8, 1/8, 1/8, 1/2, 1/8]
        prob4 = [1/2, 3/8, 1/8]
        prob5 = [1/8, 7/8]
        probs = [prob1, prob2, prob3, prob4, prob5]
        expected = 0
        for i in range(5):
            expected += stats.entropy(probs[i]) / 5
        result = aq.get_seq_divergence(msa, all_aa)
        self.assertAlmostEqual(expected, result)

    def test_get_seq_divergence1(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 15:20]
        prob1 = [6/7, 1/7]
        prob2 = [1/7, 4/7, 1/7, 1/7]
        prob3 = [1/8, 1/8, 1/8, 1/2, 1/8]
        prob4 = [1/2, 3/8, 1/8]
        prob5 = [1/8, 7/8]
        probs = [prob1, prob2, prob3, prob4, prob5]
        expected = 0
        for i in range(5):
            expected += stats.entropy(probs[i]) / 5
        result = aq.get_seq_divergence(msa, all_aa)
        self.assertAlmostEqual(expected, result)

    def test_get_seq_divergence2(self):
        msa = AlignIO.read("../data/sequence/PAB1-aybrah.fa", "fasta")[26:, 20:25]
        prob1 = [7/8, 1/8]
        prob2 = [1/8, 3/8, 1/2]
        prob3 = [1/8, 1/8, 3/4]
        prob4 = [1/4, 3/4]
        probs = [prob1, prob2, prob3, prob4]
        expected = 0
        for i in range(4):
            expected += stats.entropy(probs[i]) / 4
        result = aq.get_seq_divergence(msa, all_aa)
        self.assertAlmostEqual(expected, result)
                
if __name__ == '__main__':
    unittest.main(verbosity=3)