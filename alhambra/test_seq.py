import unittest
import alhambra.seq as seq

class SeqTests(unittest.TestCase):
    seq1 = "aggagtacc"
    seq2 = "aggn"
    seq3 = "gngc"
    seq4 = "hggc"
    def test_merge_correct(self):
        # same sequence
        self.assertEqual(seq.merge(self.seq1,self.seq1),self.seq1)
        # constraints on both sides
        self.assertEqual(seq.merge(self.seq2,self.seq4),'aggc')
        self.assertEqual(seq.merge(self.seq4,self.seq2),'aggc')
    def test_merge_conflicts(self):
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq1,self.seq2)
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq2,self.seq3)
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq3,self.seq4)
    def test_is_null(self):
        self.assertTrue( seq.is_null("") )
        self.assertTrue( seq.is_null("nnnnnnnnnnnnnnnnNNNNNNN") )
        self.assertFalse( seq.is_null("nnagaagannna") )
        self.assertFalse( seq.is_null("nnnnnnnynnnn") )
    def test_is_definite(self):
        import random
        defbases = frozenset(('a','g','t','c'))
        ambbases = frozenset(seq._L_TO_N.keys()) - defbases
        for _ in range(0,100):
            self.assertTrue(
                seq.is_definite(
                    ''.join(random.choices(tuple(defbases),k=50))) )
            self.assertFalse(
                seq.is_definite(
                    ''.join(random.choices(tuple(ambbases),k=50))) )
            self.assertFalse(
                seq.is_definite(
                    ''.join(
                        random.choices(
                            tuple(frozenset.union(defbases,ambbases)),k=50))))
    def test_space(self):
        self.assertTrue( seq.is_null('Nn nn  nn\tnn\nnn') )
        self.assertTrue( seq.is_null('       \t\n  \n') )
        self.assertFalse( seq.is_definite('    \t\n  \n'))
        self.assertFalse( seq.is_null('nnnnnnnn      \t a \n nnnnnn'))
        self.assertFalse( seq.is_definite('ag  n n\t\nn') )
        self.assertTrue( seq.is_definite('ag  c t\t\ntt') )
            
if __name__ == '__main__':
    unittest.main()
