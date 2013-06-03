# -*- coding: utf-8 -*-
"""

Created on Fri May 31 16:33:03 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

A test module for Utility library

"""

import sys
import gzip
import unittest
import Bio.Seq

sys.path.append("..")

import GClib.Utility

class TestFastaFile(unittest.TestCase):
    #To verify a seqObject
    def setUp(self):
        self.seqObj = list(Bio.SeqIO.parse(gzip.open("chr21.fa.gz"), "fasta"))[0]
        
        #Open the sequence with fastafile module
        self.test_seqObj = GClib.Utility.FastaFile("chr21.fa.gz")
    
    def test_GetNextSeq(self):
        #test getting the next_sequence
        self.assertEqual(self.seqObj, self.test_seqObj.GetNextSeq())
    
if __name__ == "__main__":
    unittest.main()   