# -*- coding: utf-8 -*-
"""

Created on Fri May 31 16:33:03 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

A test module for Utility library

"""

import sys
import Bio
import types
import StringIO
import unittest


sys.path.append("..")

import GClib
import GClib.Utility

#Setting libary verbosity
GClib.logger.threshold = 0

class TestLogger(unittest.TestCase):
    #To verify the Logger class
    def setUp(self):
        """Define a Logger class"""
        self.test_outfile = StringIO.StringIO()
        self.test_errfile = StringIO.StringIO()
        self.test_logger = GClib.Utility.Logger(threshold=1,outfile=self.test_outfile,errfile=self.test_errfile)
        
    def test_Log_lower(self):
        """Testing log lower and equal than threshold"""
        
        self.test_logger.log(0, "test")
        self.test_logger.log(1, "test")
        
        self.test_outfile.seek(0)
        data = [line.split(": ")[1] for line in self.test_outfile.readlines()]
        
        self.assertEqual(data, ['test\n', 'test\n'])
        
    def test_Log_higher(self):
        """Testing log higher than threshold (no message)"""
        
        self.test_logger.log(2, "test")
        self.test_outfile.seek(0)
        data = [line.split(": ")[1] for line in self.test_outfile.readlines()]
        self.assertEqual(data, [])
        
    def test_Err_lower(self):
        """Testing err lower and equal than threshold"""
        
        self.test_logger.err(0, "test")
        self.test_logger.err(1, "test")
        
        self.test_errfile.seek(0)
        data = [line.split(": ")[1] for line in self.test_errfile.readlines()]
        
        self.assertEqual(data, ['test\n', 'test\n'])
        
    def test_Err_higher(self):
        """Testing err higher than threshold (no message)"""
        
        self.test_logger.err(2, "test")
        self.test_errfile.seek(0)
        data = [line.split(": ")[1] for line in self.test_errfile.readlines()]
        self.assertEqual(data, [])
        

class TestFastaFile(unittest.TestCase):
    #To verify a seqObject
    def setUp(self):
        """Open the sequence with GClib.Utility.FastaFile"""
        
        #Open the sequence with fastafile module
        self.test_seqObj = GClib.Utility.FastaFile("chr21.fa.gz")
    
    def test_GetNextSeq(self):
        """Testing GetNextSeq returns a Bio.SeqRecord.SeqRecord object"""
        self.assertEqual(type(self.test_seqObj.GetNextSeq()), Bio.SeqRecord.SeqRecord)
    
    def test_GetNextSeqNone(self):
        """Testing GetNextSeq returns None (no more sequences in the file)"""
        
        #Take the first sequence
        self.test_seqObj.GetNextSeq()
        
        #No sequence left. Verify
        self.assertEqual(self.test_seqObj.GetNextSeq(), None)
        
    def test_GetSeqbyID(self):
        """Testing getting sequence by ID"""
        
        seqRecord = self.test_seqObj.GetSeqbyID("chr21")
        self.assertEqual(type(seqRecord), Bio.SeqRecord.SeqRecord)
        self.assertEqual(seqRecord.name, "chr21")
        
    def test_IterSeqs(self):
        """Testing IterSeqs"""
        
        self.assertEqual(types.GeneratorType, type(self.test_seqObj.IterSeqs()))
    

if __name__ == "__main__":
    unittest.main()
    