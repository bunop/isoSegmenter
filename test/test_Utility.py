# -*- coding: utf-8 -*-
"""

    Copyright (C) 2013-2021 ITB - CNR

    This file is part of isoSegmenter.

    isoSegmenter is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    isoSegmenter is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with isoSegmenter.  If not, see <http://www.gnu.org/licenses/>.


If you use isoSegmenter in your work, please cite this manuscript:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into
    Isochores. Evolutionary Bioinformatics. 2015;11:253-261.
    doi:10.4137/EBO.S27693

Created on Fri May 31 16:33:03 2013

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>

A test module for Utility library

"""

import os
import Bio
import types
import tempfile
import unittest

import GClib.Utility
import GClib

# getting module path
module_path = os.path.dirname(__file__)


class TestFastaFile(unittest.TestCase):
    # To verify a seqObject
    def setUp(self):
        """Open the sequence with GClib.Utility.FastaFile"""

        # Open the sequence with fastafile module
        self.test_seqObj = GClib.Utility.FastaFile(
            os.path.join(module_path, "chr21.fa.gz"))

    def test_GetNextSeq(self):
        """Testing GetNextSeq returns a Bio.SeqRecord.SeqRecord object"""
        self.assertEqual(
            type(
                self.test_seqObj.GetNextSeq()),
            Bio.SeqRecord.SeqRecord)

    def test_GetNextSeqNone(self):
        """Testing GetNextSeq returns None (no more sequences in the file)"""

        # Take the first sequence
        self.test_seqObj.GetNextSeq()

        # No sequence left. Verify
        self.assertEqual(self.test_seqObj.GetNextSeq(), None)

    def test_GetSeqbyID(self):
        """Testing getting sequence by ID"""

        seqRecord = self.test_seqObj.GetSeqbyID("chr21")
        self.assertEqual(type(seqRecord), Bio.SeqRecord.SeqRecord)
        self.assertEqual(seqRecord.name, "chr21")

    def test_IterSeqs(self):
        """Testing IterSeqs"""

        self.assertEqual(
            types.GeneratorType, type(
                self.test_seqObj.IterSeqs()))


class TestFileExists(unittest.TestCase):
    # To verify a seqObject
    def setUp(self):
        """Testing FileExists utility"""

        # create an empty temporary file
        fd, filename = tempfile.mkstemp()

        # my filename
        self.filename = filename

    def test_RaiseIfExists(self):
        """Raise Exception if file exists"""

        self.assertRaises(IOError, GClib.Utility.FileExists, self.filename)

    def test_RemoveIfExists(self):
        """Remove if file exists (no exception raised)"""

        GClib.Utility.FileExists(self.filename, remove_if_exists=True)

    def tearDown(self):
        """Removing tempfile"""

        if os.path.exists(self.filename):
            os.remove(self.filename)


if __name__ == "__main__":
    unittest.main()
