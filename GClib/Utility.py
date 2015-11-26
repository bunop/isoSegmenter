# -*- coding: utf-8 -*-
"""


    Copyright (C) 2013-2015 ITB - CNR

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

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into Isochores.
    Evolutionary Bioinformatics. 2015;11:253-261. doi:10.4137/EBO.S27693

Created on Fri May  3 15:07:25 2013

@author: Paolo Cozzi <paolo.cozzi@ptp.it>

A class in which we can find all utilities to deal with isochores

"""

import os
import sys
import gzip
import time
import GClib
import Bio.SeqIO

__author__ = "Paolo Cozzi <paolo.cozzi@ptp.it>"

from . import __copyright__, __license__, __version__

#for logging messages
class Logger:
    """To set up a level of verbosity of an application"""
    def __init__(self,threshold=1,outfile=sys.stdout,errfile=sys.stderr):
        self.threshold = threshold
        self.outfile = outfile
        self.errfile = errfile

    def log(self,level,message,outfile=None):
        '''Logs message if level lower than or equal to threshold.
        Message can even be an object - normal print functionality
        will be used'''

        #Se non specifico un oufile, allora prendo quello di default
        if outfile is None: outfile = self.outfile

        if level <= self.threshold:
            print >> outfile, "%s: %s" %(time.ctime(),message)
            #flushing outfile
            outfile.flush()

    def err(self,level,message):
        '''Logs message in stderr, instead of predefined file'''
        self.log(level,message,self.errfile)


#To deal with fasta files
class FastaFile:
    """A class to deal with fasta files"""

    def __init__(self, fasta_file=None):
        """To instantiate the class. You may give a fasta path (also compressed)"""

        self.last_idx = 0
        self.seqs_list = []
        self.seqs_ids = {}
        self.n_of_sequences = 0

        #Open a fasta file, if requested
        if fasta_file != None:
            self.Load(fasta_file)

    def Load(self,fasta_file):
        """Load a fasta file"""

        #debug
        GClib.logger.log(1, "Opening %s..." %(fasta_file))

        #verify the file extension
        extension = os.path.splitext(fasta_file)[1]

        if extension == '.gz':
            #Open handle with gzip
            fasta_fh = gzip.open(fasta_file,"rb")

        else:
            #open handle in universal mode
            fasta_fh = open(fasta_file,"rU")

        #Parsing sequences with Bio.SeqIO
        self.seqs_list = list(Bio.SeqIO.parse(fasta_fh, "fasta"))

        #How many sequences were read?
        self.n_of_sequences = len(self.seqs_list)

        #Which are sequences ids
        for idx, seq in enumerate(self.seqs_list):
            self.seqs_ids[seq.id] = idx

        #debug
        GClib.logger.log(1, "%s sequences read" %(self.n_of_sequences))

    def IterSeqs(self):
        """Iters through Bio.Seqs Objects"""

        for seq_obj in self.seq_list:
            yield seq_obj

    def GetNextSeq(self):
        """Give the next sequence"""

        if self.last_idx >= self.n_of_sequences:
            return None

        seq_obj = self.seqs_list[self.last_idx]
        self.last_idx += 1

        return seq_obj

    def GetSeqbyID(self, id):
        """Return a SeqObj by id"""

        #get the position in list
        idx = self.seqs_ids[id]

        #return seq obj
        return self.seqs_list[idx]


#a function to check file existance and remove file if needed
def FileExists(filename, remove_if_exists=False):
    """Testing for file existance and removing file if needed"""

    #return if filenames is None
    if filename == None:
        return

    if os.path.exists(filename):
        if remove_if_exists == False:
            raise IOError, "file %s exists!!!" %(filename)

        else:
            #remove the file before calculation
            GClib.logger.log(2, "file %s removed" %(filename))
            os.remove(filename)

    #this function return nothing is successful

#end of library
