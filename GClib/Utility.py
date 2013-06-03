# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:07:25 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

A class in which we can find all utilities to deal with isochores

"""

import os
import sys
import gzip
import time
import GClib
import Bio.SeqIO

#for logging messages
class Logger:
    """To set up a level of verbosity of an application"""
    def __init__(self,threshold=1,outfile=sys.stdout):
        self.threshold = threshold
        self.outfile = outfile
    
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
        self.log(level,message,sys.stderr)

#Setting a Logger for this module
Log = Logger(thresold=5)

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
        
    def old_FindGaps(self, sequenceObj):
        """Permette di analizzare una sequenza Bio::Seq
        per capire dove sono posizionati i GAP"""
        
        #queste variabili potrebbero tornarmi utili
        self.dna_alphabet = ['t', 'a', 'c', 'g', 'n']
        self.gaps_start = []
        self.gaps_end = []
        self.gaps_size = []
        self.n_of_gaps = 0
        
        #inizializzo delle variabili importanti
        gap_flag = 0
        gap_size = 0
        gap_start = 0
        i = 0
        
        #ok ciclo lungo la sequenza Bio::Seq alla ricerca di Gap
        for letter in sequenceObj:
            #attenzione ai caratteri. Diventeranno tutti minuscoli
            letter = letter.lower()
            
            #controlliamo che letter sia nell'alfabeto
            if letter not in self.dna_alphabet:
                Log.err(1, "Base sconosciuta %s in posizione %s" %(letter, i))
        
            #controlliamo lo stato dei gap
            if letter == 'n' and gap_flag == 0:
                Log.log(3, "Gap aperto in posizione %s" %(i))
                
                gap_flag = 1
                gap_size = 1
                gap_start = i
        
            elif gap_flag == 1 and letter != 'n':
                Log.log(3, "Gap chiuso in posizione %s" %(i))
                
                #Devo memorizzare questi dati
                self.gaps_start += [gap_start]
                self.gaps_end += [i]
                self.gaps_size += [gap_size]
                
                #aggiusto il GAP
                gap_flag = 0
        
            elif gap_flag == 1:
                #incremento la dimensione del GAP
                gap_size += 1
                
            #IMPORTANTISSIMO: aumentare il contatore i
            i += 1
        
        #Ok. Se finisco la sequenza e il gap_flag è ancora aperto?
        if gap_flag == 1:
            Log.log(4, "La sequenza termina con un gap")
            
            #questa volta devo chiudere il GAP
            self.gaps_start += [gap_start]
            self.gaps_end += [i]
            self.gaps_size += [gap_size]
        
        #OK adesso il numero di gaps_start, gaps_end e gaps_size
        #deve essere lo stesso
        if len(self.gaps_start) != len(self.gaps_end) or len(self.gaps_start) != len(self.gaps_size):
            raise Exception, "C'è stato un problema nella lettura dei gaps"
        
        else:
            self.n_of_gaps = len(self.gaps_start)
        
        #debug
        if self.n_of_gaps != 0:
            Log.log(2, "Sono stati trovati %s gap" %(self.n_of_gaps))
        
        #Ok, se arrivo qua ho letto le informazioni di tutti  gap
        return
        