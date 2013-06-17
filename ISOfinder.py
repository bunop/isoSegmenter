#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""

Created on Wed Jun 12 16:31:27 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

Main program for Isochore Definition

"""

import os
import argparse

#Modules for dealing with GC content and graph
import GClib
import GClib.Graphs
import GClib.Utility
import GClib.Elements

parser = argparse.ArgumentParser(description='Find Isochores in sequences')
parser.add_argument('-i', '--infile', type=str, required=True, help="Input Fasta File (also compressed)")
parser.add_argument('-o', '--outfile', type=str, required=True, help="Output isochore CSV files")
parser.add_argument('-g', '--graphfile', type=str, required=False, help="Output graph filename (PNG)")
parser.add_argument('-v', '--verbosity', type=int, required=False, default=GClib.logger.threshold, help="Verbosity level")
args = parser.parse_args()

#TODO: Setting windows_size
#TODO: Write windows in file
#TODO: Change GAP tolerance
#TODO: Setting sequence start and end
#TODO: Draw chr name
#TODO: Draw Genome Reasearch 2006 isochore profile
#TODO: Draw window
#TODO: Write GAP CSV file
#TODO: Change isochore class boundaries
#TODO: Switch to Isochore Profile and Isochore Rectangle Boxes
#TODO: Calculating a determined chromosome from a multi fasta file
#TODO: set the possibility to bypass the isochore file creation (only a graph)

if __name__ == "__main__":
    #verify verbosity level
    if args.verbosity != GClib.logger.threshold:
        #setting user defined threshold of verbosity
        GClib.logger.threshold = args.verbosity
    
    #Chromosome istance will not Dump isochore if file exist. So I can verify this 
    #before reading fasta file. Outfile is a required option
    if args.outfile != None and os.path.exists(args.outfile):
        raise Exception, "file %s esists!!!" %(args.outfile)
    
    #Checking for graph file existance
    if args.graphfile != None and os.path.exists(args.graphfile):
        raise Exception, "file %s esists!!!" %(args.graphfile)
    
    #Open the sequence file
    FastaFile = GClib.Utility.FastaFile(args.infile)
    
    #Seq Record object
    seqRecord = None
    
    #TODO: Deal with more than 1 sequence in input file
    if FastaFile.n_of_sequences == 0:
        raise Exception, "No sequence found"
        
    elif FastaFile.n_of_sequences == 1:
        seqRecord = FastaFile.GetNextSeq()
        
    else:
        raise Exception, "Cannot handle %s sequences" %(FastaFile.n_of_sequences)
        
    #Instantiating Chromosome Class with seqRecord object (gaps are determined automatically)
    Chrom = GClib.Elements.Chromosome(seqRecord)
    
    #Segmenting Sequence in windows
    #TODO: Here I can call ValueWindows with different windows sizes, or sequence coordinates
    Chrom.ValueWindows()
    
    #TODO: Here i can write windows CSV and PNG file
    
    #Finding Isochores
    Chrom.FindIsochores()
    
    #Writing Isochores in file
    Chrom.DumpIsochores(args.outfile)
    
    #Instantiating graph if it is necessary
    if args.graphfile != None:
        #Instantiating DrawChromosome Class
        Graph = GClib.Graphs.DrawChromosome()
        
        #Fixing appropriate values
        Graph.SetSequenceLength(Chrom.size)
        Graph.InitPicture()
        Graph.SetHorizontalLines([37, 41, 46, 53])
        Graph.SetColorsList(colorbyclass=True)
        Graph.DrawIsochoreRectangles(isochores=Chrom.isochores)
        Graph.DrawLegend()
        Graph.FinishPicture(drawlabels=False)
        Graph.EnlargeLabels()
        Graph.SaveFigure(args.graphfile)
    
    
