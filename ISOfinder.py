#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""


    Copyright (C) 2013 ITB - CNR

    This file is part of ISOfinder.

    ISOfinder is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ISOfinder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ISOfinder.  If not, see <http://www.gnu.org/licenses/>.


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
parser.add_argument('-i', '--infile', type=str, required=True, help="Input Fasta File (even compressed)")
parser.add_argument('-o', '--outfile', type=str, required=False, help="Output isochores CSV file")
parser.add_argument('-g', '--graphfile', type=str, required=False, help="Output graph filename (PNG)")
parser.add_argument('-w', '--windowfile', type=str, required=False, help="Output windows CSV file")
parser.add_argument('-v', '--verbose', action="count", required=False, default=0, help="Verbosity level")
parser.add_argument('--draw_legend', action='store_true', help="Draw legend on the right side of the image")
parser.add_argument('--force_overwrite', action='store_true', help="Force overwrite")
parser.add_argument('--sequence_start', type=int, required=False, default=1, help="start segmentation from this position (1-based coordinates)")
args = parser.parse_args()

#debug
#print args

#TODO: Setting windows_size
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
    #To continue work, I need almost one file to write
    if args.outfile == None and args.graphfile == None:
        raise Exception, "You must specify an output isochore file while calling this program, by graphfile or outfile option"
        
    #verify verbosity level
    if args.verbose != GClib.logger.threshold:
        #setting user defined threshold of verbosity
        GClib.logger.threshold = args.verbose
    
    #Chromosome istance will not Dump isochore if file exist. So I can verify this 
    #before reading fasta file. Outfile is a required option
    if args.outfile != None and os.path.exists(args.outfile):
        if args.force_overwrite == False:
            raise Exception, "file %s exists!!!" %(args.outfile)
        else:
            #remove the file before calculation
            os.remove(args.outfile)
    
    #Checking for graph file existance
    if args.graphfile != None and os.path.exists(args.graphfile):
        if args.force_overwrite == False:
            raise Exception, "file %s exists!!!" %(args.graphfile)
        else:
            #remove the file before calculation
            os.remove(args.graphfile)
    
    #Checking for window file existance
    if args.windowfile != None and os.path.exists(args.windowfile):
        if args.force_overwrite == False:
            raise Exception, "file %s exists!!!" %(args.windowfile)
        else:
            #remove the file before calculation
            os.remove(args.windowfile)
        
    
    #sequence_start can't be negative
    if args.sequence_start <= 0:
        raise Exception, "Sequence start must be 1-based and > 0"
    
    #Internal coordinates are 0-based, not 1-based
    args.sequence_start -= 1
    
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
    Chrom.ValueWindows(From=args.sequence_start)
    
    #Writing windows in a file (if I need it)
    if args.windowfile != None:
        Chrom.DumpWindows(args.windowfile)
        
    #Finding Isochores. This program tries to segmenting genome into isochores, and so this calculation is always done
    Chrom.FindIsochores()
    
    if args.outfile != None:
        #Writing Isochores in file
        Chrom.DumpIsochores(args.outfile)
    
    #Instantiating graph if it is necessary
    if args.graphfile != None:
        #Instantiating DrawChromosome Class. Look at sequence start (0-based sequence start, this has been fixed in the top of this main block)
        Graph = GClib.Graphs.DrawChromosome(sequence_start=args.sequence_start)
        
        #Fixing appropriate values
        Graph.SetSequenceLength(Chrom.size)
        Graph.InitPicture()
        Graph.SetHorizontalLines([37, 41, 46, 53])
        Graph.SetColorsList(colorbyclass=True)
        
        #Grep the correct values
        Graph.DrawIsochoreRectangles(isochores=Chrom.isochores)
        
        #TODO: make an option to draw windows
        #To draw window instead of isochores. 
        #Graph.DrawWindowRectangles(windows=Chrom.windows)
        
        #Draw legend or not
        if args.draw_legend == True:
            Graph.DrawLegend()
        
        Graph.FinishPicture(drawlabels=False)
        Graph.EnlargeLabels()
        Graph.SaveFigure(args.graphfile)
    
    
