#! /usr/bin/env python
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

Created on Wed Jun 12 16:31:27 2013

@author: Paolo Cozzi <paolo.cozzi@ptp.it>

Main program for Isochore Definition

"""

import argparse

#Modules for dealing with GC content and graph
import GClib
import GClib.Graphs
import GClib.Utility
import GClib.Elements

# Add epilog on bottom of help message
epilog ="""

If you use isoSegmenter in your work, please cite this manuscript:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into Isochores.
    Evolutionary Bioinformatics. 2015;11:253-261. doi:10.4137/EBO.S27693

    """

notice = """

isoSegmenter  Copyright (C) 2013-2015 ITB - CNR
This program comes with ABSOLUTELY NO WARRANTY; for details type `isoSegmenter --help'.
This is free software, and you are welcome to redistribute it
under certain conditions; show LICENSE.md for more details.

"""

parser = argparse.ArgumentParser(description='Find Isochores in sequences', epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--infile', type=str, required=True, help="Input Fasta File (even compressed)")
parser.add_argument('-o', '--outfile', type=str, required=False, help="Output isochores CSV file")
parser.add_argument('-g', '--graphfile', type=str, required=False, help="Output graph filename (PNG)")
parser.add_argument('-b', '--barfile', type=str, required=False, help="Output bar graph filename (PNG)")
parser.add_argument('-w', '--windowfile', type=str, required=False, help="Output windows CSV file")
parser.add_argument('-v', '--verbose', action="count", required=False, default=0, help="Verbosity level")
parser.add_argument('--windowgraph', type=str, required=False, help="Output windows Graph file")
parser.add_argument('--draw_legend', action='store_true', help="Draw legend on the right side of the image")
parser.add_argument('--force_overwrite', action='store_true', help="Force overwrite")
parser.add_argument('--sequence_start', type=float, required=False, default=1, help="Start segmentation from this position (1-based coordinates)")
parser.add_argument('--max_length', type=float, required=False, default=None, help="Scan for isochores until for this dimension in bp")
parser.add_argument('--draw_chname', type=str, required=False, default=None, help="Draw chromosome name in figure")
parser.add_argument('--window_size', type=int, required=False, default=GClib.WINDOW_SIZE, help="Set window size in bp (default: '%(default)s')")
parser.add_argument('--y_max', type=int, required=False, default=GClib.Graphs.GRAPH_GC_MAX, help="Set max value in graph (default: '%(default)s')")
parser.add_argument('--y_min', type=int, required=False, default=GClib.Graphs.GRAPH_GC_MIN, help="Set min value in graph (default: '%(default)s')")
parser.add_argument('--isochore_min_size', type=int, required=False, default=GClib.ISO_MIN_SIZE, help="Set how many windows an isochore need to have (default: '%(default)s')")
args = parser.parse_args()

#debug
#print args

#TODO: Change GAP tolerance
#TODO: Setting sequence start and end
#TODO: Write GAP CSV file
#TODO: Change isochore class boundaries
#TODO: Switch to Isochore Profile and Isochore Rectangle Boxes
#TODO: Calculating a determined chromosome from a multi fasta file


if __name__ == "__main__":
    #To continue work, I need almost one file to write
    if args.outfile == None and args.graphfile == None and args.barfile == None:
        raise Exception, "You must specify an output isochore file while calling this program, by graphfile, barfile or outfile option"

    #print out notice
    GClib.logger.err(0, notice)

    #verify verbosity level
    if args.verbose != GClib.logger.threshold:
        #setting user defined threshold of verbosity
        GClib.logger.threshold = args.verbose

    #Chromosome istance will not Dump isochore if file exist. So I can verify this
    #before reading fasta file. Outfile is a required option
    GClib.Utility.FileExists(args.outfile, remove_if_exists=args.force_overwrite)

    #Checking for graph file existance
    GClib.Utility.FileExists(args.graphfile, remove_if_exists=args.force_overwrite)

    #Checking for graph file existance
    GClib.Utility.FileExists(args.barfile, remove_if_exists=args.force_overwrite)

    #Checking for window file existance
    GClib.Utility.FileExists(args.windowfile, remove_if_exists=args.force_overwrite)

    #Checking for window graph file existance
    GClib.Utility.FileExists(args.windowgraph, remove_if_exists=args.force_overwrite)

    #sequence_start can't be negative
    if args.sequence_start <= 0:
        raise Exception, "Sequence start must be 1-based and > 0"

    #Internal coordinates are 0-based, not 1-based
    args.sequence_start = int(args.sequence_start) -1

    #To is the position in which isochore calculation ends. None will be threated
    #as chromosome end position
    To = None

    #Checking user max_length of sequence analysis
    if args.max_length != None:
        #coerce max_length into integer
        args.max_length = int(args.max_length)

        #Setting max_length as To coordinates, starting from sequence start
        To = args.sequence_start + args.max_length

    #Evaluaing isochore min size
    if args.isochore_min_size != GClib.ISO_MIN_SIZE:
        GClib.ISO_MIN_SIZE = args.isochore_min_size

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
    if args.window_size == GClib.WINDOW_SIZE:
        #This call only to print the warning string when we use the default window size
        Chrom.ValueWindows(From=args.sequence_start, To=To)

    else:
        #Call valuewindos with user defined window size
        Chrom.ValueWindows(window_size=args.window_size, From=args.sequence_start, To=To)

    #Writing windows in a file (if I need it)
    if args.windowfile != None:
        Chrom.DumpWindows(args.windowfile)

    #Writing the window graph file, if is needed
    if args.windowgraph != None:
        #Instantiating DrawChromosome Class. Look at sequence start (0-based sequence start, this has been fixed in the top of this main block)
        Graph = GClib.Graphs.DrawChromosome(sequence_start=args.sequence_start)

        #beware user defined min and max values
        if args.y_max != None or args.y_min != None:
            #If only one value is defined by the user, get the othert value
            if args.y_max == None:
                args.y_max = Graph.y_max

            if args.y_min == None:
                args.y_min = Graph.y_min

            #Set max and min values
            Graph.SetMinMaxValues(args.y_min, args.y_max)

        #Fixing appropriate values
        if args.max_length != None:
            #SetSequencelength needs the To position (the absolute end position)
            Graph.SetSequenceLength(To)

        else:
            Graph.SetSequenceLength(Chrom.size)

        Graph.InitPicture()
        Graph.SetHorizontalLines([37, 41, 46, 53])
        Graph.SetColorsList(colorbyclass=True)

        #Draw the correct values
        Graph.DrawWindowRectangles(windows=Chrom.windows)

        #Draw legend or not
        if args.draw_legend == True:
            Graph.DrawLegend()

        #Draw ChName
        if args.draw_chname != None:
            Graph.DrawChName(args.draw_chname)

        #Finishing picture
        Graph.FinishPicture(drawlabels=False)
        Graph.EnlargeLabels()
        Graph.SaveFigure(args.windowgraph)

    #Finding Isochores. This program tries to segmenting genome into isochores, and so this calculation is always done
    Chrom.FindIsochores()

    if args.outfile != None:
        #Writing Isochores in file
        Chrom.DumpIsochores(args.outfile)

    #Instantiating graph if it is necessary
    if args.graphfile != None:
        #Instantiating DrawChromosome Class. Look at sequence start (0-based sequence start, this has been fixed in the top of this main block)
        Graph = GClib.Graphs.DrawChromosome(sequence_start=args.sequence_start)

        #beware user defined min and max values
        if args.y_max != None or args.y_min != None:
            #If only one value is defined by the user, get the othert value
            if args.y_max == None:
                args.y_max = Graph.y_max

            if args.y_min == None:
                args.y_min = Graph.y_min

            #Set max and min values
            Graph.SetMinMaxValues(args.y_min, args.y_max)

        #Fixing appropriate values
        if args.max_length != None:
            #SetSequencelength needs the To position (the absolute end position)
            Graph.SetSequenceLength(To)

        else:
            Graph.SetSequenceLength(Chrom.size)

        Graph.InitPicture()
        Graph.SetHorizontalLines([37, 41, 46, 53])
        Graph.SetColorsList(colorbyclass=True)

        #Draw the correct values
        Graph.DrawIsochoreRectangles(isochores=Chrom.isochores)

        #Draw legend or not
        if args.draw_legend == True:
            Graph.DrawLegend()

        #Draw ChName
        if args.draw_chname != None:
            Graph.DrawChName(args.draw_chname)

        Graph.FinishPicture(drawlabels=False)
        Graph.EnlargeLabels()
        Graph.SaveFigure(args.graphfile)


    #Create bar graph isocore grap (as Schmidt and Frishman 2008) if it is necessary
    if args.barfile != None:
        #Instantiating DrawBarChromosome Class. Look at sequence start (0-based sequence start, this has been fixed in the top of this main block)
        Graph = GClib.Graphs.DrawBarChromosome(sequence_start=args.sequence_start)

        #there are no min and max values in this graph style

        #Fixing appropriate values
        if args.max_length != None:
            #SetSequencelength needs the To position (the absolute end position)
            Graph.SetSequenceLength(To)

        else:
            Graph.SetSequenceLength(Chrom.size)

        Graph.InitPicture()
        # No horyzontal lines
        Graph.SetColorsList(colorbyclass=True)

        #Draw the correct values
        Graph.DrawIsochoreRectangles(isochores=Chrom.isochores)

        #Draw legend or not
        if args.draw_legend == True:
            Graph.DrawLegend()

        #Draw ChName
        if args.draw_chname != None:
            Graph.DrawChName(args.draw_chname)

        Graph.FinishPicture(drawlabels=False)
        Graph.EnlargeLabels()
        Graph.SaveFigure(args.barfile)
