#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""


    Copyright (C) 2013-2016 ITB - CNR

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

Created on Thu Jun 13 14:53:14 2013

@author: Paolo Cozzi <paolo.cozzi@ptp.it>

This is the main program to find isochores families in a user defined directory

"""

import os
import argparse

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

isoSegmenter  Copyright (C) 2013-2016 ITB - CNR
This program comes with ABSOLUTELY NO WARRANTY; for details type `isoFamily.py --help'.
This is free software, and you are welcome to redistribute it
under certain conditions; show LICENSE.md for more details.

"""

parser = argparse.ArgumentParser(description='Find Isochores Families in a user defined directory', epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--indir', type=str, required=True, help="Input directory in which isochore have been calculated")
parser.add_argument('-o', '--outfile', type=str, required=True, help="Output families CSV files")
parser.add_argument('-g', '--graphfile', type=str, required=False, help="Output graph filename (PNG)")
parser.add_argument('-r', '--regexp', type=str, required=False, default=".csv", help="pattern for isochore file search (default: '%(default)s')")
parser.add_argument('-v', '--verbose', action="count", required=False, default=0, help="verbose level")
parser.add_argument('--force_overwrite', action='store_true', default=False, help="Force overwrite")
parser.add_argument('--x_max', type=int, required=False, default=GClib.Graphs.GRAPH_GC_MAX, help="Set X max value in graph (default: '%(default)s')")
parser.add_argument('--x_min', type=int, required=False, default=GClib.Graphs.GRAPH_GC_MIN, help="Set X min value in graph (default: '%(default)s')")
args = parser.parse_args()

#TODO: Adding option for chainging Xmax and Xmin values
#TODO: Adding option for change bin dimension
#TODO: Change label dimensions
#TODO: Adding support for user axis definition
#TODO: Adding support for adding title
#TODO: change figure resolution

if __name__ == "__main__":
    #print out notice
    GClib.logger.err(0, notice)

    #verify verbosity level
    if args.verbose != GClib.logger.threshold:
        #setting user defined threshold of verbosity
        GClib.logger.threshold = args.verbose

    #chceking for csv existance
    GClib.Utility.FileExists(args.outfile, remove_if_exists=args.force_overwrite)

    #Checking for graph file existance
    GClib.Utility.FileExists(args.graphfile, remove_if_exists=args.force_overwrite)

    #instantiate a families element
    families = GClib.Elements.Families()

    #scanning for isochore files in user defined directory
    families.Scan4Files(args.indir, pattern=args.regexp)

    #Group isochore relying their GClevels
    #TODO: here i can set Xmax and Xmin values and bin size dimension
    families.GroupByIsochores()

    #Dump families in a CSV file
    families.DumpFamilies(outfile=args.outfile)

    #Instantiating graph if it is necessary
    if args.graphfile != None:
        #Instantiating DrawFamilies Class
        graph = GClib.Graphs.DrawFamilies(families)

        #Setting Grids and axis Labels
        graph.DrawAxisLabels()
        graph.DrawGrid()

        #Override Xmax and Xmin values
        axes = [args.x_min, args.x_max, 0, 0]
        graph.SetAxisLimits(axis=axes)

        #save image
        #TODO: change figure quality
        graph.SaveFigure(filename=args.graphfile)
