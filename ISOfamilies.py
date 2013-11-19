#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""


    Copyright (C) 2013 Paolo Cozzi

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


Created on Thu Jun 13 14:53:14 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

This is the main program to find isochores families in a user defined directory

"""

import os
import argparse

import GClib
import GClib.Graphs
import GClib.Elements

parser = argparse.ArgumentParser(description='Find Isochores Families in a user defined directory')
parser.add_argument('-i', '--indir', type=str, required=True, help="Input directory in which isochore have been calculated")
parser.add_argument('-o', '--outfile', type=str, required=True, help="Output families CSV files")
parser.add_argument('-g', '--graphfile', type=str, required=False, help="Output graph filename (PNG)")
parser.add_argument('-r', '--regexp', type=str, required=False, default=".csv", help="pattern for isochore file search (default: '%(default)s')")
parser.add_argument('-v', '--verbosity', type=int, required=False, default=GClib.logger.threshold, help="Verbosity level")
args = parser.parse_args()

#TODO: Adding option for chainging Xmax and Xmin values
#TODO: Adding option for change bin dimension
#TODO: Change label dimensions
#TODO: Adding support for user axis definition
#TODO: Adding support for adding title
#TODO: change figure resolution

if __name__ == "__main__":
    #verify verbosity level
    if args.verbosity != GClib.logger.threshold:
        #setting user defined threshold of verbosity
        GClib.logger.threshold = args.verbosity
    
    #Chromosome istance will not Dump isochore if file exist. So I can verify this 
    #before reading fasta file. Outfile is a required option
    if args.outfile != None and os.path.exists(args.outfile):
        raise Exception, "file %s exists!!!" %(args.outfile)
    
    #Checking for graph file existance
    if args.graphfile != None and os.path.exists(args.graphfile):
        raise Exception, "file %s exists!!!" %(args.graphfile)

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
        
        #save image
        #TODO: change figure quality
        graph.SaveFigure(filename=args.graphfile)
        
        