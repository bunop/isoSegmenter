# -*- coding: utf-8 -*-
"""

Created on Fri May 31 17:07:40 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

"""

import os
import csv
import sys
import tempfile
import unittest
import Bio.SeqUtils

from numpy.numarray import mlab

sys.path.append("..")

import GClib
import GClib.Utility
import GClib.Elements

class test_CalcClass(unittest.TestCase):
    def setUp(self):
        """A test case to verify class assignment"""
        
        self.GClevels = ((31, "L1"),
                         (37, "L1"),
                         (37.1, "L2"),
                         (40, "L2"),
                         (41, "L2"),
                         (41.1, "H1"),
                         (46, "H1"),
                         (46.1, "H2"),
                         (53, "H2"),
                         (53.1, "H3"),
                         (70, "H3")
                        )

    def test_CalcClass(self):
        """Testing class assignment"""
        
        for GClevel, Class in self.GClevels:
            result = GClib.Elements.CalcClass(GClevel)
            self.assertEqual(result, Class)

class test_Element(unittest.TestCase):
    def setUp(self):
        """Testing Element instantiation"""
        
        self._test_Element = GClib.Elements.Element()
        
    def test_SetSize(self):
        """Testing Element SetSize"""
        
        #testing setting start and no end and viceversa
        self.assertRaises(GClib.Elements.ElementError, self._test_Element.SetSize, None, 100000)
        self.assertRaises(GClib.Elements.ElementError, self._test_Element.SetSize, 100000, None)
        
        #cal SetSize and verify dimensions and start, end coordinates
        start = 100000
        end = 500000
        size = end-start
        
        #Normal SetSize call
        self._test_Element.SetSize(start,end)
        self.assertEqual(self._test_Element.start, start)
        self.assertEqual(self._test_Element.end, end)
        self.assertEqual(self._test_Element.size, size)
        
        #Verify coordinate inversion
        self._test_Element.SetSize(end, start)
        self.assertEqual(self._test_Element.start, start)
        self.assertEqual(self._test_Element.end, end)
        self.assertEqual(self._test_Element.size, size)
        
        #verify setsize with length == 0
        self.assertRaises(GClib.Elements.ElementError, self._test_Element.SetSize, 100000, 100000)

class test_Gap(unittest.TestCase):
    def setUp(self):
        """Testing Gap instantiation"""
        
        self._test_Gap = GClib.Elements.Gap()
        
    def test_SetSize(self):
        """Testing Gap SetSize"""
        
        #testing setting start and no end and viceversa
        self.assertRaises(GClib.Elements.ElementError, self._test_Gap.SetSize, None, 100000)
        self.assertRaises(GClib.Elements.ElementError, self._test_Gap.SetSize, 100000, None)
        
        #cal SetSize and verify dimensions and start, end coordinates
        start = 100000
        end = 500000
        size = end-start
        
        #Normal SetSize call
        self._test_Gap.SetSize(start,end)
        self.assertEqual(self._test_Gap.start, start)
        self.assertEqual(self._test_Gap.end, end)
        self.assertEqual(self._test_Gap.size, size)
        
        #Verify coordinate inversion
        self._test_Gap.SetSize(end, start)
        self.assertEqual(self._test_Gap.start, start)
        self.assertEqual(self._test_Gap.end, end)
        self.assertEqual(self._test_Gap.size, size)

class test_Window(unittest.TestCase):
    def setUp(self):
        """Testing Window instatiation"""
        
        self.GClevels = ((31, "L1"),
                         (37, "L1"),
                         (37.1, "L2"),
                         (40, "L2"),
                         (41, "L2"),
                         (41.1, "H1"),
                         (46, "H1"),
                         (46.1, "H2"),
                         (53, "H2"),
                         (53.1, "H3"),
                         (70, "H3")
                        )
        
        self._test_Window = GClib.Elements.Window()
        
    def test_SetSize(self):
        """Testing window SetSize"""
        
        #testing setting start and no end and viceversa
        self.assertRaises(GClib.Elements.ElementError, self._test_Window.SetSize, None, 100000)
        self.assertRaises(GClib.Elements.ElementError, self._test_Window.SetSize, 100000, None)
        
        #cal SetSize and verify dimensions and start, end coordinates
        start = 100000
        end = 500000
        size = end-start
        
        #Normal SetSize call
        self._test_Window.SetSize(start,end)
        self.assertEqual(self._test_Window.start, start)
        self.assertEqual(self._test_Window.end, end)
        self.assertEqual(self._test_Window.size, size)
        
        #Verify coordinate inversion
        self._test_Window.SetSize(end, start)
        self.assertEqual(self._test_Window.start, start)
        self.assertEqual(self._test_Window.end, end)
        self.assertEqual(self._test_Window.size, size)
        
    def test_SetGClevel(self):
        """Testing SetGClevel (assigning GC class)"""
        
        for GClevel, Class in self.GClevels:
            self._test_Window.SetGClevel(GClevel)
            self.assertEqual(self._test_Window.Class, Class)

class test_Isochore(unittest.TestCase):
    def setUp(self):
        """testing Isochore instantiation"""
        
        #Starting an Isochore element from a Window Element
        window = GClib.Elements.Window(start=200000, end=300000, GClevel=38)
        self._test_Isochore = GClib.Elements.Isochore(window=window)
        
    def test_AddWindow(self):
        """Testing adding window to an isochore"""
        
        #test Adding only a Window element
        element = GClib.Elements.Element(start=0,end=100000)
        gap = GClib.Elements.Gap(start=0,end=100000)
        isochore = GClib.Elements.Isochore(window=GClib.Elements.Window(start=200000, end=300000, GClevel=38))
        
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, element)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, gap)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, isochore)
        
        #Adding a simple window
        window = GClib.Elements.Window(start=300000, end=400000, GClevel=40)
        self._test_Isochore.AddWindow(window)
        self.assertEqual(self._test_Isochore.GClevels, [38,40])
        self.assertEqual(self._test_Isochore.start, 200000)
        self.assertEqual(self._test_Isochore.end, 400000)
        self.assertEqual(self._test_Isochore.avg_GClevel, round(mlab.mean([38,40]),6))
        self.assertEqual(self._test_Isochore.stddev_GClevel, round(mlab.std([38,40]),6))
        
        #Adding a window before this element (why one can do this?)
        window = GClib.Elements.Window(start=100000, end=200000, GClevel=39)
        self._test_Isochore.AddWindow(window)
        self.assertEqual(self._test_Isochore.GClevels, [38,40, 39])
        self.assertEqual(self._test_Isochore.start, 100000)
        self.assertEqual(self._test_Isochore.end, 400000)
        self.assertEqual(self._test_Isochore.avg_GClevel, round(mlab.mean([38,40,39]),6))
        self.assertEqual(self._test_Isochore.stddev_GClevel, round(mlab.std([38,40,39]),6))
        
        #Adding a non contiguous window
        window = GClib.Elements.Window(start=500000, end=600000, GClevel=39)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, window)
        
        #Adding a windows already in this isochore
        window = GClib.Elements.Window(start=100000, end=400000, GClevel=38)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, window)
        window = GClib.Elements.Window(start=250000, end=400000, GClevel=38)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, window)
        window = GClib.Elements.Window(start=0, end=150000, GClevel=38)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddWindow, window)
        
    def test_AddIsochore(self):
        """Testing adding isochore to current one"""
        
        #test Adding only a isochore element
        element = GClib.Elements.Element(start=0,end=100000)
        gap = GClib.Elements.Gap(start=0,end=100000)
        window = GClib.Elements.Window(start=200000, end=300000, GClevel=38)
        
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, element)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, gap)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, window)
        
        #defining an Isochore
        window = GClib.Elements.Window(start=300000, end=400000, GClevel=40)
        self._test_Isochore.AddWindow(window)
        
        window = GClib.Elements.Window(start=100000, end=200000, GClevel=39)
        self._test_Isochore.AddWindow(window)
        
        #Creating a new Isochore
        window = GClib.Elements.Window(start=400000, end=500000, GClevel=39)
        isochore = GClib.Elements.Isochore(window)
        window = GClib.Elements.Window(start=500000, end=600000, GClevel=38)
        isochore.AddWindow(window)
        window = GClib.Elements.Window(start=600000, end=700000, GClevel=40)
        isochore.AddWindow(window)
        
        #Add the new Isochore and verify values
        self._test_Isochore.AddIsochore(isochore)
        self.assertEqual(self._test_Isochore.GClevels, [38,40, 39, 39, 38, 40])
        self.assertEqual(self._test_Isochore.start, 100000)
        self.assertEqual(self._test_Isochore.end, 700000)
        self.assertEqual(self._test_Isochore.avg_GClevel, round(mlab.mean([38,40,39, 39, 38, 40]),6))
        self.assertEqual(self._test_Isochore.stddev_GClevel, round(mlab.std([38,40,39, 39, 38, 40]),6))
        
        #Adding a non contiguous isochore
        window = GClib.Elements.Window(start=1000000, end=1100000, GClevel=39)
        isochore = GClib.Elements.Isochore(window)
        window = GClib.Elements.Window(start=1100000, end=1200000, GClevel=39)
        isochore.AddWindow(window)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, isochore)

        #Adding and Isochore already in this Isochore
        window = GClib.Elements.Window(start=200000, end=300000, GClevel=38)
        isochore = GClib.Elements.Isochore(window)
        window = GClib.Elements.Window(start=300000, end=400000, GClevel=40)
        isochore.AddWindow(window)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, isochore)
        
        window = GClib.Elements.Window(start=0, end=100000, GClevel=38)
        isochore = GClib.Elements.Isochore(window)
        window = GClib.Elements.Window(start=100000, end=200000, GClevel=40)
        isochore.AddWindow(window)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, isochore)
        
        window = GClib.Elements.Window(start=600000, end=700000, GClevel=38)
        isochore = GClib.Elements.Isochore(window)
        window = GClib.Elements.Window(start=700000, end=800000, GClevel=40)
        isochore.AddWindow(window)
        self.assertRaises(GClib.Elements.IsochoreError, self._test_Isochore.AddIsochore, isochore)

class test_Chromosome(unittest.TestCase):
    #maybe loads SeqRecord object once
    fastafile = GClib.Utility.FastaFile("chr21.fa.gz")
    seqRecord = fastafile.GetNextSeq()
    
    #Open a local Gap CSV file, obtained by the precendent version of this software
    old_filein = open("old_chr21_hg19_gaps.csv", "rU")
    old_csvin = csv.reader(old_filein)
    
    #Discarding headers
    old_csvin.next()
    
    #Here's the list of gaps tests
    old_gaps = []
    
    #the coordinatees are 0 based in the old gap file version
    for line in old_csvin:
        start, end, size = [int(col) for col in line]
        
        old_gaps += [GClib.Elements.Gap(start=start, end=end)]
        
    old_filein.close()
    
    #Open a local Windows CSV file, obtained by the evaluation study of this program
    test_filein = open("test_windows_chr21.csv", "rU")
    test_csvin = csv.reader(test_filein)
    
    #Discarding headers
    test_csvin.next()
    
    #Here's the list of windows tests
    test_windows = []
    
    #the coordinatees are 1 based and have to be fixed
    for line in test_csvin:
        start, end, size, Class, GClevel = int(line[0]), int(line[1]), int(line[2]), line[3], line[4]
        
        if Class == "gap":
            gap = GClib.Elements.Gap()
            gap.start = start - 1 #beware to the 1 based coordinates
            gap.end = end
            gap.size = size
            gap.Class = Class
            
            #add this gap to window list
            test_windows += [gap]
            
        else:
            window = GClib.Elements.Window()
            window.start = start - 1 #beware to the 1 based coordinates
            window.end = end
            window.size = size
            window.Class = Class
            window.GClevel = float(GClevel)
        
            #Add this window to window list
            test_windows += [window]
    
    #close windows file
    test_filein.close()
    
    #Opena a local Isochore file, obtained bt the evaluation study of this program
    test_filein = open("test_isochores3_chr21.csv", "rU")
    test_csvin = csv.reader(test_filein)
    
    #Discarding headers
    test_csvin.next()
    
    #here's the list of isocores test
    test_isochores = []
    
    #the coordinatees are 1 based and have to be fixed
    for line in test_csvin:
        start, end, size, Class, avg_GClevel, stddev_GClevel = int(line[0]), int(line[1]), int(line[2]), line[3], line[4], line[5]
        
        if Class == "gap":
            gap = GClib.Elements.Gap()
            gap.start = start - 1 #beware to the 1 based coordinates
            gap.end = end
            gap.size = size
            gap.Class = Class
            
            #add this gap to isochore list
            test_isochores += [gap]
            
        else:
            isochore = GClib.Elements.Isochore()
            isochore.start = start - 1 #beware to the 1 based coordinate
            isochore.end = end
            isochore.size = size
            isochore.Class = Class
            isochore.avg_GClevel = float(avg_GClevel)
            isochore.stddev_GClevel = float(stddev_GClevel)
            
            #finally I have to set the GC_levels attribute class, by reading those values from
            #windows that compone isochores
            for window in test_windows:
                #pay attention to coordinates
                if window.start >= start-1 and window.end <= end:
                    #print "%s added to %s" %(window, isochore)
                    isochore.GClevels += [window.GClevel]
            
            #add this isochore to isochore list
            test_isochores += [isochore]
            
            #debug
            #break
        
    #close isochores file
    test_filein.close()
    
    #In order to instantiate this class from a (i)python terminal
    def runTest(self):
        pass
    
    def setUp(self):
        """Testing chromosome instantiation"""
        
        self._test_Chromosome = GClib.Elements.Chromosome(self.seqRecord)
        
        #Shut down logger
        self._default_threshold = GClib.logger.threshold
        GClib.logger.threshold = 0
        
    def tearDown(self):
        """A method called after each test"""
        
        #restart logger with the default verbosity
        GClib.logger.threshold = self._default_threshold
        
    def test_Scan4Gaps(self):
        """Testing Scan4Gaps"""
        
        new_gaps = self._test_Chromosome.gaps
        
        #Now testing gaps legth
        self.assertEqual(len(self.old_gaps), len(new_gaps))
        
        #And now thes all the gap element
        for i in range(len(self.old_gaps)):
            new_gap = new_gaps[i]
            test_gap = self.old_gaps[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_gap, test_gap)
        
    
    def test_LoadDumpGaps(self):
        """Testing Dump and Load Gaps"""
        
        testfile = tempfile.mktemp()
        
        #Passing file name to Dump and load Gap
        self._test_Chromosome.DumpGaps(outfile=testfile)
        
        chromosome = GClib.Elements.Chromosome()
        chromosome.LoadGaps(infile=testfile)
        
        #Now testing gaps legth
        self.assertEqual(len(self.old_gaps), len(chromosome.gaps))
        
        #And now thes all the gap element
        for i in range(len(self.old_gaps)):
            new_gap = chromosome.gaps[i]
            test_gap = self.old_gaps[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_gap, test_gap)
        
        #deleting the old file
        os.remove(testfile)
        
        #Now passing an open file handle
        outfile = open(testfile, "w")
        self._test_Chromosome.DumpGaps(outfile=outfile)
        
        #Closing file
        outfile.close()
        
        #Open a file for reading
        infile = open(testfile, "rU")
        chromosome = GClib.Elements.Chromosome()
        chromosome.LoadGaps(infile=infile)
        
        #closing inputfile
        infile.close()
        
        #Now testing gaps legth
        self.assertEqual(len(self.old_gaps), len(chromosome.gaps))
        
        #And now thes all the gap element
        for i in range(len(self.old_gaps)):
            new_gap = chromosome.gaps[i]
            test_gap = self.old_gaps[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_gap, test_gap)
            
        #remove testfile
        os.remove(testfile)
        
    def test_ValueWindows(self):
        """Testing Value windows"""
        
        self._test_Chromosome.ValueWindows()
        new_windows = self._test_Chromosome.windows
        
        #Now testing windows legth
        self.assertEqual(len(self.test_windows), len(new_windows))
        
        #And now thes all the gap element
        for i in range(len(self.test_windows)):
            new_window = new_windows[i]
            test_window = self.test_windows[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_window, test_window)
    
    def test_LoadDumpWindows(self):
        """Testing Dump and Load Windows"""
        
        #calculating windows
        self._test_Chromosome.ValueWindows()
        
        #A test file for reading and writing
        testfile = tempfile.mktemp()
        
        #Passing file name to Dump and load Gap
        self._test_Chromosome.DumpWindows(outfile=testfile)
        
        chromosome = GClib.Elements.Chromosome()
        chromosome.LoadWindows(infile=testfile)
        
        #Now testing windows legth
        self.assertEqual(len(self.test_windows), len(chromosome.windows))
        
        #And now thes all the windows element
        for i in range(len(self.test_windows)):
            new_window = chromosome.windows[i]
            test_window = self.test_windows[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_window, test_window)
        
        #deleting the old file
        os.remove(testfile)

    def test_ValueWindowsWithDifferentStartEnd(self):
        """Testing Value windows with different start - end coordinates"""
        
        #This is the position of first GAP in chromosome sequence start:1,end:9411193,size:9411193,Class:gap. So:
        gap0 = GClib.Elements.Gap(start=9000000,end=9411193)
        
        #Setting a different start point
        self._test_Chromosome.ValueWindows(From=9000000)
        
        #Asserting the first gap
        self.assertEqual(gap0, self._test_Chromosome.windows[0])
        
        #Now setting the start site immediately after the GAP. The coordinates are 0 based
        window0 = GClib.Elements.Window(start=9411193,end=9511193,GClevel=36.522)
        
        self._test_Chromosome.ValueWindows(From=9411193)
        
        #Asserting first window
        self.assertEqual(window0, self._test_Chromosome.windows[0])
        
        #Starting windows from a different position
        window1 = GClib.Elements.Window(start=15000000,end=15100000,GClevel=42.56)
        
        self._test_Chromosome.ValueWindows(From=15000000)
        
        #Asserting casual window
        self.assertEqual(window1, self._test_Chromosome.windows[0])
        
        #Now testing windows until a determined position
        last_gap = GClib.Elements.Gap(start=48119895,end=48120000)
        
        #limiting analysis to a position inside a GAP
        self._test_Chromosome.ValueWindows(To=48120000)
        
        self.assertEqual(last_gap, self._test_Chromosome.windows[-1])
        
        #Limiting analysis immediatly after a GAP
        last_windows = GClib.Elements.Window(start=48082664,end=48119895,GClevel=43.482582)
        self._test_Chromosome.ValueWindows(To=48119895)
        
        self.assertEqual(last_windows, self._test_Chromosome.windows[-1])
        
        #Now evaluate a segment completely inside a GAP
        whole_gap = GClib.Elements.Gap(start=1000000,end=4000000)
        self._test_Chromosome.ValueWindows(From=1000000,To=4000000)
        
        #I want only a gap
        self.assertEqual(len(self._test_Chromosome.windows),1)
        self.assertEqual(self._test_Chromosome.windows[0], whole_gap)
        
        #And now we evaluate a region without gaps
        From, To = 20000000,21000000
        window_size = 100000
        
        #test for GC levels in this regions (look how it simple when no gap is considered ando supposing no isochore detection
        GC_values = [round(Bio.SeqUtils.GC(self.seqRecord[step:step+window_size].seq),6) for step in range(From, To, window_size)]
        
        #determining windows
        self._test_Chromosome.ValueWindows(From=From, To=To)
        
        #inspect windows positions and values:
        self.assertEqual(len(self._test_Chromosome.windows), len(GC_values))
        
        for i, window in enumerate(self._test_Chromosome.windows):
            start = From + i*window_size
            end = start + window_size
            self.assertEqual(window.start, start)
            self.assertEqual(window.end, end)
            self.assertEqual(window.GClevel, GC_values[i])
        
    def test_FindIsochores(self):
        """Testing FindIsochores"""
        
        self._test_Chromosome.ValueWindows()
        self._test_Chromosome.FindIsochores()
        
        new_isochores = self._test_Chromosome.isochores
        
        #Now testing windows legth
        self.assertEqual(len(self.test_isochores), len(new_isochores))
        
        #And now thes all the gap element
        for i in range(len(self.test_isochores)):
            new_isochore = new_isochores[i]
            test_isochore = self.test_isochores[i]
            
            #Assert equality with __eq__ method
            self.assertEqual(new_isochore, test_isochore)
            
    def test_LoadDumpIsochores(self):
        """Testing Dump and Load Isochores"""
        
        #calculating isochores
        self._test_Chromosome.ValueWindows()
        self._test_Chromosome.FindIsochores()
        
        #A test file for reading and writing
        testfile = tempfile.mktemp()
        
        #Passing file name to Dump and load Gap
        self._test_Chromosome.DumpIsochores(outfile=testfile)
        
        chromosome = GClib.Elements.Chromosome()
        chromosome.LoadIsochores(infile=testfile)
        
        #Now testing windows legth
        self.assertEqual(len(self.test_isochores), len(chromosome.isochores))
        
        #And now thes all the windows element
        for i in range(len(self.test_isochores)):
            new_isochore = chromosome.isochores[i]
            test_isochore = self.test_isochores[i]
            
            #the new load isochore leaks from GClevels attribute, since it isn't present
            #in CSV file. So I fix this term before comparison, if other class attributes
            #will be different, the test will fail
            if new_isochore.Class != "gap" and test_isochore.Class != "gap":
                new_isochore.GClevels = test_isochore.GClevels
            
            #Assert equality with __eq__ method
            self.assertEqual(new_isochore, test_isochore)
        
        #deleting the old file
        os.remove(testfile)

if __name__ == "__main__":
    unittest.main()
