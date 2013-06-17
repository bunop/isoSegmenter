# -*- coding: utf-8 -*-
"""

Created on Wed Jun  5 10:16:54 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

A module to test Graphs module

"""

import os
import sys
import StringIO
import tempfile
import unittest


sys.path.append("..")

import GClib
import GClib.Graphs
import GClib.Elements

__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

class test_BaseGraph(unittest.TestCase):
    def setUp(self):
        self._test_BaseGrap = GClib.Graphs.BaseGraph()
        
        #Setting known values to obtaining known results
        self._test_BaseGrap.scale = 30000 #17500 #higher values shrink images
        self._test_BaseGrap.border = 90 #the white space on the left and on the right of the figure
        self._test_BaseGrap.top = 70 #the upper space before the X axis
        self._test_BaseGrap.y = 385 #the height of the graphic area (in which isocore are printed, not image height)
        
    def test_SetMinMaxValues(self):
        """Testing set max and min values"""
        
        #passing Max and Min value instead of Min and Max
        self._test_BaseGrap.SetMinMaxValues(65,30)
        
        self.assertEqual(self._test_BaseGrap.y_min, 30)
        self.assertEqual(self._test_BaseGrap.y_max, 65)
        self.assertEqual(self._test_BaseGrap.py, 9)
        
    def test_SetSequenceLength(self):
        """Testing set sequence length"""
        
        #set sequence length to 1Mb
        self._test_BaseGrap.SetSequenceLength(1e6)
        self.assertEqual(self._test_BaseGrap.x, 213)
        
    def test_SetFontSize(self):
        """Testing font size over maximum values"""
        
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetFontSize, (5,))
        
    def test_InitPicture(self):
        """Testing InitPicture (GD have to work properly)"""
        
        #testing exception when no sequence length is provided
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.InitPicture)
        
        #Set the right sequence length
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.InitPicture()
        
        imagesize = self._test_BaseGrap.graph.size()
        self.assertEqual(imagesize, (213, 420))
    
    def test_SetHorizontalLines(self):
        """Testing SetHorizontaLines"""
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetColorsList)
        
        #correct call of SetHorizontalLines
        self._test_BaseGrap.SetMinMaxValues(65,30)
        
        #testing a wrong parameter
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetHorizontalLines, None)
        
        #Setting a list
        self._test_BaseGrap.SetHorizontalLines([37, 41, 46, 53])
        self.assertEqual(self._test_BaseGrap.h_lines, [37, 41, 46, 53])
        
        #Setting a integer
        self._test_BaseGrap.SetHorizontalLines(5)
        self.assertEqual(self._test_BaseGrap.n_of_h_lines, 5)
        
    
    def test_SetColorsList(self):
        """Testing color palette"""
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetColorsList)
        
        #testing methods before InitPicture
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetColorsList)
        
        #testing colorbyclass. Setting prerequisities:
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.SetColorsList(colorbyclass=True)
        
        #Some values to control
        isochore_label = ['L1', 'L2', 'H1', 'H2', 'H3']
        isochore_values = [37, 41, 46, 53, 100]
        n_of_colors = 5
        
        self.assertEqual(self._test_BaseGrap.isochore_label, isochore_label)
        self.assertEqual(self._test_BaseGrap.isochore_values, isochore_values)
        self.assertEqual(self._test_BaseGrap.n_of_colors, n_of_colors)
        
        #testing color gradient
        self._test_BaseGrap.SetColorsList(colorbyclass=False)
        
        #Some values to control
        isochore_label = ['31.75', '33.5', '35.25', '37.0', '38.75', '40.5', '42.25', '44.0', '45.75', '47.5', '49.25', '51.0', '52.75', '54.5', '56.25', '58.0', '59.75', '61.5', '63.25', '65.0']
        isochore_values = [31.75, 33.5, 35.25, 37.0, 38.75, 40.5, 42.25, 44.0, 45.75, 47.5, 49.25, 51.0, 52.75, 54.5, 56.25, 58.0, 59.75, 61.5, 63.25, 65.0]
        n_of_colors = 20
        
        self.assertEqual(self._test_BaseGrap.isochore_label, isochore_label)
        self.assertEqual(self._test_BaseGrap.isochore_values, isochore_values)
        self.assertEqual(self._test_BaseGrap.n_of_colors, n_of_colors)
        
    def test_GetColorByGClevel(self):
        """Testing color assignment (by Class)"""
        
        GClevels = ((0,3),
                     (31, 3),
                     (37, 3),
                     (37.1, 4),
                     (40, 4),
                     (41, 4),
                     (41.1, 5),
                     (46, 5),
                     (46.1, 6),
                     (53, 6),
                     (53.1, 7),
                     (70, 7),
                     (100, 7)
                    )
        
        #testing function before InitPicture
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.GetColorByGClevel, 30)
        
        #Instantiating the right values
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.SetColorsList(colorbyclass=True)
        
        #testing the function
        for GClevel, value in GClevels:
            color = self._test_BaseGrap.GetColorByGClevel(GClevel)
            self.assertEqual(color, value)
            
        #testing a value returning a warning
        old_errfile = GClib.logger.errfile
        GClib.logger.errfile = StringIO.StringIO()
        
        color = self._test_BaseGrap.GetColorByGClevel(101)
        self.assertEqual(color, 7)
        
        GClib.logger.errfile.seek(0)
        message = GClib.logger.errfile.read()
        
        self.assertTrue("GClevel higher than Maximum value" in message)
        
        GClib.logger.errfile = old_errfile
        
    def test_DrawChName(self):
        """Testing DrawChName"""
        
        #testing function before InitPicture
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawChName, "21")
        
        #instantiating the right values
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.InitPicture()
        
        #Draw Chromosome name
        self._test_BaseGrap.DrawChName("21")
        
    def test_DrawMinMaxValues(self):
        """Testing DrawMinMaxValues"""
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetColorsList)
        
        #testing methods before InitPicture
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SetColorsList)
        
        #Instatiating the class
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.DrawMinMaxValues()
        
    def test_DrawXaxes(self):
        """Testing DrawXaxis"""
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawXaxes)
        
        #testing methods before InitPicture
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawXaxes)
        
        #Call the function in the correct way
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.DrawXaxes()
        
        
    def test_DrawHorizontalLines(self):
        """Testing DrawHorizontalLines"""
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawXaxes)
        
        #testing methods before InitPicture
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawXaxes)
        
        #Call the function in the correct way
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetHorizontalLines([37, 41, 46, 53])
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.DrawHorizontalLines()
        
    def test_Enlargelabels(self):
        """Testing EnlargeLabels"""
        
        #disable log temporarily
        old_threshold = GClib.logger.threshold
        GClib.logger.threshold = 0
        
        #testing methods before SetMaxMinValues
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.DrawXaxes)
        
        #Call the function in the correct way
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetHorizontalLines([37, 41, 46, 53])
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.DrawHorizontalLines()
        self._test_BaseGrap.EnlargeLabels()
        
        #resetting threshold
        GClib.logger.threshold = old_threshold
        
    def test_SaveFigure(self):
        """Testing SaveImage"""
        
        #Call the function in the correct way
        self._test_BaseGrap.SetMinMaxValues(65,30)
        self._test_BaseGrap.SetSequenceLength(1e6)
        self._test_BaseGrap.SetHorizontalLines([37, 41, 46, 53])
        self._test_BaseGrap.InitPicture()
        self._test_BaseGrap.DrawHorizontalLines()
        
        #Get a temporary filename for testing
        testfile = tempfile.mktemp()
        
        #disable log temporarily
        old_threshold = GClib.logger.threshold
        GClib.logger.threshold = 0
        
        #Save the figure
        self._test_BaseGrap.SaveFigure(testfile)
        
        #check no file are overwritten
        self.assertRaises(GClib.Graphs.BaseGraphError, self._test_BaseGrap.SaveFigure, testfile)
        
        #remove temporary file
        if os.path.exists(testfile):
            os.remove(testfile)
            
        #resetting threshold
        GClib.logger.threshold = old_threshold
    
#The testing methods for DrawChromosome classes
class test_DrawChromosome(unittest.TestCase):
    #To test image representation, I need an useful test case. Element.Chromosome
    #is tested by apposited methods
    chromosome = GClib.Elements.Chromosome()
    
    #read isochore from isochores list
    chromosome.LoadIsochores("test_isochores3_chr21.csv")
    isochores = chromosome.isochores
    
    #read windows from windows list
    chromosome.LoadWindows("test_windows_chr21.csv")
    windows = chromosome.windows
    
    #determing sequence length from isochore coordinates
    start = isochores[0].start
    end = isochores[-1].end
    
    sequence_length = end-start
    
    def setUp(self):
        self._test_DrawChromosome = GClib.Graphs.DrawChromosome()
        self._test_DrawChromosome.SetSequenceLength(self.sequence_length)
        self._test_DrawChromosome.InitPicture()
        self._test_DrawChromosome.SetHorizontalLines([37, 41, 46, 53])
        
    def test_DrawIsochoreProfile(self):
        """Testing DrawIsochoreProfile"""
        
        self._test_DrawChromosome.DrawIsochoreProfile(isochores=self.isochores)
        
    def test_DrawWindowProfile(self):
        """Testing DrawIsochoreProfile"""
        
        self._test_DrawChromosome.DrawWindowProfile(windows=self.windows)
        
    def test_DrawIsochoreRectangles(self):
        """Testing DrawIsochoreRectangles"""
        
        self._test_DrawChromosome.SetColorsList(colorbyclass=True)
        self._test_DrawChromosome.DrawIsochoreRectangles(isochores=self.isochores)
        
    def test_DrawWindowRectangles(self):
        """Testing DrawWindowRectangles"""
        
        self._test_DrawChromosome.SetColorsList(colorbyclass=True)
        self._test_DrawChromosome.DrawWindowRectangles(windows=self.windows)
        
    def test_DrawLegend(self):
        """Testing DrawLegend"""
        
        #assert raise when colorbyclass=False
        self._test_DrawChromosome.SetColorsList(colorbyclass=False)
        self._test_DrawChromosome.DrawIsochoreRectangles(isochores=self.isochores)
        self.assertRaises(GClib.Graphs.DrawChromosomeError, self._test_DrawChromosome.DrawLegend)
        
        self._test_DrawChromosome.SetColorsList(colorbyclass=True)
        self._test_DrawChromosome.DrawIsochoreRectangles(isochores=self.isochores)
        self._test_DrawChromosome.DrawLegend()

#TODO: Define test code for drawing graphs

if __name__ == "__main__":
    unittest.main()
