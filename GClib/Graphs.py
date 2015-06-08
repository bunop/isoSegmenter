# -*- coding: utf-8 -*-
"""


    Copyright (C) 2013-2015 ITB - CNR

    This file is part of isochoreFinder.

    isochoreFinder is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    isochoreFinder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with isochoreFinder.  If not, see <http://www.gnu.org/licenses/>.


Created on Tue Jun  4 11:05:34 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

A module to deal with isochore graphs. This module is inspired from draw_chromosome.pl
(2002 Jan Paces, http://genomat.img.cas.cz)

"""

import gd
import os
import GClib
import shutil
import types
import tempfile
import GClib.Elements

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

from matplotlib import pyplot

__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

from . import __copyright__, __license__, __version__

#exception definition
class BaseGraphError(Exception) : pass
class DrawChromosomeError(BaseGraphError) : pass
class DrawFamiliesError(Exception): pass

#The maximum and minumu values of DrawChromosome graph (in percentage)
GRAPH_GC_MAX = 65
GRAPH_GC_MIN = 30

class BaseGraph():
    """A base class to make graps like draw chromosomes"""
    
    def __init__(self, sequence_start=0):
        #default values in points (pixel)
        self.scale = 30000 #17500 #higher values shrink images
        self.border = 90 #the white space on the left and on the right of the figure
        self.top = 70 #the upper space before the X axis
        self.bottom = 45 # the bottom size border  
        self.y = 385 #the height of the graphic area (in which isocore are printed, not image height)
        
        #image heigth = self.y + self.bottom
        
        #Other default values
        self.sequence_start = sequence_start
        self.SetFontSize() #default characters size
        
        #A series of values to be defined before printing isochores
        self.y_min = None #the minumum y value printable
        self.y_max = None #the maximum y value printable
        self.py = None #scale y value to be printed in the graph
        self.sequence_length = None #sequence length (affect the x size)
        self.x = None #the width of the image
        self.n_of_h_lines = None #the number of horizontal lines to draw (horizontal grid)
        self.h_lines = None #the values in which horizontal line will be drawn
        self.isochore_label = None #isochores names used to draw legend
        self.isochore_values = None #isocores values used to print horizontal lines and their values
        self.colorslist = None #a list of 
        self.n_of_colors = None #the number of colors
        self.colored_by_class = False #Set True or False if SetColorsList is called with colorbyclass=True or not
        
        #To control image labels with PIL
        self.drawn_labels = False #if labels are drawn by DrawXaxes, it will be True
        self.imagefile = None #if labels are drawn by PIL, this will be the path of image file
        
        #Defined by InitPicture method
        self.graph = None #GD graph instance will be put here
        self.white = None
        self.black = None
        self.gray = None
        
    def __str__(self):
        """A method useful for debugging"""
        
        myclass = str(self.__class__)
        myattributes = self.__dict__
        
        #the returned string
        message = "\n %s Object\n\n" %(myclass)
        
        for key, value in myattributes.iteritems():
            message += "\t%s -> %s\n" %(key, value)
            
        return message
        
    def __repr__(self):
        return self.__str__()
        
    def __del__(self):
        """Delete temporary file if exists"""
        
        if self.imagefile != None:
            if os.path.exists(self.imagefile):
                os.remove(self.imagefile)
        
    def SetMinMaxValues(self, min_value, max_value):
        """Set the maximum and minimum values printable in the graphs"""
        
        if min_value > max_value:
            min_value, max_value = max_value, min_value
        
        self.y_min = min_value
        self.y_max = max_value
        
        #All y values have to be scaled by this value
        self.py = int(round(float(self.y - self.top) / (max_value-min_value))) 
    
    #TODO: define a function in order to set image width in pixels
    def SetSequenceLength(self, sequence_length):
        """Set sequence length and image width"""
        
        if type(sequence_length) != types.IntType:
            sequence_length = int(sequence_length)
        
        self.sequence_length = sequence_length - self.sequence_start
        
        #Now we can calculate the X extension (width)
        self.x = 2 * self.border + int(self.sequence_length / self.scale)
        
        #debug
        GClib.logger.log(3, "Sequence %s bp; image width %s pixels" %(self.sequence_length, self.x))

    def SetFontSize(self, fontsize=gd.gdFontGiant):    
        """Define the labels characters size"""
        
        #gd.gdFontGiant == 4. Higher values are not supported
        if fontsize > 4:
            raise BaseGraphError, "Characters size %s not supported by GD" %(fontsize)
        
        self.fontsize = fontsize
        
    def InitPicture(self):
        """Initialize the figure and make a GD object"""
        
        #Sequence length must be defined to determine the width of the image
        if self.x == None:
            raise BaseGraphError, "Sequence length must be defined by SetSequenceLength"
        
        GClib.logger.log(2, "Starting figure...")
        
        #A gd image. Note that image height is equal to graphic height (self.y) plus
        #upper space (self.top) in which put labels and bottom space (self.bottom) for aesthetic
        self.graph = gd.image((self.x, self.y + self.bottom))
        
        #Allocate base colors
        self.white = self.graph.colorAllocate((255, 255, 255))
        self.black = self.graph.colorAllocate((0, 0, 0))
        self.gray = self.graph.colorAllocate((230, 230, 230))
        
    def SetHorizontalLines(self, arg):
        """Defines the number of horizontal lines. It is possible to specify both
        the number of the horizontal line or a list with the values in which to draw 
        lines"""
        
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
        
        #deal with integer and list in a different way
        if type(arg) == types.IntType:
            #here is the number of lines value
            self.n_of_h_lines = arg
            
            #determining the distance between lines
            step = float(self.y_max - self.y_min) / (self.n_of_h_lines + 1)
            value = self.y_min + step
            
            #Here is the list in which we want to put lines
            self.h_lines = []
            
            while value < self.y_max:
                self.h_lines += [value]
                value += step
        
        elif type(arg) == types.ListType:            
            #The number of horizontal lines to draw is equal to array length
            self.n_of_h_lines = len(arg)
            
            #and the array contains the values in which lines will be drawn
            self.h_lines = arg
            
        else:
            raise BaseGraphError, "Value %s not supported/recognized"

    def SetColorsList(self, colorbyclass=True):
        """Set the possible colors that can be used in an image. The colorbyclass
        values set a distinct color for each class"""
        
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors
            raise BaseGraphError, "InitPicture must be called before this method"

        if colorbyclass == True:
            #setting the flag value in class attribute
            self.colored_by_class = True
            
            #One color for each class. Getting the possible values sorted by GClevel
            #note that GClib.CLASS_TO_LEVEL is a dictionary, so keys could be in random order
            items = sorted(GClib.CLASS_TO_LEVEL.items(), key=lambda x: x[1])
            
            #items are somethin like this: [('L1', 37), ('L2', 41), ('H1', 46), ('H2', 53), ('H3', 100)]
            self.isochore_label = [element[0] for element in items] #isochores names used to draw legend
            self.isochore_values = [element[1] for element in items] #isocores values used to print horizontal lines and their values
            
            #setting the color list
            mycolorslist = []
            
            #Each self.graph.colorAllocate call assign one of 256 true 
            #colors in PNG image and return an ordinal integer of the
            #color already instantiated
            mycolorslist += [self.graph.colorAllocate((0, 100, 255))]
            mycolorslist += [self.graph.colorAllocate((0, 200, 255))]
            mycolorslist += [self.graph.colorAllocate((255, 255, 0))]
            mycolorslist += [self.graph.colorAllocate((255, 130, 0))]
            mycolorslist += [self.graph.colorAllocate((255, 0, 0))]
            
        else:
            #try to define a continue color palette
            mycolorslist = [None for i in range(5*4)]
            
            for i in range(5):
                #base color 0, 0->200, 255
                mycolorslist[i] = self.graph.colorAllocate((0, 50*i, 255))
                
                #base color 0->200, 255, 255
                mycolorslist[i+5] = self.graph.colorAllocate((50*i, 255, 255))
                
                #base color 255, 255, 200->0
                mycolorslist[i+10] = self.graph.colorAllocate((255, 255, 200-50*i))
                
                #base color 255, 200->0, 0 
                mycolorslist[i+15] = self.graph.colorAllocate((255, 200-50*i, 0))
            
            #Now I will define isochore label and values, starting from min and max values
            interval = self.y_max - self.y_min
            step = float(interval) / len(mycolorslist) #20
            self.isochore_values = [self.y_min+step*i for i in range(1, 21)]
            self.isochore_label = [str(value) for value in self.isochore_values]
        
        #memorizzo la lista dei colori
        self.colorslist = mycolorslist
        self.n_of_colors = len(mycolorslist)
        
    def GetColorByGClevel(self, GClevel):
        """Starting from a GClevel values, returns a GD color ID"""
        
        if self.colorslist == None:
            raise BaseGraphError, "SetColorsList must be called before retrive color by GClevel"
        
        #It is possible that we have defined a Maximum value and that I could have and
        #isochore higher than this value. In this case, this element will be plotted outside
        #the box, and I will trow a warning message. In such case the color is the
        #color assigned to the higher value
        color = self.colorslist[-1]
        flag_assigned = False
        
        for i in range(self.n_of_colors-1, -1, -1):
            if GClevel <= self.isochore_values[i]:
                color = self.colorslist[i]
                flag_assigned = True
                
        #Controlling the color assignment
        if flag_assigned == False:
            GClib.logger.err(0, "GClevel higher than Maximum value (%s > %s). Maybe the maximum value have to be raised with SetMinMaxValues" %(GClevel, self.isochore_values[-1]))
        
        #return a GD color ID
        return color
        
    def GetLabelByGClevel(self, GClevel):
        """Starting from a GClevel values, returns the class using GClib.Elements.CalcClass"""
        
        if self.colorslist == None:
            raise BaseGraphError, "SetColorsList must be called before retrive class by GClevel"
            
        return GClib.Elements.CalcClass(GClevel)

    def DrawChName(self,chname):
        """Draws chromosome name inside the graph"""
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors and label
            raise BaseGraphError, "InitPicture must be called before this method"
            
        #chname must be string
        if type(chname) != types.StringType:
            try:
                chname = str(chname)
            except:
                raise BaseGraphError, "Chromosome name must be a string, or something converible in string"
        
        #believe in these values. I could express these values in points, but when you
        #will change self.scale, all these values have to be changed. So express all the
        #coordinates relying on self parameters
        [x1,y1] = [self.border / 6 * 3, int(self.top/6*3)] #the center of the label
        self.graph.arc((x1,y1), (55,40), 0, 360, self.black)
        self.graph.fill((x1,y1),self.black)
        self.graph.string(gd.gdFontGiant,(x1-len(chname)*4,y1-8),chname,self.white);

    def DrawXaxes(self, drawlabels=False):
        """Draw X axis and graduated scale"""
        
        #This function has prerequisites
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
            
        if self.graph == None:
            raise BaseGraphError, "InitPicture must be called before this method"
        
        tick = 100000 * 5
        bigtick = tick * 2
        label = bigtick
        
        #y1 is the y coordinate (in the final image) in which pass the ruled line
        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py))
        
        #this affects temporarily the thickness of all the line. Later in the code I will
        #reset this values to the original values (which I think to be 1)
        self.graph.setThickness(2)
        
        #the ruler line
        self.graph.line((self.border, y1), (self.x-self.border, y1),self.black)
        
        #A big notch on the ruler every "bigtick" bp. If I will start from a sequence position
        #different from 0, I have to put the tick on absolute positions
        for i in range(0, self.sequence_length + self.sequence_start, bigtick):
            if i < self.sequence_start:
                continue
                
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-14),(position,y1),self.black)
            
    
        #a line on the bottom of the graph
        self.graph.line((self.border, int(self.y +1)), (self.x-self.border, int(self.y +1)), self.black)
        
        #reset the thickness of all the line
        self.graph.setThickness(1)
        
        #Now put a small notch on the ruler every "tick" bp
        for i in range(0, self.sequence_length + self.sequence_start, tick):
            if i < self.sequence_start:
                continue
                
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-7),(position,y1),self.black)
        
        #Pheraps it's better to add labels via Python Image Library, because we can enlarge character dimension
        if drawlabels == True:
            #Setting the proper flag
            self.drawn_labels = True
            
            for i in range(0, self.sequence_length + self.sequence_start, label * 2):
                if i < self.sequence_start:
                    continue
                
                position = int((i - self.sequence_start) / self.scale + self.border)
                self.graph.string(self.fontsize, (position-6,y1-30), str(i/label),self.black)
        
            #Write MB at bottom of the ruler
            position = self.x - self.border/5*4
            self.graph.string(self.fontsize,(position,y1-30),"Mb",self.black)
        
    def DrawHorizontalLines(self, drawlabels=True):
        """Draw Horyzontal lines and their value on the left of the graph"""
        
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors
            raise BaseGraphError, "InitPicture must be called before this method"
        
        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py))
        
        y_max, y_min = None, None
        
        #pay attention to self.y_max and self.y_min
        if type(self.y_max) == types.FloatType or type(self.y_min) == types.FloatType:
            y_max = "%.3f" %(self.y_max) + "%"
            y_min = "%.3f" %(self.y_min) + "%"
            
        else:
            y_max = str(self.y_max)+"%"
            y_min = str(self.y_min)+"%"
            
        if drawlabels == True:
            self.graph.string(self.fontsize,(int(self.border / 3), y1-8), y_max, self.black)
            self.graph.string(self.fontsize,(int(self.border / 3),self.y-7), y_min, self.black)
        
        #this is the line style
        self.graph.setStyle((self.black, gd.gdTransparent))
        
        #Draw percentage on the rigth side and horizontal lines
        for i in range(self.n_of_h_lines):
            #don't write a line outside y_max and y_min
            if self.h_lines[i] > self.y_max or self.h_lines[i] < self.y_min:
                continue
            
            y1 = int(round(self.y - (self.h_lines[i]-self.y_min) * self.py))
            label = self.h_lines[i]
            
            if type(label) == types.FloatType:
                label = "%.2f" %(label)
            
            else:
                label = str(label)
            
            self.graph.line((self.border,y1), (self.x-self.border,y1), gd.gdStyled)
            
            if drawlabels == True:
                #Write the value on the left and a dotted line
                self.graph.string(self.fontsize, (int(self.border / 3), y1-8), label + "%",self.black)
            
    def FinishPicture(self,drawlabels=True):
        """Call functions for x,y axis and horizontal lines. Drawlabels flag specifies
        if labels are drawn or not"""
        
        self.DrawXaxes(drawlabels=drawlabels)
        self.DrawHorizontalLines(drawlabels=drawlabels)
        
    def EnlargeLabels(self):
        """Enlarge labels in picture"""
        
        #This function has prerequisites
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
        
        if self.drawn_labels == True:
            raise BaseGraphError, "Labels were drawn by DrawXaxes, and cannot be overwritten by this function"
        
        #determining a temp file for image
        imagefile = tempfile.mktemp(suffix=".png")
        
        #Save the image for the first time
        self.SaveFigure(imagefile)
        
        #Setting the proper attribute to file position
        self.imagefile = imagefile
        
        #Open the temporary image in order to modify it
        im = Image.open(imagefile)
        
        #These are fonts used to draw images. Ensure thay you have the file specified in
        #GClib/__init__.py module
        myfont = ImageFont.truetype(GClib.graph_font_type, 30)
        
        #Questo oggetto mi serve per scriverci dentro
        draw = ImageDraw.Draw(im)
        
        #Determining left size point y1
        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py)) - 45
        
        #the interval (in bp) in which labels will be drawn
        label = 1000000
        iteration = 0
        
        #Wrinting labels. Pay attention to sequence start, no label before sequence starts. 
        for i in range(0, self.sequence_length + self.sequence_start, label):
            if i < self.sequence_start:
                continue
            
            #Write a label in correspondance to thicks, every two "label" distance
            if iteration % 2 == 0:
                position = int(round((i - self.sequence_start) / self.scale + self.border))
                
                #a different X position for different label precision (1, 10, 100)
                if i/label < 10:
                    draw.text((position-7,y1), str(i/label), font=myfont, fill=1)
                elif i/label < 100:
                    draw.text((position-15,y1), str(i/label), font=myfont, fill=1)
                else:
                    draw.text((position-23,y1), str(i/label), font=myfont, fill=1)
                    
            #Next step
            iteration += 1
            
        #For the Mb text
        position = self.x - self.border/5*4
        draw.text((position+5,y1), "Mb", font=myfont, fill=1)
        
        #Now write the percentage labels:
        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py))
        
        y_max, y_min = None, None
        
        #pay attention to self.y_max and self.y_min
        if type(self.y_max) == types.FloatType or type(self.y_min) == types.FloatType:
            y_max = "%.3f" %(self.y_max) + "%"
            y_min = "%.3f" %(self.y_min) + "%"
            
        else:
            y_max = str(self.y_max)+"%"
            y_min = str(self.y_min)+"%"
        
        #Draw max and min values with PIL
        draw.text((int(self.border / 3)-8, y1-12), y_max, font=myfont, fill=1)
        draw.text((int(self.border / 3)-8, self.y-12), y_min, font=myfont, fill=1)
         
        #Draw percentage on the rigth side and horizontal lines
        for i in range(self.n_of_h_lines):
            #don't write a line outside y_max and y_min
            if self.h_lines[i] > self.y_max or self.h_lines[i] < self.y_min:
                continue
         
            y1 = int(round(self.y - (self.h_lines[i]-self.y_min) * self.py))
            label = self.h_lines[i]
            
            if type(label) == types.FloatType:
                label = "%.3f" %(label)
            
            else:
                label = str(label)
            
            draw.text((int(self.border / 3)-8, y1-12), label + "%", font=myfont, fill=1)
            
        #save the new figure
        im.save(imagefile)

    def SaveFigure(self, filename):
        """Draw the image in a new file"""
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors
            raise BaseGraphError, "InitPicture must be called before this method"
        
        #checking for file existance
        if os.path.exists(filename):
            raise BaseGraphError, "File %s exists!!!" %(filename)
        
        #Determing if the Image is already drawn in temporary files
        if self.imagefile == None:
            #write a new image
            self.graph.writePng(filename)
            
            GClib.logger.log(1, "Image written in %s" %(filename))
            
        else:
            #move the temporary image in user files
            shutil.move(self.imagefile, filename)
        
            GClib.logger.log(1, "Image moved in %s" %(filename))
        

#The main class which simulates the behaviour of draw_chromsome.pl
class DrawChromosome(BaseGraph):
    """The main class which simulates the behaviour of draw_chromsome.pl"""
    
    def __init__(self, sequence_start=0):
        """Instantiate the class"""
        
        #Instantiate the base methods and the default attribute class
        BaseGraph.__init__(self, sequence_start=sequence_start)
        
        #The y max and min values are decided by graph type. In this case, GClevel values
        #comprised by 30 and 65 are expected
        self.SetMinMaxValues(GRAPH_GC_MIN, GRAPH_GC_MAX)
    
    def DrawGenericProfile(self, elements, attribute, color, myshift):
        """This function draw elements with a line, which heigth is equal to 'attribute'
        defined by the user. This function is called to draw a window or isochore profile"""
        
        #The x1 grap position is critical to be determined. It depends from first element
        #position and user shift value.
        x1 = int(self.border + round((elements[0].start - self.sequence_start + myshift) / self.scale) - 1)
        
        #This is the lower point drawn in image. It will be needed when representing GAPS
        y2 = self.y
        old_y1 = None
        
        #Allocate color for this profile
        color = self.graph.colorAllocate(color)
        
        #changing thickness for profile
        self.graph.setThickness(2)
        
        #cicling through isochore list
        for element in elements:
            #the length of drawn isochore is derived considering image size and user shift
            x2 = int(self.border + round((element.end - self.sequence_start + myshift) / self.scale) - 1)
        
            GClib.logger.log(5, "Considering %s" %(element))
            
            #When reading from 2006 profile I have more data than what I have to draw
            if x2 > self.border:
                #fix x1 if needed
                if  x1 < self.border:
                    x1 = self.border
                    
                #resize x2 if needed
                if x2 > self.x-self.border:
                    x2 = self.x-self.border
                    
                if element.Class == 'gap':
                    #a grey rectangle which height is image height and lenght is gap length (x2-x1)
                    y1 = int(self.y - (self.y_max-self.y_min) * self.py)
                    
                    #draw the rectangle
                    self.graph.filledRectangle((x1,y1), (x2,y2), self.gray)
                    
                    #each horizontal line is merged to previous line by a vertical
                    #line. In case of gap, no line is needed
                    old_y1 = None
                    
                else:
                    #This is the true isochore
                    y1 = int(self.y - (getattr(element, attribute)-self.y_min) * self.py)
                    
                    #the drawn horizontal line
                    self.graph.line((x1,y1), (x2,y1), color)
                    
                    #draw a vertical line if it is needed
                    if old_y1 != None:
                        self.graph.line((x1,y1), (x1,old_y1), color)
                        
                    #update old_y1 to draw the next element vertical line
                    old_y1 = y1
                
            #this is the new starting point for the new element
            x1 = x2 + 1
            
            #condition: breack if the left size is reached
            if x1 > self.x-self.border:
                break
            
        #resetting the original thickness
        self.graph.setThickness(1)
        
    def DrawIsochoreProfile(self, isochores, color=(0, 0, 0), myshift=0):
        """This function draw isochores with a line, which heigth is equal to class 
        avg_GClevel. This representation can be considered like a profile of isochores.
        User havo to define an isochore list and the color of the line. Graphs can
        be shifted by 'myshift' value"""
        
        #call the generic function
        self.DrawGenericProfile(elements=isochores, attribute="avg_GClevel", color=color, myshift=myshift)
        
        GClib.logger.log(2, "Isochores profile drawn")
        
    def DrawWindowProfile(self, windows, color=(0, 0, 0), myshift=0):
        """This function draw windows with a line, which heigth is equal to window
        GClevel. This representation can be considered like a profile of window.
        User havo to define a windows list and the color of the line. Graphs can
        be shifted by 'myshift' value"""
        
        #call the generic function
        self.DrawGenericProfile(elements=windows, attribute="GClevel", color=color, myshift=myshift)
        
        GClib.logger.log(2, "Windows profile drawn")
    
    def DrawGenericRectangles(self, elements, attribute, myshift):
        """Draw a rectangle in correspondence to elements found. This function needs
        a list of elements to represent and the class attribute to determine rectanlges
        hight. Called by DrawIsochoreRectangles and DrawWindowRectangles"""
    
        #Determining the left coordinataes of boxes. The x1 grap position is critical 
        #to be determined. It depends from first element position and user shift value.
        x1 = int(self.border + round((elements[0].start - self.sequence_start + myshift) / self.scale) - 1)
        
        #This is the lower point drawn in image. It will be needed when representing GAPS
        y2 = self.y
        
        for element in elements:
            #the length of drawn element is derived considering image size and user shift.
            #This value is derived starting from left side of the image every time to
            #avoid that errors on position of first element will be added to last isochores
            x2 = int(self.border + round((element.end - self.sequence_start + myshift) / self.scale) - 1)
            
            if element.Class == 'gap':
                #a grey rectangle which height is image height and lenght is gap length (x2-x1)
                y1 = int(self.y - (self.y_max-self.y_min) * self.py)
                
                #draw the rectangle
                self.graph.filledRectangle((x1,y1), (x2,y2), self.gray)
                
            else:
                color = self.GetColorByGClevel(getattr(element, attribute))
                y1 = int(self.y - (getattr(element, attribute)-self.y_min) * self.py)
                
                #Warning when drawing objects outside max and min values. 
                #getattr(element, attribute) returns GClevel of windows or isochores
                #depending on the type of the class
                if getattr(element, attribute) > self.y_max:
                    GClib.logger.err(0, "Maximum graph point reached (%s > %s). Increase picture max value" %(getattr(element, attribute), self.y_max))
                    
                    #don't write a line outside graph
                    y1 = int(self.y - (self.y_max-self.y_min) * self.py)
                    
                if getattr(element, attribute) < self.y_min:
                    GClib.logger.err(0, "Minimum graph point reached (%s < %s). Decrease picture min value" %(getattr(element, attribute), self.y_min))
                    
                    #don't write a line outside graph
                    y1 = self.y
                
                #draw a colored filled rectangle
                self.graph.filledRectangle((x1,y1), (x2,y2), color)
                
            #updating x1
            x1 = x2 + 1
    
    def DrawIsochoreRectangles(self, isochores, myshift=0):
        """This function draw isochores with filled rectangles, which heigth is equal to 
        class avg_GClevel. This representation can be considered similar to draw_chromosome.pl
        User havo to define an isochore list and evetually provide a value 'myshift'
        to shift the representation"""
        
        #call the generic function
        self.DrawGenericRectangles(elements=isochores, attribute="avg_GClevel", myshift=myshift)
        
        GClib.logger.log(2, "Isochores rectangles drawn")
        
    def DrawWindowRectangles(self, windows, myshift=0):
        """This function draw windows with filled rectangles, which heigth is equal to 
        class GClevel. This representation can be considered similar to draw_chromosome.pl
        User havo to define a window list and evetually provide a value 'myshift'
        to shift the representation"""
        
        #call the generic function
        self.DrawGenericRectangles(elements=windows, attribute="GClevel", myshift=myshift)
        
        GClib.logger.log(2, "Windows rectangles drawn")
    
    def DrawLegend(self):
        """Draw a small legend on the right side of the graph (use with DrawIsochoreRectangles 
        or DrawWindowRectangles and SetColorsList(colorbyclass=True)"""
        
        if self.colored_by_class == False:
            raise DrawChromosomeError, "DrawLegend must be called after SetColorsList(colorbyclass=True)"
        
        #draw a colored box. This is the upper box in legend
        y1 = self.y - (self.y_max - self.y_min) * self.py
        
        #don't write anything outside margins. Draw the first correct colour
        self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.GetColorByGClevel(self.y_max))
        
        #for all the other boxes
        indexes = range(self.n_of_colors-1)
        indexes.reverse()
        
        for i in indexes:
            #don't write anything outside margins
            if self.isochore_values[i] > self.y_max or self.isochore_values[i] < self.y_min :
                continue
        
            #determining where I have to draw
            y1 = self.y - (self.isochore_values[i]-self.y_min)*self.py
            self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.colorslist[i])
        
        #draw labels on legend. The upper label:
        self.graph.string(gd.gdFontGiant, (self.x-self.border/5*3, self.y - (self.y_max - self.y_min) * self.py + 5), self.GetLabelByGClevel(self.y_max), self.black)
        
        #all the remaining labels
        for i in indexes:
            #don't write anything outside margins
            if self.isochore_values[i] > self.y_max or self.isochore_values[i] < self.y_min:
                continue
            
            y1 = self.y-(self.isochore_values[i]-self.y_min) * self.py + 5
            self.graph.string(gd.gdFontGiant, (self.x-self.border/5*3, y1), self.isochore_label[i], self.black)
            
    
# End of DrawChromosome class

# A class to draw isochores like isobase (Schmidt and Frishman 2008 style)
class DrawBarChromosome(BaseGraph):
    """The main class which simulates the behaviour of draw_chromsome.pl"""
    
    def __init__(self, sequence_start=0):
        """Instantiate the class"""
        
        #Instantiate the base methods and the default attribute class
        BaseGraph.__init__(self, sequence_start=sequence_start)

        #re modulate bottom border        
        self.bottom = 60
        
        #Shrink image
        self.scale = 40000
        
    #override basegrap set colour list
    def SetColorsList(self, colorbyclass=True):
        """Set the possible colors that can be used in an image. The colorbyclass
        values set a distinct color for each class"""
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors
            raise BaseGraphError, "InitPicture must be called before this method"

        if colorbyclass == True:
            #setting the flag value in class attribute
            self.colored_by_class = True
            
            #One color for each class. Getting the possible values sorted by GClevel
            #note that GClib.CLASS_TO_LEVEL is a dictionary, so keys could be in random order
            items = sorted(GClib.CLASS_TO_LEVEL.items(), key=lambda x: x[1])
            
            #items are somethin like this: [('L1', 37), ('L2', 41), ('H1', 46), ('H2', 53), ('H3', 100)]
            self.isochore_label = [element[0] for element in items] #isochores names used to draw legend
            self.isochore_values = [element[1] for element in items] #isocores values used to print horizontal lines and their values
            
            #setting the color list
            mycolorslist = []
            
            #Each self.graph.colorAllocate call assign one of 256 true 
            #colors in PNG image and return an ordinal integer of the
            #color already instantiated
            mycolorslist += [self.graph.colorAllocate((0, 100, 255))]
            mycolorslist += [self.graph.colorAllocate((0, 200, 255))]
            mycolorslist += [self.graph.colorAllocate((255, 255, 0))]
            mycolorslist += [self.graph.colorAllocate((255, 130, 0))]
            mycolorslist += [self.graph.colorAllocate((255, 0, 0))]
            
        else:
            #try to define a continue color palette
            raise BaseGraphError, "Not yet implemented"
        
        #memorizzo la lista dei colori
        self.colorslist = mycolorslist
        self.n_of_colors = len(mycolorslist)
    
    #Overriding X axes
    def DrawXaxes(self, drawlabels=False):
        """Draw X axis and graduated scale"""
            
        if self.graph == None:
            raise BaseGraphError, "InitPicture must be called before this method"
        
        tick = 100000 * 5
        bigtick = tick * 2
        label = bigtick
        
        #y1 is the y coordinate (in the final image) in which pass the ruled line
        y1 = self.top
        
        #this affects temporarily the thickness of all the line. Later in the code I will
        #reset this values to the original values (which I think to be 1)
        self.graph.setThickness(2)
        
        #the ruler line
        self.graph.line((self.border, y1), (self.x-self.border, y1),self.black)
        
        #A big notch on the ruler every "bigtick" bp. If I will start from a sequence position
        #different from 0, I have to put the tick on absolute positions
        for i in range(0, self.sequence_length + self.sequence_start, bigtick):
            if i < self.sequence_start:
                continue
                
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-14),(position,y1),self.black)
            
    
        #a line on the bottom of the graph
        self.graph.line((self.border, int(self.y +1)), (self.x-self.border, int(self.y +1)), self.black)
        
        #reset the thickness of all the line
        self.graph.setThickness(1)
        
        #Now put a small notch on the ruler every "tick" bp
        for i in range(0, self.sequence_length + self.sequence_start, tick):
            if i < self.sequence_start:
                continue
                
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-7),(position,y1),self.black)
        
        #Pheraps it's better to add labels via Python Image Library, because we can enlarge character dimension
        if drawlabels == True:
            #Setting the proper flag
            self.drawn_labels = True
            
            for i in range(0, self.sequence_length + self.sequence_start, label * 2):
                if i < self.sequence_start:
                    continue
                
                position = int((i - self.sequence_start) / self.scale + self.border)
                self.graph.string(self.fontsize, (position-6,y1-30), str(i/label),self.black)
        
            #Write MB at bottom of the ruler
            position = self.x - self.border/5*4
            self.graph.string(self.fontsize,(position,y1-30),"Mb",self.black)
            
    def FinishPicture(self,drawlabels=True):
        """Call functions for x,y axis and horizontal lines. Drawlabels flag specifies
        if labels are drawn or not"""
        
        self.DrawXaxes(drawlabels=drawlabels)
        #self.DrawHorizontalLines(drawlabels=drawlabels)

        
    def DrawGenericRectangles(self, elements, attribute, myshift):
        """Draw a rectangle in correspondence to elements found. This function needs
        a list of elements to represent and the class attribute to determine rectanlges
        hight. Called by DrawIsochoreRectangles"""
    
        #Determining the left coordinataes of boxes. The x1 grap position is critical 
        #to be determined. It depends from first element position and user shift value.
        x1 = int(self.border + round((elements[0].start - self.sequence_start + myshift) / self.scale) - 1)
        
        #This is the lower point drawn in image
        y2 = self.y
        
        for element in elements:
            #the length of drawn element is derived considering image size and user shift.
            #This value is derived starting from left side of the image every time to
            #avoid that errors on position of first element will be added to last isochores
            x2 = int(self.border + round((element.end - self.sequence_start + myshift) / self.scale) - 1)
            
            #When reading from 2006 profile I have more data than what I have to draw
            if x2 > self.border:
                #fix x1 if needed
                if  x1 < self.border:
                    x1 = self.border
                    
                #resize x2 if needed
                if x2 > self.x-self.border:
                    x2 = self.x-self.border
            
                if element.Class == 'gap':
                    #a grey rectangle which height is image height and lenght is gap length (x2-x1)
                    y1 = self.top
                    
                    #draw the rectangle
                    self.graph.filledRectangle((x1,y1), (x2,y2), self.gray)
                    
                else:
                    #Set the color of isochore
                    color = self.GetColorByGClevel(getattr(element, attribute))
                    
                    #a colored rectangle which height is image height and lenght is isochore length (x2-x1)
                    y1 = self.top
                    
                    #draw a colored filled rectangle
                    self.graph.filledRectangle((x1,y1), (x2,y2), color)
                
            #updating x1
            x1 = x2 + 1
            
            #condition: breack if the left size is reached
            if x1 > self.x-self.border:
                break
            
    def DrawIsochoreRectangles(self, isochores, myshift=0):
        """This function draw isochores with filled rectangles, which heigth is equal to 
        class avg_GClevel. This representation can be considered similar to draw_chromosome.pl
        User havo to define an isochore list and evetually provide a value 'myshift'
        to shift the representation"""
        
        #call the generic function
        self.DrawGenericRectangles(elements=isochores, attribute="avg_GClevel", myshift=myshift)
        
        GClib.logger.log(2, "Isochores bar rectangles drawn")
        
    #Overriding enlarge labels
    def EnlargeLabels(self):
        """Enlarge labels in picture"""
        
        if self.drawn_labels == True:
            raise BaseGraphError, "Labels were drawn by DrawXaxes, and cannot be overwritten by this function"
        
        #determining a temp file for image
        imagefile = tempfile.mktemp(suffix=".png")
        
        #Save the image for the first time
        self.SaveFigure(imagefile)
        
        #Setting the proper attribute to file position
        self.imagefile = imagefile
        
        #Open the temporary image in order to modify it
        im = Image.open(imagefile)
        
        #These are fonts used to draw images. Ensure thay you have the file specified in
        #GClib/__init__.py module
        myfont = ImageFont.truetype(GClib.graph_font_type, 30)
        
        #Questo oggetto mi serve per scriverci dentro
        draw = ImageDraw.Draw(im)
        
        #Determining left size point y1
        y1 = self.top - 45
        
        #the interval (in bp) in which labels will be drawn
        label = 1000000
        iteration = 0
        
        #Wrinting labels. Pay attention to sequence start, no label before sequence starts. 
        for i in range(0, self.sequence_length + self.sequence_start, label):
            if i < self.sequence_start:
                continue
            
            #Write a label in correspondance to thicks, every two "label" distance
            if iteration % 2 == 0:
                position = int(round((i - self.sequence_start) / self.scale + self.border))
                
                #a different X position for different label precision (1, 10, 100)
                if i/label < 10:
                    draw.text((position-7,y1), str(i/label), font=myfont, fill=1)
                elif i/label < 100:
                    draw.text((position-15,y1), str(i/label), font=myfont, fill=1)
                else:
                    draw.text((position-23,y1), str(i/label), font=myfont, fill=1)
                    
            #Next step
            iteration += 1
            
        #For the Mb text
        position = self.x - self.border/5*4
        draw.text((position+5,y1), "Mb", font=myfont, fill=1)
            
        #save the new figure
        im.save(imagefile)

    
class DrawFamilies:
    """A class to plot isochores families in histograms"""
    
    def __init__(self,families=None):
        """Instantiate the class starting from Families Element"""
        
        if families.__class__ != GClib.Elements.Families or families.data == {} :
            raise DrawFamiliesError, "This class must be instantiated only by a valid Families Element Class"""
        
        #setting families element
        self.families = families
        
        #data will be plotted using bar plot
        self.all_bar = []
        
        #Image proportions
        scale = 20.0 / 12

        self.fig = pyplot.figure(figsize=(13,13/scale))
        self.fontsize=30 #"x-large"

        #instantiate a bar graph
        for bin, bin_data in self.families.data.iteritems():
            #print "%s:%s" %(bins[i], data[i]),
            length = bin_data["size"]
            self.all_bar += [pyplot.bar(bin-families.bin_size*0.2, int(round(length/1e6,0)), width=0.4*families.bin_size, bottom=0)]
        
        #setting axes
        self.x_min = int(round(families.min_value / families.precision,0))
        self.x_max = int(round(families.max_value / families.precision,0))

    def __str__(self):
        """A method useful for debugging"""
        
        myclass = str(self.__class__)
        myattributes = self.__dict__
        
        #the returned string
        message = "\n %s Object\n\n" %(myclass)
        
        for key, value in myattributes.iteritems():
            message += "\t%s -> %s\n" %(key, value)
            
        return message
        
    def __repr__(self):
        return self.__str__()

    def DrawAxisLabels(self):
        """Set axis label"""
        
        myticks = []
        mylabels = []
        
        #To evitate labels in axis origin
        flag_ticks = False
        
        #To draw ticks on X axis, set a list for ticks and labels
        for i in range(self.families.n_of_bins):
            if i % self.families.precision == 0:
                #A tick on each bin
                myticks += [self.families.bins[i]]
                
                #a label every 2 bins
                if i % (self.families.precision*2) == 0:
                    if flag_ticks == False:
                        flag_ticks = True
                        mylabels += ['']
                    
                    else:
                        mylabels += [int(self.families.bins[i])] 
                    
                else:
                    mylabels += ['']
        
        #debug
        #print myticks, len(myticks)
        #print mylabels, len(mylabels)
        
        #draw the thicks and labels
        pyplot.xticks([myticks[i] for i in range(0,len(myticks))],[mylabels[i] for i in range(0,len(myticks))],size=self.fontsize)
        pyplot.yticks(size=self.fontsize)

    def SetAxisLimits(self, axis=[0,0,0,0]):
        """Setting axis [Xmin, Xmax, Ymin, Ymax]"""
    
        #recording the internal axis values
        tmp_axis = list(pyplot.axis())
    
        #function call without parameters. Setting axis relying on x_max, x_min if both are equal to 0
        if axis[0] == 0 and axis[1] == 0:
            #changing x values
            axis[0] = self.x_min
            axis[1] = self.x_max
        
        #setting Y value if they are not defined by user   
        if axis[2] == 0 and axis[3] == 0:
            axis[2] = tmp_axis[2]
            axis[3] = tmp_axis[3]
        
        #shift the axis down
        #axis[2] -= int(round(axis[3]*(10.0/300),0))
        
        #Applying the axis
        pyplot.axis(axis)

    def DrawGrid(self):
        """Draw grid in grapsh"""
        
        #this draws the grid in the graph
        self.axes = pyplot.axes()
        self.axes.grid(linewidth=2)

    def DrawTitle(self, title):
        """Draws title in graph"""
        
        pyplot.text(0.8, 0.8, title,
            horizontalalignment='center',
            verticalalignment='center',
            transform = pyplot.gca().transAxes,
            size=40
        )
        
    def SaveFigure(self, filename, dpi=100):
        """Draw the image in a new file. DPI quality can be specified"""
        
        #checking for file existance
        if os.path.exists(filename):
            raise DrawFamiliesError, "File %s exists!!!" %(filename)
            
        #save picture in file
        pyplot.savefig(filename, dpi=dpi)
        
        GClib.logger.log(1, "Image written in %s" %(filename))

#debug: define a test function to works on BaseGraph
def test_BaseGraph(filename="test.png"):
    graph = BaseGraph()
    graph.SetMinMaxValues(30,65)
    graph.SetSequenceLength(5e7)
    graph.InitPicture()
    #graph.SetHorizontalLines([37, 41, 46, 53])
    graph.SetHorizontalLines(5)
    graph.SetColorsList(colorbyclass=True)
    graph.DrawChName("21")
    graph.DrawXaxes(drawlabels=True)
    graph.DrawHorizontalLines()
    
    #Draw the image
    graph.SaveFigure(filename)
    
    #return the object for testing
    return graph
    
    
    
