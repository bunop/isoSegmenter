# -*- coding: utf-8 -*-
"""
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
import Image
import ImageDraw
import ImageFont
import tempfile

__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

#exception definition
class BaseGraphError(Exception) : pass
class DrawChromosomeError(BaseGraphError) : pass

class BaseGraph():
    """A base class to make graps like draw chromosomes"""
    
    def __init__(self, sequence_start=0):
        #default values in points (pixel)
        self.scale = 30000 #17500 #higher values shrink images
        self.border = 90 #the white space on the left and on the right of the figure
        self.top = 70 #the upper space before the X axis
        self.y = 385 #the height of the graphic area (in which isocore are printed, not image height)
        
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
        #upper space (self.top) in which put labels and bottom space (self.top/2) for aesthetic
        self.graph = gd.image((self.x, self.y + self.top / 2))
        
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

    def DrawMinMaxValues(self):
        """Draw labels for Min and Max values"""
        
        if self.y_min == None or self.y_max == None:
            raise BaseGraphError, "Max and Min y values must be defined by SetMinMaxValues"
        
        if self.graph == None:
            #if GD image isn't instantiated yed, I couldn't instantiate colors
            raise BaseGraphError, "InitPicture must be called before this method"

        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py))
        
        #pay attention to self.y_max and self.y_min
        if type(self.y_max) == types.FloatType or type(self.y_min) == types.FloatType:
            self.graph.string(self.fontsize,(int(self.border / 4), y1-8), "%.3f" %(self.y_max) + "%", self.black)
            self.graph.string(self.fontsize,(int(self.border / 4),self.y-7), "%.3f" %(self.y_min) + "%",self.black)
            
        else:
            #Si suppone che siano degli interi o stringhe            
            self.graph.string(self.fontsize,(int(self.border / 3),y1-8),str(self.y_max)+"%",self.black)
            self.graph.string(self.fontsize,(int(self.border / 3),self.y-7),str(self.y_min)+"%",self.black)

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
        
        #A big notch on the ruler every "bigtick" bp
        for i in range(self.sequence_start, self.sequence_length + self.sequence_start, bigtick):
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-14),(position,y1),self.black)
            
        #a line on the bottom of the graph
        self.graph.line((self.border, int(self.y +1)), (self.x-self.border, int(self.y +1)), self.black)
        
        #reset the thickness of all the line
        self.graph.setThickness(1)
        
        #Now put a small notch on the ruler every "tick" bp
        for i in range(self.sequence_start, self.sequence_length + self.sequence_start, tick):
            position = int((i - self.sequence_start) / self.scale + self.border)
            self.graph.line((position,y1-7),(position,y1),self.black)
        
        #Pheraps it's better to add labels via Python Image Library, because we can enlarge character dimension
        if drawlabels == True:
            #Setting the proper flag
            self.drawn_labels = True
            
            for i in range(self.sequence_start, self.sequence_length + self.sequence_start, label * 2):
                position = int((i - self.sequence_start) / self.scale + self.border)
                self.graph.string(self.fontsize, (position-6,y1-30), str(i/label),self.black)
        
            #Write MB at bottom of the ruler
            position = self.x - self.border/5*4
            self.graph.string(self.fontsize,(position,y1-30),"Mb",self.black)
        
    def DrawHorizontalLines(self):
        """Draw Horyzontal lines and their value on the left of the graph"""
        
        #per le linee
        self.graph.setStyle((self.black, gd.gdTransparent))
        
        #Sono le percentuali a SX dell'immagine e le loro linee orizzontali (nuova versione)
        for i in range(self.n_of_h_lines):
            y1 = int(round(self.y - (self.h_lines[i]-self.y_min) * self.py))
            label = self.h_lines[i]
            
            if type(label) == types.FloatType:
                label = "%.2f" %(label)
            
            else:
                label = str(label)
            
            #Write the value on the left and a dotted line
            self.graph.string(self.fontsize, (int(self.border / 3), y1-8), label + "%",self.black)
            self.graph.line((self.border,y1), (self.x-self.border,y1), gd.gdStyled)
            
    def FinishPicture(self,drawlabels=True):
        """Call functions for x,y axis and horizontal lines. Drawlabels flag specifies
        if labels are drawn or not"""
        
        self.DrawMinMaxValues()
        self.DrawXaxes(drawlabels=drawlabels)
        self.DrawHorizontalLines()
        
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
        
        #Una volta salvato il grafico, è il momento di tirarsi le storie per la dimensione delle scritte
        im = Image.open(imagefile)
        
        #Carico i font con cui scrivere dentro l'immagine
        myfont = ImageFont.truetype(GClib.graph_font_type, 30)
        
        #Questo oggetto mi serve per scriverci dentro
        draw = ImageDraw.Draw(im)
        
        #Determining left size point y1
        y1 = int(round(self.y - (self.y_max - self.y_min) * self.py)) - 45
        
        #the interval (in bp) in which labels will be drawn
        label = 1000000
        
        #Wrinting labels
        for i in range(self.sequence_start, self.sequence_length + self.sequence_start, label*2):
            position = int(round((i - self.sequence_start) / self.scale + self.border))
            
            #a different X position for different label precision (1, 10, 100)
            if i/label < 10:
                draw.text((position-7,y1), str(i/label), font=myfont, fill=1)
            elif i/label < 100:
                draw.text((position-15,y1), str(i/label), font=myfont, fill=1)
            else:
                draw.text((position-23,y1), str(i/label), font=myfont, fill=1)
            
        #For the Mb text
        position = self.x - self.border/5*4
        draw.text((position+5,y1), "Mb", font=myfont, fill=1)
        
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
        self.SetMinMaxValues(30,65)
    
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
        
        #draw a colored box
        y1 = self.y - (self.y_max - self.y_min) * self.py
        
        #this is the upper box in legend
        self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.colorslist[-1])
        
        #for all the other boxes
        indexes = range(self.n_of_colors-1)
        indexes.reverse()
        
        for i in indexes:
            y1 = self.y - (self.isochore_values[i]-self.y_min)*self.py
            self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.colorslist[i])
        
        #draw labels on legend. The upper label:
        self.graph.string(gd.gdFontGiant, (self.x-self.border/5*3, self.y - (self.isochore_values[-1]-self.y_max) * self.py + 5), self.isochore_label[-1], self.black)
        
        #all the remaining labels
        for i in range(self.n_of_colors-1):
            y1 = self.y-(self.isochore_values[i]-self.y_min) * self.py + 5
            self.graph.string(gd.gdFontGiant, (self.x-self.border/5*3, y1), self.isochore_label[i], self.black)
            
    

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
    graph.DrawMinMaxValues()
    graph.DrawXaxes(drawlabels=True)
    graph.DrawHorizontalLines()
    
    #Draw the image
    graph.SaveFigure(filename)
    
    #return the object for testing
    return graph
    
    
    