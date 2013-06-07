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
import types

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
        
        #Defined by InitPicture method
        self.graph = None #GC graph instance will be put here
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
            

    def SaveFigure(self, filename):
        """Draw the image in a new file"""
        
        #checking for file existance
        if os.path.exists(filename):
            raise BaseGraphError, "File %s exists!!!" %(filename)
        
        self.graph.writePng(filename)
        
        GClib.logger.log(1, "Image written in %s" %(filename))
        

#The main class which simulates the behaviour of draw_chromsome.pl
class DrawChromosome(BaseGraph):
    """The main class which simulates the behaviour of draw_chromsome.pl"""
    
    def __init__(self, sequence_start=0):
        """Instantiate the class"""
        
        #Istanzio le variabili della classe base
        BaseGraph.__init__(self, sequence_start=sequence_start)
        
        #Le scala dell'asse y (dipende dal tipo di grafico)
        self.SetMinMaxValues(30,65)
        
    #da sistemare!
    def DrawColorGraph(self,data):
        """Disegna una linea colorata per ogni valore letto"""
        points = len(data)
        points2 = int(self.seq_len / self.window)
        px = int((self.x - self.border*2) / points2)
        y_base = self.y
        
        for i in range(points):
            x1 = self.border + (px * i)
            y1 = self.y - ((data[i] - self.y_min) * self.py)
            color = self.ResolveColor(data[i])
            self.graph.line((x1, y1), (x1, y_base), color)
    
    #Questa funzione vorrei che segnasse le isocore come una linea continua
    #con altezza pari alla percentuale delle iscore
    def DrawGraph(self, first_gruppo, media_attributo, color=(0, 0, 0), myshift=0):
        """Disegna una linea che rappresenta la percentuale in GC lungo il genoma.
        E' possibile specificare il colore come terna RGB. Inoltre e possibile fissare
        a posteriori uno shift per cui spostare i gruppi"""
        
        #Come per DrawRectangleIsocore, mi serve leggere la prima delle strutture
        gruppo = first_gruppo
        
        #La posizione di inizio è critica, se non ci sono dati
        x1 = int(self.border + round((gruppo.start - self.seq_start + myshift) / self.scale) - 1)
        
        #questo è il punto più basso dell'immagine. Serve solo per disegnare i gap
        y2 = self.y
        old_y1 = None
        
        #devo costruire il colore
        color = self.graph.colorAllocate(color)
        
        #per le linee. Aumento lo spessore
        self.graph.setThickness(2)
        
        #parto con il ciclo
        while gruppo != None:
            #Calcolo la lunghezza delle regione in questione. Anziche sommarè ogni volta
            #lo ricalcolo ogni volta cosi come DrawRectangleIsocore
            x2 = int(self.border + round((gruppo.end - self.seq_start + myshift) / self.scale) - 1)
        
            logger.log(5, "Processo la regione %s : %s" %(gruppo.start, gruppo.end))
            
            if gruppo.classe == 'gap':
                #come per draw_chromosome, il gap crea un rettangolo alto come
                #tutto il grafico
                y1 = int(self.y - (self.y_max-self.y_min) * self.py)
                
                #ok, adesso posso picchiarci dentro un rettangolo
                self.graph.filledRectangle((x1,y1), (x2,y2), self.gray)
                
                #se c'è il gap, non voglio disegnare linee verticali
                old_y1 = None
                
            else:
                #se non sto trattando un Gap, allora devo disegnare una linea
                #print gruppo.classe
                
                y1 = int(self.y - (getattr(gruppo, media_attributo)-self.y_min) * self.py)
                
                #questa è la linea orizzontale
                self.graph.line((x1,y1), (x2,y1), color)
                
                #sarebbe bello disegnare anche le linee verticali
                if old_y1 != None:
                    self.graph.line((x1,y1), (x1,old_y1), color)
                    
                #aggiorno il valore di old_y1
                old_y1 = y1
                
            #aggiorno i valori e passo alla isocora successiva
            x1 = x2 + 1
            
            gruppo = gruppo.GetNextGruppo()
        
        #per le linee. Ripristino lo spessore
        self.graph.setThickness(1)
        
        #debug
        logger.log(2, "Isocore disegnate")
        
        return
        
    def DrawColoredRectangles(self, elements):
        """permette di disegnare dei rettangoli in corrispondenza delle
        isocore individuate. Necessita di una struttura gruppo di partenza e
        dell'attributo da graficare (media del GC del gruppo o TN)"""
        
        #ok posso procedere con le mie analisi. La prima cosa che mi manca
        #Quindi devo partire dalla prima struttura e cominciare a fare i rettangoli
        #in base alla classe decido se fare un rettangolo per un gap, oppure per un isocora
        x1 = self.border
        y2 = self.y
        
        for element in elements:
            #devo determinare dove si posiziona questa regione. La larghezza della
            #regione è valida allo stesso modo per i gap e per le isocore. Start e
            #end devono essere scalati in modo opportuno
            #x2 = x1 + math.ceil(gruppo.size / self.scale)-1
            
            #se aggingo un rettangolo a quello precedente, ogni volta che commetto un errore,
            #lo trascino nella prossima isocora. Allora la coordinata x2 deve essere calcolata
            #solamente in base alla lunghezza dell'isocora
            x2 = int(self.border + round( (element.end - self.seq_start) / self.scale) - 1)
            
            if element.Class == 'gap':
                #come per draw_chromosome, il gap crea un rettangolo alto come
                #tutto il grafico
                y1 = int(self.y - (self.y_max-self.y_min) * self.py)
                
                #ok, adesso posso picchiarci dentro un rettangolo
                self.graph.filledRectangle((x1,y1), (x2,y2), self.gray)
                
            else:
                #in questo caso, la media dell'attributo che voglio misurare mi darà
                #l'altezza del rettangolo, e pure il colore delle storie
                try:
                    color = self.ResolveColor(getattr(element, "avg_GClevel"))
                    y1 = int(self.y - (getattr(element, "avg_GClevel")-self.y_min) * self.py)
                
                except AttributeError:
                    color = self.ResolveColor(getattr(element, "GClevel"))
                    y1 = int(self.y - (getattr(element, "GClevel")-self.y_min) * self.py)
                
                #provo a metterci un rettangolo colorato
                self.graph.filledRectangle((x1,y1), (x2,y2), color)
                
            #aggiorno x1 e passo alla isocora successiva
            x1 = x2 + 1
    
    def DrawLegend(self):
        """Serve per segnare sul grafico la legenda delle isocore."""
        
        #per i box
        y1 = self.y - (self.y_max - self.y_min) * self.py
        self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.colorslist[-1])
        
        lista = range(self.n_of_colors-1)
        lista.reverse()
        
        for i in lista:
            y1 = self.y - (self.isochore_proc[i]-self.y_min)*self.py
            self.graph.filledRectangle((self.x-self.border/5*4, y1), (self.x-self.border/5, self.y), self.colorslist[i])
        
        #per le etichette. Se non ho isochore label è perchè non ho colorato per famiglie
        self.graph.string(gd.gdFontGiant, (self.x-self.border/5*3, self.y - (self.isochore_proc[-1]-self.y_max) * self.py + 5), self.isochore_label[-1], self.black)
        
        for i in range(self.n_of_colors-1):
            y1 = self.y-(self.isochore_proc[i]-self.y_min) * self.py + 5
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
    
    
    