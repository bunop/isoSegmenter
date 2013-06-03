# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:04:39 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

In this file are defined all the Elements needed to deal with isocores. All internal
coordinates are 0 based, with the last coordinate excluded (like python sequences). 
Printed coordinates are 1 based and the last printed coordinate is INCLUDED.

"""

import re
import Bio
import GClib

from numpy.numarray import mlab

#Exceptions for each methods definitions
class ElementError(Exception) : pass
class WindowError(ElementError) : pass
class ChromosomeError(Exception) : pass


class Element:
    """A basic class for windows, gaps and isochores"""
    
    def __init__(self, start=None, end=None):
        self.start = start
        self.end = end
        
        #default values
        self.size = None
        self.Class = None
        
        if start != None or end != None:
            self.SetSize(self.start, self.end)
        
    def __repr__(self):
        try:
            return self.__str__()
        
        except AttributeError:
            raise ElementError, "Element Base Class: you have to define a __str__() method to override this one"
        
    def SetSize(self, start, end):
        """Define the size of the element by specifing the positions"""
        
        if start == None or end == None:
            raise ElementError, "You have to define BOTH start and end, in order to calculate the size of this element"
        
        if start > end:
            start, end = end, start
        
        self.start = start
        self.end = end
        self.size = end - start

class Window(Element):
    """The window class"""
    
    def __init__(self, start=None, end=None, GClevel=None):
        try:
            Element.__init__(self, start, end)
            
        except ElementError, message:
            #if you don't defined BOTH start and end, you will have an error
            raise WindowError, message

        #default value
        self.GClevel = 0
        
        if (GClevel != None):
            self.SetGClevel(GClevel)
        
    def __str__(self):
        return "Window: start:%s,end:%s,size:%s,GClevel:%.3f,Class:%s" %(self.start+1, self.end, self.size, self.GClevel, self.Class)

    def SetGClevel(self, GClevel):
        """To set GClevel and class for a window"""
        self.GClevel = GClevel
        self.Class = CalcClass(GClevel)

#need to define a window class? (which doesn't contain values like average GC, and so on)?

#A class for dealing isochores
class Isochore(Element):
    def __init__(self, window=None):
        """Can instantiate an Isochore from a window"""
        
        self.avg_GClevel = None
        self.stddev_GClevel = None
        self.GClevels = []
        
        if window != None:
            self.GClevels = [window.GClevel]
            self.avg_GClevel = mlab.mean(self.GClevels)
            self.start = window.start
            self.end = window.end
            self.size = window.size
            self.Class = window.Class
    
    def __str__(self):
        return "Isochore: start:%s,end:%s,size:%s,Class:%s,avg_GClevel:%s,stddev_GClevel:%s" %(self.start+1, self.end, self.size, self.Class, self.avg_GClevel, self.stddev_GClevel)
        
    #how many windows define an isochore?
    def __len__(self):
        """return the number of windows in this isochore"""
        return len(self.GClevels)

    def AddWindow(self, window):
        self.GClevels += [window.GClevel]
        self.avg_GClevel = mlab.mean(self.GClevels)
        
        if len(self.GClevels) > 1:
            self.stddev_GClevel = mlab.std(self.GClevels)

        #now I have to find the position of this elements
        if window.start < self.start:
            self.start = window.start
            
        if window.end > self.end:
            self.end = window.end
            
        #TODO: raise an exception while passing a windows containing this element
        self.size = self.end - self.start
        
    def AddIsochore(self, isochore):
        #merging the windows
        self.GClevels += isochore.GClevels
        
        #calculating the isochore parameters
        self.avg_GClevel = mlab.mean(self.GClevels)
        self.stddev_GClevel = mlab.std(self.GClevels)
        
        if isochore.start < self.start:
            self.start = isochore.start
            
        if isochore.end > self.end:
            self.end = isochore.end
            
        #TODO: raise an exception while passing a windows containing this element
        self.size = self.end - self.start
        
        #The class may change:
        self.Class = CalcClass(self.avg_GClevel)
        
    def Test(self,isochore1,isochore2=None):
        """Testing the stddev if this isochores will be added to another one (or two)"""
        
        GClevels = self.GClevels + isochore1.GClevels
        
        if isochore2 != None:
            GClevels += isochore2.GClevels
        
        return mlab.std(GClevels)
        
        
class Gap(Element):
    def __init__(self, start=None, end=None):
        Element.__init__(self, start, end)
        
        #default value
        self.Class = "gap"
        
    def __str__(self):
        return "Gap: start:%s,end:%s,size:%s,class:%s" %(self.start+1, self.end, self.size, self.Class)

#A generic chromosome Class
class Chromosome:
    """A class to deal with chromosomes. This class need a Bio.Seq object to work
    on Gaps, Windows and Isochores. It is better to pass a SeqRecord object at the
    time of Chromosome instantiation, in order to set useful values and scan for
    gaps presence over the sequence."""
    def __init__(self, seqRecord=None):
        self.seqRecord = seqRecord
        self.gaps = []
        self.GClevel = 0
        self.name = None
        self.size = 0
        self.isochores = []
        self.windows = []
        
        if self.seqRecord != None:
            self.gaps = self.Scan4Gaps()
            self.name = seqRecord.name
            self.size = len(seqRecord)
            self.GClevel = self.WholeGCcontent()
            
    def __str__(self):
        if self.seqRecord != None and self.name != self.seqRecord.name:
            self.name = self.seqRecord.name
            
        if self.seqRecord != None and self.size != len(self.seqRecord):
            self.size = len(self.seqRecord)
            
        if self.seqRecord != None and self.GClevel == 0:
            self.GClevel = self.WholeGCcontent()
            
        return "Chromosome: name:%s, size:%s, whole GCcontent:%.3f" %(self.name, self.size, self.GClevel)
    
    def __repr__(self):
        return self.__str__()
    
    def WholeGCcontent(self,seqRecord=None):
        """Calculate the GC counter for the whole chromosome"""
        
        if seqRecord == None:
            seqRecord = self.seqRecord
            
        if type(seqRecord) != Bio.SeqRecord.SeqRecord:
            raise ChromosomeError, "I can calculate the whole GClevel only for a Bio.SeqRecord.SeqRecord object!"
            
        return Bio.SeqUtils.GC(seqRecord.seq)
    
    def Scan4Gaps(self, seqRecord=None):
        """Scan sequence in order to find Gaps"""
        
        if seqRecord == None:
            seqRecord = self.seqRecord
        
        if type(seqRecord) != Bio.SeqRecord.SeqRecord:
            raise ChromosomeError, "I can search for gaps only on Bio.SeqRecord.SeqRecord object!"
        
        #A gaps list
        gaps = []
        
        #A very quick method to find gaps on chromosome
        for match in re.finditer("N+", seqRecord.seq.tostring(), flags=re.IGNORECASE):
            GClib.logger.log(2, "Gap found from %s to %s" %(match.start(), match.end()))
            
            #Adding this gap to the gaps list
            gaps += [Gap(start=match.start(), end=match.end())]
            
        return gaps
        
    def ValueWindows(self, window_size=None, From=None, To=None, gap_tolerance=None):
        """Segments the sequence in non-overlapping windows of fixed size, and calculate
        the GC countent on each window"""
        
        #No window divisions if there is no sequence
        if self.seqRecord == None or type(self.seqRecord) != Bio.SeqRecord.SeqRecord:
            raise ChromosomeError, "The seqRecord class attribute must be instantiated with a valid Bio.Seqrecord object in order to divide sequence in windows"
            
        if self.size == 0 and len(self.seqRecord) == 0:
            raise ChromosomeError, "It makes no sense to divide in windows a sequence of 0 length"
            
        #Setting self.size if it is different from sequence size
        if self.size != len(self.seqRecord):
            self.size = len(self.seqRecord)
            
        if window_size == None:
            GClib.logger.log(3, "No window size provided. Setting the default value to %s bp" %(GClib.WINDOW_SIZE))
            window_size = GClib.WINDOW_SIZE
        
        #Setting gap tolerance if not provided
        if gap_tolerance == None:
            GClib.logger.log(3, "No gap tolerance provided. Setting the default value to %s bp" %(GClib.GAP_TOLERANCE))
            gap_tolerance = GClib.GAP_TOLERANCE
        
        #The user may want to analize sequence between two coordinates. Otherwise I will set the default values
        if From == None:
            From = 0
            
        if To == None:
            To = self.size
            
        #resetting self.windows if any
        self.windows = []
        
        #get a seq object
        seqObj = self.seqRecord.seq
        
        #cicling over the sequence
        start = From
        
        while start < To:
            end = start + window_size
            
            #I cannot take a end coordinate higher than sequence length, or user final coordinate
            if end > To:
                end = To
            
            #cheching gap presence in windows
            for gap in self.gaps:
                if gap.size <= gap_tolerance:
                    #I can ignore this gap
                    continue
                
                if start >= gap.start and start < gap.end:
                    #add this gap to windows list
                    self.windows += [gap]
                    
                    #resizing start and end. The end gap coordinate is excluded (python behavior)
                    new_start = gap.end
                    new_end = new_start + window_size
                    
                    #I cannot take a end coordinate higher than sequence length, or user final coordinate
                    if new_end > To:
                        new_end = To
                        
                    #debug
                    GClib.logger.log(4, "%s found in this window (start:%s,end:%s,size:%s). Setting window to (start:%s,end:%s,size:%s)" %(gap, start, end, end-start, new_start, new_end, new_end-new_start))
                    
                    #Updating the windows coordinates
                    start = new_start
                    end = new_end
                    
                if end >= gap.start and end < gap.end:
                    #resize window in order to avoid the next gap
                    new_end = gap.start
                    
                    #debug
                    GClib.logger.log(4, "%s found in this window (start:%s,end:%s,size:%s). Setting window end to %s (start:%s,end:%s,size:%s)" %(gap, start, end, end-start, new_end, start, new_end, new_end-start))
                    
                    #Updating windows coordinates
                    end = new_end
                    
            #calculate the GClevel of this windows
            GClevel = Bio.SeqUtils.GC(seqObj[start:end])
            
            #Instantiate a new window
            new_window = Window(start=start, end=end, GClevel=GClevel)
            
            #debug
            GClib.logger.log(3, new_window)
            
            #add this windows to the windows list
            self.windows = [new_window]
            
            #update start coordinate for next step
            start = end
            
        
            
        

#A function to define the class of a sequence window
def CalcClass(GClevel):
    """Returns the isochore class of a %GC"""
    
    Class = None
    
    #L1 = 37
    if GClevel <= 37:
        Class = "L1"
    
    #L2 = 41
    elif GClevel > 37 and GClevel <= 41:
        Class = "L2"
        
    #H1 = 46
    elif GClevel > 41 and GClevel <= 46:
        Class = "H1"

    #H2 = 53
    elif GClevel > 46 and GClevel <= 53:
        Class = "H2"
        
    #Otherwise is an H3 class
    else:
        Class = "H3"
        
    return Class
    


