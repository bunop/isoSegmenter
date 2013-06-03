# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:04:39 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

In this file are defined all the Elements needed to deal with isocores

"""

from numpy.numarray import mlab

class Element:
    """A basic class for windows, gaps and isochores"""
    
    def __init__(self):
        self.start = None
        self.end = None
        self.size = None
        self.Class = None
        
    def __repr__(self):
        return self.__str__()
        
    def SetSize(self, start, end):
        """Define the size of the element by specifing the positions"""
        
        if start > end:
            start, end = end, start
        
        self.start = start
        self.end = end
        self.size = end - start

class Window(Element):
    """The window class"""
    
    def __init__(self):
        Element.__init__(self)
        self.GClevel = 0
        
    def __str__(self):
        return "Window: start:%s,end:%s,size:%s,GClevel:%.3f,Class:%s" %(self.start, self.end, self.size, self.GClevel, self.Class)

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
        return "Isochore: start:%s,end:%s,size:%s,Class:%s,avg_GClevel:%s,stddev_GClevel:%s" %(self.start, self.end, self.size, self.Class, self.avg_GClevel, self.stddev_GClevel)
        
    def __repr__(self):
        return self.__str__()
        
    #how many windows define an isochore?
    def __len__(self):
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
    def __init__(self):
        Element.__init__(self)
        self.Class = "gap"
        
    def __str__(self):
        return "Gap: start:%s,end:%s,size:%s,class:%s" %(self.start, self.end, self.size, self.Class)

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
    



