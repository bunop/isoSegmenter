# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:04:39 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

In this file are defined all the Elements needed to deal with isocores. All internal
coordinates are 0 based, with the last coordinate excluded (like python sequences). 
Printed coordinates are 1 based and the last printed coordinate is INCLUDED.

"""

import os
import re
import csv
import sys
import Bio
import numpy
import GClib
import types
import Bio.SeqUtils


from numpy.numarray import mlab

__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

#Exceptions for each methods definitions
class ElementError(Exception) : pass
class WindowError(ElementError) : pass
class IsochoreError(Exception) : pass
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
            
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def SetSize(self, start, end):
        """Define the size of the element by specifing the positions"""
        
        if start == None or end == None:
            raise ElementError, "You have to define BOTH start and end, in order to calculate the size of this element"
        
        if start > end:
            start, end = end, start
        
        #No element of length 0 is admintted
        if start == end:
            raise ElementError, "start and end cannot be equal"
        
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
        """Print Window in 0-Based coordinates"""
        
        return "Window: start:%s,end:%s,size:%s,GClevel:%.6f,Class:%s" %(self.start, self.end, self.size, self.GClevel, self.Class)

    def SetGClevel(self, GClevel):
        """To set GClevel and class for a window"""
        self.GClevel = GClevel
        self.Class = CalcClass(GClevel)

#need to define a window class? (which doesn't contain values like average GC, and so on)?

#A class for dealing isochores
class Isochore():
    def __init__(self, window=None):
        """Can instantiate an Isochore from a window"""
        
        #default values
        self.start = None
        self.end = None
        self.size = None
        self.Class = None
        self.avg_GClevel = None
        self.stddev_GClevel = None
        self.GClevels = []
        
        if window != None:
            self.GClevels = [window.GClevel]
            self.avg_GClevel = round(mlab.mean(self.GClevels),6)
            self.start = window.start
            self.end = window.end
            self.size = window.size
            self.Class = window.Class
    
    def __str__(self):
        """Print Isochore in 0-Based coordinates"""
        
        return "Isochore: start:%s,end:%s,size:%s,Class:%s,avg_GClevel:%s,stddev_GClevel:%s" %(self.start, self.end, self.size, self.Class, self.avg_GClevel, self.stddev_GClevel)
        
    def __repr__(self):
        return self.__str__()
        
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
        
    #how many windows define an isochore?
    def __len__(self):
        """return the number of windows in this isochore"""
        return len(self.GClevels)

    def AddWindow(self, window):
        """Add a window to the current isochore. Windows must to be contiguous and
        non overlapped"""
        
        #Add only a instatiated Window Element
        if window.__class__ != Window or window.Class == None:
            raise IsochoreError, "You can add only an instantiated window to an isochore with AddWindow"
        
        #Raising an exception when adding a non contiguous window. This condition
        #ensures that no windows already inside this isochore will be added. Note that
        #coordinates are in python format
        if window.start != self.end and window.end != self.start:
            raise IsochoreError, "Windows must be contiguous to be added to current isochore"
        
        self.GClevels += [window.GClevel]
        self.avg_GClevel = round(mlab.mean(self.GClevels),6)
        
        if len(self.GClevels) > 1:
            self.stddev_GClevel = round(mlab.std(self.GClevels),6)

        #now I have to find the position of this elements
        if window.start < self.start:
            self.start = window.start
            
        if window.end > self.end:
            self.end = window.end
        
        #determining the new size for this element
        self.size = self.end - self.start
        
        #The class may change when adding a window
        old_Class = self.Class
        self.Class = CalcClass(self.avg_GClevel)
        
        if old_Class != self.Class:
            GClib.logger.err(1, "The Class changed between %s and %s for %s" %(old_Class, self.Class, self))
        
    def AddIsochore(self, isochore):
        """Add and Isochore to the current one. Isochore must be contiguous and 
        non overlapped"""
        
        #Add only a instatiated Isochore Element
        if isochore.__class__ != Isochore or isochore.Class == None:
            raise IsochoreError, "You can add only an instantiated isochore to an isochore with AddIsochore"
            
        #Raising an exception when adding a non contiguous isochore. This condition
        #ensures that no isochores already inside this isochore will be added. Note that
        #coordinates are in python format
        if isochore.start != self.end and isochore.end != self.start:
            raise IsochoreError, "Windows must be contiguous to be added to current isochore"
        
        #merging the windows
        self.GClevels += isochore.GClevels
        
        #calculating the isochore parameters
        self.avg_GClevel = round(mlab.mean(self.GClevels),6)
        self.stddev_GClevel = round(mlab.std(self.GClevels),6)
        
        if isochore.start < self.start:
            self.start = isochore.start
            
        if isochore.end > self.end:
            self.end = isochore.end
            
        #determining the new size for this element
        self.size = self.end - self.start
        
        #The class may change when adding and isochore:
        old_Class = self.Class
        self.Class = CalcClass(self.avg_GClevel)
        
        #Pheraps this event isn't so significant
        if old_Class != self.Class:
            GClib.logger.err(5, "The Class changed between %s and %s for %s" %(old_Class, self.Class, self))
        
        
    def TestHypoSTD(self,isochore1,isochore2=None):
        """Testing the stddev if this isochores is added to another one (or two)"""
        
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
        """Print Gap in 0-Based coordinates"""
        
        return "Gap: start:%s,end:%s,size:%s,Class:%s" %(self.start, self.end, self.size, self.Class)
    

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
            self.Scan4Gaps()
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
        
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False
    
    def WholeGCcontent(self):
        """Calculate the GC counter for the whole chromosome"""
        
        if type(self.seqRecord) != Bio.SeqRecord.SeqRecord:
            raise ChromosomeError, "I can calculate the whole GClevel only for a Bio.SeqRecord.SeqRecord object!"
            
        return Bio.SeqUtils.GC(self.seqRecord.seq)
    
    def Scan4Gaps(self):
        """Scan sequence in order to find Gaps"""
        
        if type(self.seqRecord) != Bio.SeqRecord.SeqRecord:
            raise ChromosomeError, "I can search for gaps only on Bio.SeqRecord.SeqRecord object!"

        GClib.logger.log(2, "Starting GAPs calculation")

        #A gaps list
        gaps = []
        
        #A very quick method to find gaps on chromosome
        for match in re.finditer("N+", self.seqRecord.seq.tostring(), flags=re.IGNORECASE):
            GClib.logger.log(3, "Gap found from %s to %s" %(match.start(), match.end()))
            
            #Adding this gap to the gaps list
            gaps += [Gap(start=match.start(), end=match.end())]
            
        self.gaps = gaps
        
        GClib.logger.log(2, "GAPs calculation finished")
        
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
            GClib.logger.log(1, "No window size provided. Setting the default value to %s bp" %(GClib.WINDOW_SIZE))
            window_size = GClib.WINDOW_SIZE
        
        #Setting gap tolerance if not provided
        if gap_tolerance == None:
            GClib.logger.log(1, "No gap tolerance provided. Setting the default value to %s bp" %(GClib.GAP_TOLERANCE))
            gap_tolerance = GClib.GAP_TOLERANCE
        
        #The user may want to analize sequence between two coordinates. Otherwise I will set the default values
        if From == None:
            From = 0
        
        #Coordinate are 0 based, so
        if To == None:
            To = self.size
        
        GClib.logger.log(2, "Starting window calculation")
        
        #resetting self.windows if any
        self.windows = []
        
        #get a seq object
        seqObj = self.seqRecord.seq
        
        #cicling over the sequence
        start = From
        
        while start < To:
            #determining end coordinates
            end = start + window_size
            
            #I cannot take a end coordinate higher than sequence length, or user final coordinate
            if end > To:
                end = To
            
            #Debug
            GClib.logger.log(5, "Evaluation of window (start:%s, end:%s, size:%s)" %(start, end, end-start))
            
            #cheching gap presence in windows
            for gap in self.gaps:
                if gap.size <= gap_tolerance:
                    #I can ignore this gap
                    continue
                
                #Case 1: Start is inside a GAP
                if start >= gap.start and start < gap.end:
                    #defin a new gap istance with this position. It will be useful when modifing coordinates by user limits
                    new_gap = Gap(start=gap.start, end=gap.end)
                    
                    #Ensure that if the user started from a different position, Gap will start with the same positions
                    if new_gap.start < From:
                        GClib.logger.log(4, "Resizing %s to user start:%s coordinates" %(new_gap, From))
                        new_gap.SetSize(From, new_gap.end) #Start the gap from the same user positions
                    
                    #Maybe the user want to terminate windows calculations before gap ending
                    if new_gap.end > To:
                        GClib.logger.log(4, "Resizing %s to user end:%s coordinates" %(new_gap, To))
                        new_gap.SetSize(new_gap.start, To)
                        
                        #This assignments will interrupt the while cicle
                        start = To
                        end = To
                        
                        #Adding this gap to window list before break the cicle
                        self.windows += [new_gap]
                        
                        #Breaking the cicle
                        break
                    
                    #resizing start and end. The end gap coordinate is excluded (python behavior)
                    #Ensure if the user ended to a determined position, Gap will end with the same position
                    new_start = new_gap.end
                    new_end = new_start + window_size
                    
                    #I cannot take a end coordinate higher than sequence length, or user To coordinates
                    if new_end > To:
                        new_end = To
                        
                    #debug
                    GClib.logger.log(4, "Case 1: %s found in this window (start:%s,end:%s,size:%s). Setting window to (start:%s,end:%s,size:%s)" %(gap, start, end, end-start, new_start, new_end, new_end-new_start))
                    
                    #add this gap to windows list
                    self.windows += [new_gap]
                    
                    #Updating the windows coordinates
                    start = new_start
                    end = new_end
                
                #Case 2: End is inside a GAP
                if end >= gap.start and end < gap.end:
                    #resize window in order to avoid the next gap
                    new_end = gap.start
                    
                    #debug
                    GClib.logger.log(4, "Case 2: %s found in this window (start:%s,end:%s,size:%s). Setting window end to %s (start:%s,end:%s,size:%s)" %(gap, start, end, end-start, new_end, start, new_end, new_end-start))
                    
                    #Updating windows coordinates
                    end = new_end
                    
                #Case 3: A gap inside a window
                if start < gap.start and end >= gap.end:
                    #resize window in order to avoid the next gap
                    new_end = gap.start
                    
                    #debug
                    GClib.logger.log(4, "Case 3: %s found in this window (start:%s,end:%s,size:%s). Setting window end to %s (start:%s,end:%s,size:%s)" %(gap, start, end, end-start, new_end, start, new_end, new_end-start))
                    
                    #Updating windows coordinates
                    end = new_end
                    
                
                
            #In case of sequence terminate with gaps, GAP Case1 will be called, start will be gap.end and end will be equal to To, since end = start+window_size cannot be higher than to.
            #However it makes no sense to calculate a window with 0 length, because sequence is terminated
            if start == end:
                GClib.logger.log(2, "Windows calculation finished")
                break
                
            #calculate the GClevel of this windows
            GClevel = Bio.SeqUtils.GC(seqObj[start:end])
            
            #Round GClevel to first 6 decimal digits
            GClevel = round(GClevel, 6)
            
            #Instantiate a new window
            new_window = Window(start=start, end=end, GClevel=GClevel)
            
            #debug
            GClib.logger.log(3, "New %s defined" %(new_window))
            
            #add this windows to the windows list
            self.windows += [new_window]
            
            #update start coordinate for next step
            start = end
        
    def FindIsochores(self):
        """A function for calculating isochores for this chromosome. Windows must be
        calculated to call this function (call ValueWindows())"""
        
        if self.windows == []:
            raise ChromosomeError, "Windows must be calculated to call this function"
            
        #The isochore calculation is performed in three step. In the first one, we put
        #in the same isochore all adjacent windows with the same class. Resetting isochores if any
        self.isochores = []
        
        #processing all windows
        GClib.logger.log(2, "Starting isochores calculation")
        GClib.logger.log(4, "Isochores calculations. Step (1)...")
        
        for window in self.windows:
            GClib.logger.log(5, "Current %s" %(window))
            if window.Class == "gap":
                #in this Case, i will add a new element
                self.isochores += [window]
                GClib.logger.log(4, "%s added to isochore list" %(window))
            
            #Add a window to an isochore if I have already seen a window
            elif len(self.isochores) > 0 and window.Class == self.isochores[-1].Class:
                GClib.logger.log(4, "Adding %s to %s" %(window, self.isochores[-1]))
                self.isochores[-1].AddWindow(window)
                GClib.logger.log(4, "%s updated" %(self.isochores[-1]))
                
            else:
                self.isochores += [Isochore(window=window)]
                GClib.logger.log(4, "New %s defined" %(self.isochores[-1]))
        
        #Merge isocore below a certain limit (1 window)
        #TODO: define a parameter for a minimun size of an isochore
        GClib.logger.log(4, "Step (1) completed. Starting Step (2)...")
        
        #Now we can merge isochore under a certain size. Since we are modifing the
        #isochore list by removing indexes, it is safer to procede by reversing array. To obtain the -2 value,
        #I have to consider len(isochores) -2 because the length of an array is +1
        #higher than the last index
        for i in range(len(self.isochores)-2,1,-1):
            if self.isochores[i].Class == "gap" or len(self.isochores[i]) > 1:
                GClib.logger.log(5, "Ignoring %s" %(self.isochores[i]))
                continue
            
            else:
                #debug
                GClib.logger.log(5, "Cicle %s. Considering %s:%s, %s:%s and %s:%s" %(i, i, self.isochores[i], i-1, self.isochores[i-1], i+1, self.isochores[i+1]))
                
                #Three test in order to evaluate the most reliable isochore. 
                T1 = None
                T2 = None
                T3 = None
                
                #T1: adding this isochore to the next. Mind the gaps
                if self.isochores[i+1].Class != "gap":
                    T1 = self.isochores[i].TestHypoSTD(self.isochores[i+1])
                
                #T2: adding this isochore to the previous one
                if self.isochores[i-1].Class != "gap":
                    T2 = self.isochores[i].TestHypoSTD(self.isochores[i-1])
                
                #T3: adding this isochore to the previous and the next one
                if self.isochores[i+1].Class != "gap" and self.isochores[i-1].Class != "gap":
                    T3 = self.isochores[i].TestHypoSTD(self.isochores[i+1], self.isochores[i-1])
                    
                #Sorting the three test
                T = dict(T1=T1,T2=T2,T3=T3)
                T = sorted(T.items(), key=lambda x: x[1])
                
                #debug
                GClib.logger.log(5, "Sorted test STDDEV %s" %(T))
                
                #getting the first element != None
                for case, value in T:
                    if value != None:
                        if case == "T1":
                            GClib.logger.log(4, "%s hypothesis selected. Adding %s to %s" %(case, i,i+1))
                            self.isochores[i].AddIsochore(self.isochores[i+1])
                            
                            #deleting i+1
                            del(self.isochores[i+1])
                            
                        elif case == "T2":
                            GClib.logger.log(4, "%s hypothesis selected. Adding isochore %s to %s" %(case, i-1,i))
                            self.isochores[i-1].AddIsochore(self.isochores[i])
                            
                            #deleting i
                            del(self.isochores[i])
                            
                        else:
                            #case T3
                            GClib.logger.log(4, "%s hypothesis selected. Adding isochore %s to %s and %s" %(case, i-1,i,i+1))
                            self.isochores[i-1].AddIsochore(self.isochores[i])
                            self.isochores[i-1].AddIsochore(self.isochores[i+1])
                            
                            #deleting i+1
                            del(self.isochores[i+1])
                            
                            #deleting i
                            del(self.isochores[i])
                        
                        #Last
                        break
                
                #Finding the better case
            
            #debug
            #break
        
        GClib.logger.log(4, "Step (2) completed. Starting Step (3)...")
        
        #Now we could two distinct isochore with the same class, and we want to merge them
        for i in range(len(self.isochores)-2,0,-1):
            if self.isochores[i].Class == self.isochores[i+1].Class:
                #debug
                GClib.logger.log(4, "Merging %s to %s" %(self.isochores[i],self.isochores[i+1]))
                
                #cathing the old class
                old_Class = self.isochores[i].Class
                
                #merge the isochores
                self.isochores[i].AddIsochore(self.isochores[i+1])
                
                #deleting i+1
                del(self.isochores[i+1])
                
                #debug
                GClib.logger.log(4, "%s updated" %(self.isochores[i]))
                
                #May the class change?
                if self.isochores[i].Class != old_Class:
                    GClib.logger.err(1, "The class has changed for %s" %(self.isochores[i]))
                    #at this moment, I want to see this event
                    raise ChromosomeError, "The class has changed for %s" %(self.isochores[i])
            
            #cicle i
        
        GClib.logger.log(4, "Step (3) completed.")
        
        #Print out each isochore instantiation, like windows and gaps
        for isochore in self.isochores:
            GClib.logger.log(3, "New %s defined" %(isochore))
        
        GClib.logger.log(2, "Isochores calculation finished")
    
    def _handle_output(self,outfile):
        """This function open a file for writing if necessary"""
        
        #A flag to determine if I have to close the file (don't close stdout)
        flag_close = False
        
        if type(outfile) == types.StringType:
            #testing for file existance
            if os.path.exists(outfile):
                raise ChromosomeError, "File %s exists. I cannot overwrite it"
                
            #else
            outfile = open(outfile, "w")
            
            #I have to close this file once I've finished
            flag_close = True
            
        elif type(outfile) != types.FileType:
            raise ChromosomeError, "I don't know ho to handle %s : %s" %(outfile, type(outfile))
            
        return outfile, flag_close
        
    def _handle_input(self,infile):
        """This function open a file for reading"""
        
        #A flag to determine if I have to close the file (don't close stdout)
        flag_close = False
        
        if type(infile) == types.StringType:
            #open file
            infile = open(infile, "rU")
            
            #I have to close this file once I've finished
            flag_close = True
            
        elif type(infile) != types.FileType:
            raise ChromosomeError, "I don't know ho to handle %s : %s" %(infile, type(infile))
            
        return infile, flag_close
    
    def DumpGaps(self, outfile=sys.stdout):
        """Dumps gaps in CSV. The output could be an open file handle or
        a filename to write on. Coordinates are 1 based"""
        
        if self.gaps == []:
            raise ChromosomeError, "Gaps must be calculated to call this function"
            
        #Assuming to work with a open filehandle
        outfile, flag_close = self._handle_output(outfile)
        
        #Here, I must have an open file type
        csv_writer = csv.writer(outfile, lineterminator="\n")
        csv_writer.writerow(["Start", "End", "Size"])
        
        for gap in self.gaps:
            csv_writer.writerow([gap.start+1, gap.end, gap.size])
            outfile.flush()
            
        #closing file if necessary
        if flag_close == True: outfile.close()
    
    def LoadGaps(self, infile):
        """Load Gaps from file into chromosome istance. You can pass also a filename
        or an open file hanlde"""
        
        #verify to work with an open file handle
        infile, flag_close = self._handle_input(infile)
        
        #discarding the header
        csv_reader = csv.reader(infile)
        csv_reader.next()
        
        #resetting gaps
        self.gaps = []
        
        for line in csv_reader:
            #beware to integer values
            start, end, size = [int(col) for col in line]
            
            #reset coordinates in python internal coordinates (start 0 based, end excluded)
            self.gaps += [Gap(start=start-1,end=end)]
        
        #closing file if necessary
        if flag_close == True: infile.close()
    
    def DumpWindows(self, outfile=sys.stdout):
        """Dumps windows data in CSV. The output could be an open file handle or
        a filename to write on. Coordinates are 1 based"""
        
        if self.windows == []:
            raise ChromosomeError, "Windows must be calculated with ValueWindows to call this function"
        
        #Assuming to work with a open filehandle
        outfile, flag_close = self._handle_output(outfile)
        
        #Here, I must have an open file type
        csv_writer = csv.writer(outfile, lineterminator="\n")
        csv_writer.writerow(["Start", "End", "Size", "Class", "GClevel"])
        
        for window in self.windows:
            #mind the gap element
            if window.Class == "gap":
                csv_writer.writerow([window.start+1, window.end, window.size, window.Class, None])
                
            else:
                csv_writer.writerow([window.start+1, window.end, window.size, window.Class, "%.6f" %(window.GClevel)])
                
            outfile.flush()
            
        #closing file if necessary
        if flag_close == True: outfile.close()  
        
    def LoadWindows(self, infile):
        """Load windows from file into chromosome istance. Filename or open file 
        handle are accepted"""
        
        #verify to work with an open file handle
        infile, flag_close = self._handle_input(infile)
        
        #discarding the header
        csv_reader = csv.reader(infile)
        csv_reader.next()
        
        #resetting windows
        self.windows = []
        
        for line in csv_reader:
            #beware to integer values. Pay attention to GAPs
            start, end, size, Class, GClevel = int(line[0]), int(line[1]), int(line[2]), line[3], line[4]
            
            #reset coordinates in python internal coordinates (start 0 based, end excluded)
            if Class == "gap":
                self.windows += [Gap(start=start-1,end=end)]
            
            else:
                #GClevel must be float
                GClevel = round(float(GClevel),6)
                
                self.windows += [Window(start=start-1,end=end, GClevel=GClevel)]
                
                #Forcing class assegnation read from file
                self.windows[-1].Class = Class
        
        #closing file if necessary
        if flag_close == True: infile.close()
    
    def DumpIsochores(self, outfile=sys.stdout):
        """Dumps isochores data in CSV. The output could be an open file handle or
        a filename to write on"""
        
        if self.isochores == []:
            raise ChromosomeError, "Isochores must be calculated to call this function"
        
        #Assuming to work with a open filehandle
        outfile, flag_close = self._handle_output(outfile)
        
        #Here, I must have an open file type
        csv_writer = csv.writer(outfile, lineterminator="\n")
        csv_writer.writerow(["Start", "End", "Size", "Class", "AVG_GClevel", "STDDEV_GClevel"])
        
        for isochore in self.isochores:
            #mind the gap element. Coordinates are 1-based
            if isochore.Class == "gap":
                csv_writer.writerow([isochore.start+1, isochore.end, isochore.size, isochore.Class, None, None])
                
            else:
                #If isochore is composed by one element, its stddev will be None
                if len(isochore) > 1:
                    csv_writer.writerow([isochore.start+1, isochore.end, isochore.size, isochore.Class, "%.6f" %(isochore.avg_GClevel), "%.6f" %(isochore.stddev_GClevel)])
                    
                else:
                    csv_writer.writerow([isochore.start+1, isochore.end, isochore.size, isochore.Class, "%.6f" %(isochore.avg_GClevel), None])
                
            outfile.flush()
            
        #closing file if necessary
        if flag_close == True: outfile.close()
        
    def LoadIsochores(self, infile):
        """Load isochores from file into chromosome istance. Filename or open file 
        handle are accepted"""
        
        #verify to work with an open file handle
        infile, flag_close = self._handle_input(infile)
        
        #discarding the header
        csv_reader = csv.reader(infile)
        csv_reader.next()
        
        #resetting isochores
        self.isochores = []
        
        for line in csv_reader:
            #beware to integer values. Pay attention to GAPs
            start, end, size, Class, avg_GClevel, stddev_GClevel = int(line[0]), int(line[1]), int(line[2]), line[3], line[4], line[5]
            
            #reset coordinates in python internal coordinates (start 0 based, end excluded)
            if Class == "gap":
                self.isochores += [Gap(start=start-1,end=end)]
            
            else:
                #avg_GClevel, stddev_GClevel  must be float
                avg_GClevel = round(float(avg_GClevel),6)
                
                #Stddev will be none for isochore of length 0. Here I don't know isochore length (in windows)
                if stddev_GClevel != "": stddev_GClevel = round(float(stddev_GClevel),6)
                
                #Instantiating a new isochore
                isochore = Isochore()
                isochore.start = start - 1 #beware to the 1 based coordinate
                isochore.end = end
                isochore.size = size
                isochore.Class = Class
                isochore.avg_GClevel = avg_GClevel
                isochore.stddev_GClevel = stddev_GClevel
                
                #adding this isochore to isochore list
                self.isochores += [isochore]
        
        #closing file if necessary
        if flag_close == True: infile.close()
    

#A function to define the class of a sequence window
def CalcClass(GClevel):
    """Returns the isochore class of a %GC"""
    
    #Calculating class starting from integer or float
    if type(GClevel) not in [types.IntType, types.FloatType, numpy.float64] :
        raise Exception, "GClevel must be integer or float %s:%s" %(GClevel, type(GClevel))
    
    Class = None
    
    #L1 = 37
    if GClevel <= GClib.CLASS_TO_LEVEL["L1"]:
        Class = "L1"
    
    #L2 = 41
    elif GClevel > GClib.CLASS_TO_LEVEL["L1"] and GClevel <= GClib.CLASS_TO_LEVEL["L2"]:
        Class = "L2"
        
    #H1 = 46
    elif GClevel > GClib.CLASS_TO_LEVEL["L2"] and GClevel <= GClib.CLASS_TO_LEVEL["H1"]:
        Class = "H1"

    #H2 = 53
    elif GClevel > GClib.CLASS_TO_LEVEL["H1"] and GClevel <= GClib.CLASS_TO_LEVEL["H2"]:
        Class = "H2"
        
    #Otherwise is an H3 class
    else:
        Class = "H3"
        
    return Class
    


