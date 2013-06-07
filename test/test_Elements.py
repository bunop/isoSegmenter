# -*- coding: utf-8 -*-
"""

Created on Fri May 31 17:07:40 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

"""

import sys
import unittest

sys.path.append("..")

import GClib
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
        """Testing element instantiation"""
        
        self._test_Element = GClib.Elements.Element()
        
    def test_SetSize(self):
        """Testing element set size"""
        
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
        
        
        

if __name__ == "__main__":
    unittest.main()
