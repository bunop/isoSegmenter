# -*- coding: utf-8 -*-
"""

Created on Fri May 31 17:07:40 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.it>

"""

import sys
import unittest

sys.path.append("..")
import GClib.Elements

class test_CalcClass(unittest.TestCase):
    def setUp(self):
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
        """Testing class assignment..."""
        
        for GClevel, Class in self.GClevels:
            result = GClib.Elements.CalcClass(GClevel)
            self.assertEqual(result, Class)
            
if __name__ == "__main__":
    unittest.main()