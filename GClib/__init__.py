# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:49:03 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

"""

__all__ = ["Elements", "Utility" ]
__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

import Utility
import Elements

#default level of vebosity
logger = Utility.Logger(threshold=1)

#the default window size in bp
WINDOW_SIZE = 100000

#The gap tolerance in bp
GAP_TOLERANCE = 5000

#The GClevel upper limit for each class
CLASS_TO_LEVEL = { "L1":37, "L2":41, "H1":46, "H2":53, "H3":100 }
