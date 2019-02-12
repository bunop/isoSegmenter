#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 19:18:15 2019

@author: Paolo Cozzi <paolo.cozzi@ptp.it>
"""

import os

# the default window size in bp
WINDOW_SIZE = 100000

# The gap tolerance in bp
GAP_TOLERANCE = 5000

# The GClevel upper limit for each class
CLASS_TO_LEVEL = {"L1": 37, "L2": 41, "H1": 46, "H2": 53, "H3": 100}

# The true type font used by EnlargeLabels(Graphs module)
module_path = os.path.dirname(__file__)
graph_font_type = os.path.join(module_path, "FreeSerifBold.ttf")

# The minimum size of an isochore
ISO_MIN_SIZE = 2
