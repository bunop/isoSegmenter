# -*- coding: utf-8 -*-
"""


    Copyright (C) 2013 ITB - CNR

    This file is part of ISOfinder.

    ISOfinder is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ISOfinder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ISOfinder.  If not, see <http://www.gnu.org/licenses/>.


Created on Fri May  3 14:49:03 2013

@author: Paolo Cozzi <paolo.cozzi@tecnoparco.org>

"""

__all__ = ["Elements", "Utility", "Graphs" ]
__author__ = "Paolo Cozzi <paolo.cozzi@tecnoparco.org>"

import Graphs
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

#The true type font used by EnlargeLabels(Graphs module)
graph_font_type = "/usr/share/fonts/truetype/freefont/FreeSerifBold.ttf"