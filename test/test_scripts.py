#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Copyright (C) 2013-2019 ITB - CNR

    This file is part of isoSegmenter.

    isoSegmenter is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    isoSegmenter is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with isoSegmenter.  If not, see <http://www.gnu.org/licenses/>.


If you use isoSegmenter in your work, please cite this manuscript:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into
    Isochores. Evolutionary Bioinformatics. 2015;11:253-261.
    doi:10.4137/EBO.S27693

Created on Tue Feb 12 11:43:35 2019

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
"""

from __future__ import print_function

import os
import shlex
import unittest
import tempfile
import subprocess

# getting module path
module_path = os.path.dirname(__file__)


class IsoSegmenterTestCase(unittest.TestCase):
    """A class to test isoSegmenter scripts"""

    def setUp(self):
        # create temporary file names
        self.outfile = tempfile.mktemp()
        self.graphfile = tempfile.mktemp()

        self.infile = os.path.join(module_path, "chr21.fa.gz")

    def tearDown(self):
        # clean up stuff if exists
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

        if os.path.exists(self.graphfile):
            os.remove(self.graphfile)

    def test_isosegmenter(self):
        """Test isoSegmenter.py script"""

        cmd = (
            "isoSegmenter.py --infile {0} --outfile {1} --graphfile {2} "
            "--draw_legend --verbose").format(
                    self.infile,
                    self.outfile,
                    self.graphfile)

        cmds = shlex.split(cmd)

        # call script
        status = subprocess.check_call(cmds)

        self.assertEqual(status, 0)


class IsoFamilyTestCase(unittest.TestCase):
    """A class to test isoFamily scripts"""

    def setUp(self):
        # create temporary file names
        self.outfile = tempfile.mktemp()
        self.graphfile = tempfile.mktemp()

        self.indir = module_path

    def tearDown(self):
        # clean up stuff if exists
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

        if os.path.exists(self.graphfile):
            os.remove(self.graphfile)

    def test_isofamily(self):
        """Test isoFamily.py script"""

        cmd = (
            "isoFamily.py --indir {0} --outfile {1} --graphfile {2} "
            "--verbose --regexp test_isochores3_chr21.csv").format(
                    self.indir,
                    self.outfile,
                    self.graphfile)

        cmds = shlex.split(cmd)

        # call script
        status = subprocess.check_call(cmds)

        self.assertEqual(status, 0)


class TileImageTestCase(unittest.TestCase):
    """A class to test isoFamily scripts"""

    def setUp(self):
        # create temporary file names
        self.outfile = tempfile.mktemp(suffix=".png")

        # get a test image
        self.image_file = os.path.join(module_path, "chr21.isochores.png")

    def tearDown(self):
        # clean up stuff if exists
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

    def test_tileimages(self):
        """Test rileImages.py script"""

        cmd = (
            "tileImages.py --image_files {0} {1} -o {2}").format(
                    self.image_file,
                    self.image_file,
                    self.outfile)

        cmds = shlex.split(cmd)

        # call script
        status = subprocess.check_call(cmds)

        self.assertEqual(status, 0)


if __name__ == "__main__":
    unittest.main()
