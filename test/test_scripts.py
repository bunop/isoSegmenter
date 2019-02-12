#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 11:43:35 2019

@author: Paolo Cozzi <paolo.cozzi@ptp.it>
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
    """A class to test iSosegmenter scripts"""

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


if __name__ == "__main__":
    unittest.main()
