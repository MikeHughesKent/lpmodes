# -*- coding: utf-8 -*-
"""
Tests CSV writing and reading functions of lpmodes

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import os
import context

from lpmodes import *

import unittest

class TestCSVFunctions(unittest.TestCase):

    def setUp(self):
        # Set up some test modes
        self.modes = [
            Mode(1.0, 2.0, 0, 1, 1.5, 2, 4.5, 1.55, 1.45, 1.0),
            Mode(1.1, 2.1, 1, 2, 1.6, 2, 4.6, 1.56, 1.46, 1.1)
        ]
        self.test_file = 'test_modes.csv'

    def tearDown(self):
        # Clean up: remove the file after the test
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_modes_to_and_from_csv(self):
        
        # Write modes to CSV
        modes_to_csv(self.modes, self.test_file)

        # Read modes back from CSV
        loaded_modes = modes_from_csv(self.test_file, Mode)

        # Check if the loaded modes are the same as the original modes
        self.assertEqual(self.modes, loaded_modes)
        

if __name__ == '__main__':
    unittest.main()