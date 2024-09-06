# -*- coding: utf-8 -*-
"""
Tests the find_mode_idx function.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import unittest

import context

from lpmodes import *


class TestFindModeIdx(unittest.TestCase):

    def setUp(self):
        # Set up some test modes
        self.modes = find_modes(5,1.4,1.38,0.5)

    def test_find_existing_mode(self):
        # Test finding an existing mode
        idx = find_mode_idx(self.modes, 1, 2)
        self.assertEqual(1, self.modes[idx].l)
        self.assertEqual(2, self.modes[idx].m)

    def test_mode_not_found(self):
        # Test behavior when the mode is not found
        idx = find_mode_idx(self.modes, 300, 3)
        self.assertIsNone(idx)

    def test_empty_modes_list(self):
        # Test behavior with an empty modes list
        idx = find_mode_idx([], 0, 1)
        self.assertIsNone(idx)
        

if __name__ == '__main__':
    unittest.main()