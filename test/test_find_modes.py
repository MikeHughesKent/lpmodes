# -*- coding: utf-8 -*-
"""
Tests that find_modes returns approx number of modes predicted by V number.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import unittest

import context

import lpmodes


class TestFindModes(unittest.TestCase):

    def setUp(self):
        # Set up some test modes
        self.modes = lpmodes.find_modes(5,1.4,1.38,0.5)

    def test_num_modes(self):
        # Divide by 2 because V number includes polarisation
        est_modes = lpmodes.est_num_modes(5,1.4,1.38,0.5) / 2        
        rot_modes = lpmodes.num_rotated_modes(self.modes)
        self.assertAlmostEqual(est_modes, rot_modes, delta = 10)
        
    
if __name__ == '__main__':
    unittest.main()