# -*- coding: utf-8 -*-
"""
Tests that find_modes returns modes that are correctly normalised.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import unittest

import numpy as np

import context

import lpmodes



class TestNorm(unittest.TestCase):

   

    def test_norm(self):
        
        n_core = 1.40
        n_cladding = 1.38
        wavelength = 0.5      # microns
        core_radius = 5       # microns
        grid_size = 100
        max_plot_radius = core_radius * 2
        
        # Find all the LP modes and plot them
        modes = lpmodes.find_modes(core_radius, n_core, n_cladding, wavelength)
        solution = lpmodes.Solution(modes, grid_size, max_plot_radius)

        
        self.assertAlmostEqual(1, np.sum(solution.mode_sin[0] * solution.mode_sin[0]), places = 4)
        self.assertAlmostEqual(1, np.sum(solution.mode_sin[23] * solution.mode_sin[23]), places = 4)
        self.assertAlmostEqual(0, np.sum(solution.mode_sin[23] * solution.mode_cos[23]), places = 4)
        self.assertAlmostEqual(0, np.sum(solution.mode_sin[20] * solution.mode_sin[23]), places = 4)
        
    
if __name__ == '__main__':
    unittest.main()