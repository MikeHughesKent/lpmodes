# -*- coding: utf-8 -*-
"""
Example of using lpmodes package to find modes supported by fibre and create 
a single plot showing them all.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import sys
import os

import matplotlib.pyplot as plt
import numpy as np

# Necessary if lpmodes hasn't been installed
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import lpmodes

# We define the parameters of the fibre and the light wavelength
core_radius = 5
n_core = 1.43
n_cladding = 1.42
wavelength = 0.5

# Find the allowed modes
modes = lpmodes.find_modes(core_radius, n_core, n_cladding, wavelength)

print(f"Found {len(modes)} modes, {lpmodes.num_rotated_modes(modes)} with rotations")

# Generate amplitude plot
lpmodes.vis_all_modes(modes)

# Generate intensity plot
lpmodes.vis_all_modes(modes, intensity = True)
