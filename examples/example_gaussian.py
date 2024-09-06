# -*- coding: utf-8 -*-
"""
Demonstrates simulation of coupling a Gaussian beam into a multimode fibre.

A fibre is modelled, and a Gaussian field is generated. The Gaussian field
is coupled into the fibre, and this is then propagated along the fibre, with
the intensity plotted at intervals.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import sys
import os

import matplotlib.pyplot as plt
import numpy as np

# Necessary if lpmodes hasn't been installed
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import lpmodes

# Define the fibre characteristics and wavelength
nCore = 1.40
nCladding = 1.38
wavelength = 0.5     # microns
coreRadius = 10      # microns
fibreLength = 200       # microns

na = lpmodes.fibre_na(nCore, nCladding)

est_num_modes = round(lpmodes.est_num_modes(coreRadius, nCore, nCladding, wavelength) / 2)
print(f"Estimating {est_num_modes} LP modes from V number (not incl. polarisation).")

     
# Define the plot area
max_plot_radius = coreRadius * 2  # Lets us see the power in the cladding
grid_size = 200                     # pixels

# Find all the LP modes
modes = lpmodes.find_modes(coreRadius, nCore, nCladding, wavelength)

print(f"Found {len(modes)} modes, {lpmodes.num_rotated_modes(modes)} with rotations")

# Calculate the 2D field amplitudes of the LP modes
solution = lpmodes.Solution(modes, grid_size, max_plot_radius)

# Generate a Gaussian field
g_field = lpmodes.gaussian_field(grid_size, max_plot_radius * 2, (10,0), 5)

# Couple into fibre
solution.couple_field(g_field)

# Get an array containing the power in each mode
coupling_eff = round(solution.power_coupled, 3)
print(f"Total coupled power: {coupling_eff * 100} %")

# And then propagate the modes along the fibre
solution.propagate(fibreLength, rotations = False)
  
# Plot Gaussian field
plt.figure(dpi=150); plt.title(f"Input Gaussian power")
plt.imshow(np.abs(g_field)**2)

# Plot amplitude at proximal end
plt.figure(dpi=150); plt.title(f"Intensity Pattern at input")
plt.imshow(np.abs(solution.get_in_field())**2)

# Plot intensity at distal end
plt.figure(dpi=150); plt.title(f"Intensity Pattern after {fibreLength} microns")
plt.imshow(solution.get_out_power())

