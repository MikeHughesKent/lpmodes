# -*- coding: utf-8 -*-
"""
Demonstrates simulation of a speckle pattern produce by multimode fibre.

A fibre is modelled, and random amplitudes and phases are set for each
mode. This is then propagated along the fibre, and the 
field and intensity recovered at the far end by summing the modes.

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
nCore = 1.44
nCladding = 1.42
wavelength = 0.5     # microns
coreRadius = 15      # microns
fibreLength = 10000  # microns

na = lpmodes.fibre_na(nCore, nCladding)

est_num_modes = round(lpmodes.est_num_modes(coreRadius, nCore, nCladding, wavelength) / 2)
print(f"Estimating {est_num_modes} LP modes from V number (not incl. polarisation).")

     
# Define the plot area
max_plot_radius = coreRadius * 1.2  # Lets us see the power in the cladding
grid_size = 100                     # pixels

# Find all the LP modes
modes = lpmodes.find_modes(coreRadius, nCore, nCladding, wavelength)

print(f"Found {len(modes)} modes, {lpmodes.num_rotated_modes(modes)} with rotations")

# Calculate the 2D field amplitudes of the LP modes
solution = lpmodes.Solution(modes, grid_size, max_plot_radius)

# First we couple into each mode with random complex amplitude
solution.set_random_amplitudes()

# And then propagate the modes along the fibre
solution.propagate(fibreLength, rotations = False)
  
# Plot amplitude at proximal end
plt.figure(dpi=150); plt.title(f"Amplitude Pattern at input")
plt.imshow(np.abs(solution.get_in_field()))

# Plot phase at proximal end
plt.figure(dpi=150); plt.title(f"Phase Pattern at input")
plt.imshow(np.angle(solution.get_in_field()))

# Plot intensity at proximal end
plt.figure(dpi=150); plt.title(f"Intensity Pattern at input")
plt.imshow(solution.get_in_intensity())

# Plot amplitude at distal end
plt.figure(dpi=150); plt.title(f"Amplitude Pattern after {fibreLength} microns")
plt.imshow(np.abs(solution.get_out_field()))

# Plot intensity at distal end
plt.figure(dpi=150); plt.title(f"Intensity Pattern after {fibreLength} microns")
plt.imshow(solution.prop_intensity)

# Plot phase at distal end
plt.figure(dpi=150); plt.title(f"Phase Pattern at after {fibreLength} microns")
plt.imshow(np.angle(solution.get_in_field()))