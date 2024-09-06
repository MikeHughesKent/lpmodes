# -*- coding: utf-8 -*-
"""
Demonstrates coupling a speckle pattern into a multimode fibre.

A fibre is modelled, and a random input field with different
Gaussian smoothing filters applied (to change the speckle grain size)
is injected into the fibre.

The field is represented in the fibre modes and propagated along the fibre, the 
intensity recovered at the far end is displayed. The change in the speckle
pattern exiting the fibre can be seen as the input speckle grain size
increases.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import sys
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

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

# Generate a random speckle pattern
field_amp = np.random.rand(grid_size, grid_size)
field_phase =  np.random.rand(grid_size, grid_size) * 2 * np.pi
field = field_amp * np.exp(1j * field_phase)

# Trial sizes of the smoothing filter used to change the speckle grain size
sigmas = np.arange(0.5, 3, 0.5)

fig, axs = plt.subplots(len(sigmas), 2, figsize=(10, 5 * len(sigmas)))

for i, sigma in enumerate(sigmas):
    
    field_filtered = gaussian_filter(np.real(field), sigma) + 1j * gaussian_filter(np.imag(field), sigma) 

    # Plot intensity of speckle pattern
    axs[i, 0].imshow(np.abs(field_filtered)**2, cmap='gray')
    axs[i, 0].set_title("Input speckle pattern intensity")
    axs[i, 0].set_ylabel(f"σ = {sigma} px", fontsize=12, rotation=90, labelpad=10, ha='center')

    # Couple the field into the fibre
    solution.couple_field(field_filtered)

    # And then propagate the modes along the fibre
    solution.propagate(fibreLength, rotations = False)
    
    # Plot intensity patterns at the end of the fibre
    axs[i, 1].imshow(solution.prop_intensity, cmap='gray')
    axs[i, 1].set_title("Output intensity")
 
