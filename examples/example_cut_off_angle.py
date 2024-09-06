# -*- coding: utf-8 -*-
"""
Demonstrates use of lpmodes package to determine power coupled into fibre
from tilted plane wave at various angles.

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
core_radius = 10
n_core = 1.435
n_cladding = 1.415
wavelength = 0.5
grid_size = 100
max_plot_radius = core_radius * 2.5

est_cut_off = np.round(lpmodes.acceptance_angle(n_core, n_cladding),1)
print(f"Estimated cut-off angle: {est_cut_off} degrees")

# Find the allowed modes
modes = lpmodes.find_modes(core_radius, n_core, n_cladding, wavelength)
print(f"Found {len(modes)} modes, {lpmodes.num_rotated_modes(modes)} with rotations")

# When we create an instance of Solution, this plots all the modes
solution = lpmodes.Solution(modes, grid_size, max_plot_radius)

# We will couple in beams tilted at these angles
tilts = np.arange(0, 31, 1)

coupledPower = np.zeros((len(tilts)))

for idx, tilt in enumerate(tilts):
       
    # Generate beam with required tilt
    tiltRad = tilt * np.pi / 180
    inField = lpmodes.tilted_field(grid_size, max_plot_radius * 2, wavelength, tiltRad)

    # Couple beam into fibre
    solution.couple_field(inField)
    
    # Total power coupled at this angle
    coupledPower[idx] = np.sum(solution.mode_power)
    
# Normalise so that couple power is 1 for 0 degrees
coupledPower = coupledPower / np.max(coupledPower)

plt.figure(dpi = 150)
plt.title("Power coupled into fibre")
plt.xlabel("Angle (degrees)")
plt.ylabel("Relative power coupled")
plt.plot(tilts, np.sqrt(coupledPower))
plt.axline( (est_cut_off, 0), (est_cut_off, 1), ls = "--")
