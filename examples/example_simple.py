# -*- coding: utf-8 -*-
"""
Example of using lpmodes package to find modes supported by fibre and create 
and display an example mode.

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
core_radius = 15
n_core = 1.43
n_cladding = 1.42
wavelength = 0.5


# We first find the modes supported by this fibre at the specified wavelength.
# This returns a list of instances of the Mode class
modes = lpmodes.find_modes(core_radius, n_core, n_cladding, wavelength)

print(f"Found {len(modes)} modes, {lpmodes.num_rotated_modes(modes)} with rotations")


# From the list of modes we can extract one with the required l and m
l = 2
m = 3
mode = lpmodes.find_mode(modes, l, m)


# We can then plot this mode
grid_size = 100
max_plot_radius = core_radius * 1.1
cos_plot = mode.plot_amplitude(grid_size, max_plot_radius)
sin_plot = mode.plot_amplitude_rotated(grid_size, max_plot_radius)

# We can use ampcol() to get a nice colormap to display amplitudes. We have
# to set the vmin and vmax correctly so that 0 values are plotted white 
plt.figure(dpi = 200)
plt.title(f"Amplitude, l = {mode.l}, m = {mode.m}")
vLim = np.max(np.abs(cos_plot))
plt.imshow(cos_plot, cmap = lpmodes.ampcol(), vmin = -vLim, vmax = vLim)

plt.figure(dpi = 200)
plt.title(f"Amplitude, rotated, l = {mode.l}, m = {mode.m}")
vLim = np.max(np.abs(sin_plot))
plt.imshow(sin_plot, cmap = lpmodes.ampcol(), vmin = -vLim, vmax = vLim)

# We can also plot the intensity of this mode
plt.figure(dpi = 200)
plt.title(f"Intensity, l = {mode.l}, m = {mode.m}")
plt.imshow(np.abs(cos_plot)**2)
