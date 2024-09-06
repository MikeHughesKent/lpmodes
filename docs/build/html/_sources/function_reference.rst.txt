----------------------------------
Function Reference
----------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^
Classes
^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: Mode()

Stores properties of a single mode, and has methods for plotting mode. See `Mode class <mode.html>`_ for details.

.. py:function:: Solution()

Stores a complete fibre solution, including modes and plot of modes, with various methods
for coupling in beams and propagating along the fibre. See `Solution class <solution.html>`_ for details.



^^^^^^^^^^^^^^^^^^^^^^^^^
Mode Finding and Handling
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.find_modes

.. autofunction:: lpmodes.power_in_core

.. autofunction:: lpmodes.find_mode

.. autofunction:: lpmodes.find_mode_idx

.. autofunction:: lpmodes.num_rotated_modes




^^^^^^^^^^^^^^^^^^^^^^^^^
Fibre Properties
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.est_num_modes

.. autofunction:: lpmodes.fibre_na

.. autofunction:: lpmodes.acceptance_angle


^^^^^^^^^^^^^^^^^^^^^^^^^
Mode and Field Plotting
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.plot_modes_amplitude

.. autofunction:: lpmodes.plot_modes_intensity

.. autofunction:: lpmodes.plot_field

.. autofunction:: lpmodes.vis_all_modes

.. autofunction:: lpmodes.ampcol


^^^^^^^^^^^^^^^^^^^^^^^^^
Coupling and Propagation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.couple_beam

.. autofunction:: lpmodes.propagate_through_fibre

.. autofunction:: lpmodes.mode_coupling_plot



^^^^^^^^^^^^^^^^^^^^^^^^^
Field Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.tilted_field

.. autofunction:: lpmodes.gaussian_field


^^^^^^^^^^^^^^^^^^^^^^^^^
Import and Export
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lpmodes.modes_to_csv

.. autofunction:: lpmodes.modes_from_csv




