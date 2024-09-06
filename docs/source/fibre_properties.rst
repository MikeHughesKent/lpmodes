-----------------
Fibre Properties
-----------------

The package contains several utility functions determining properties of the 
fibre. If the core and cladding refractive indices are defined, for example::
    
    n_core = 1.4
    n_cladding = 1.38
    
then the fibre numerical aperture is given by::

    na = lpmodes.fibre_na(n_core, n_cladding)
    
and the cut-off angle (largest coupled incident ray angle) by::

    cut_off = lpmodes.acceptance_angle(n_core, n_cladding)    
         
The V number of the fibre at a specific wavelength is given by, for example::

    core_radius = 10  # microns
    wavelength = 0.5 # microns
    V = lpmodes.v_number(core_radius, n_core, n_cladding, wavelength)
    
An estimate of the number of modes supported (from the V number) is found using::

    num_modes = lpmodes.est_num_modes(core_radius, n_core, n_cladding, wavelength)
    
This includes two linear polarization directions as well as the sine and cosine 
orientations of modes with l > 0. The number returned by ``est_num_modes``
should therefore be approximately twice the number of solved modes including 
rotations, i.e. that returned by ``num_rotated.modes()``.
