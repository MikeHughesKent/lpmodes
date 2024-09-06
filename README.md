# lpmodes: Python Package for Simulation Multimode Fibres

The lpmodes package is a toolkit for finding, simulating, 
and manipulating linearly polarized (LP) modes in step index multimode optical fibers. 

Functions are provided to determine the allowed modes, compute their properties, 
and simulate coupling and propagation of beams.

A GUI, lpmodes_gui, allows modes to be solved, visualised and exported.

![Example modes](https://github.com/MikeHughesKent/lpmodes/raw/main/docs/source/modes.png)

Full documentation is on [Read the docs](https://lpmodes.readthedocs.io/en/latest/index.html). 

Also see the examples in the [examples folder](https://github.com/MikeHughesKent/lpmodes/tree/main/examples).


## Package Features

Functions are provided to:
* Find all allowed modes for step index fibres (weakly guided approximation)
* Export mode parameters as CSV
* Calculate fraction of power in core for all modes
* Plot amplitude or intensity of modes as numpy arrays
* Generate a single plot of all modes as a matplotlib figure
* Simulate coupling arbitary beams into the fibre (functions for generation of plane waves and Gaussian beams are provided)
* Calculate coupling efficiency
* Propagate coupled fields along a length of fibre, and visualise output amplitude and intensity


## GUI

A graphical user interface allow modes to be solved, visualised and exported.

![lpmodes GUI screenshot](https://github.com/MikeHughesKent/lpmodes/raw/main/docs/source/modes.png)

Parameters of the fibre (NAs, core radius) and wavelength are defined in the 
top left. Clicking 'Find Modes' then finds solutions. The table in the centre
is updated with the mode details. Power in core is not calculated immediately 
as this is more time consuming, click 'Find Power in Core' button to update
this. Select a line of the table to display an image of that mode on the 
right. There are options on the left to choose the image grid size, the 
physical size of the image as a multiple of the core size, whether to display
the amplitude or the intensity, and whether to display a circle showing the 
core size. The data can be exported as CSV, or the current mode image or all
mode images saved as tifs using buttons on the left.


## Contributing

Development is currently mainly by [Mike Hughes](https://research.kent.ac.uk/applied-optics/hughes/)' lab in the 
[Applied Optics Group](https://research.kent.ac.uk/applied-optics), Physics & Astronomy, University of Kent. 
I'm happy to collaborate with academic users to help your use case, and if you would like help using lpmodes for 
commercial purposes, consultancy is available, please contact [Mike Hughes](mailto:m.r.hughes@kent.ac.uk). 

Help testing and developing the package is welcome, and I'm happy for this to become a more collaborate effort and share credit,
please [get in touch](mailto:m.r.hughes@kent.ac.uk).


## Requirements

* Numpy
* Scipy
* Matplotlib
* PyQT5 (GUI only)
* PIL (GUI only)
