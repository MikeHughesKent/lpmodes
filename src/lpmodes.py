# -*- coding: utf-8 -*-
"""
lpmodes

Package for working with optical fibre LP modes, including finding all allowed
modes, plotting, and simple manipulations.

@author: Mike Hughes, Applied Optics Group, University of Kent
"""

import math
import csv

import matplotlib.pyplot as plt
import numpy as np

from scipy.special import jv, kn
from scipy.optimize import minimize_scalar


class Mode():
    """ Class used to store parameters of a mode, and to plot it.
    
        Arguments:
            u           : float
                          Radial propagation constant in the core.
            w           : float
                          Radial propagation constant in the cladding.
            l           : int
                          Azimuthal mode number.
            m           : int
                          Radial mode number.
            beta        : float
                          Propagation constant of the mode.
            n_eff       : float
                          effective refractive index of the mode              
            core_radius : float
                          Radius of the fibre core.
            wavelength  : float
                          Wavelength of the light in the mode.
            n_core      : float
                          Refractive index of the fibre core.
            n_cladding  : float
                          Refractive index of the cladding.

    """

    def __init__(self, u, w, l, m, beta, n_eff, core_radius, wavelength, n_core, n_cladding):
        
        self.u = u
        self.w = w
        self.l = l
        self.m = m
        self.beta = beta
        self.n_eff = n_eff
        self.core_radius = core_radius
        self.wavelength = wavelength
        self.n_core = n_core
        self.n_cladding = n_cladding
        
        
        
    def __eq__(self, other):
        """ Returns True if this mode is identical to mode 'other'.
        """
        
        if not isinstance(other, Mode):
            return False
        return (self.u == other.u and self.w == other.w and self.l == other.l and
                self.m == other.m and self.beta == other.beta and
                self.n_eff == other.n_eff and
                self.core_radius == other.core_radius and
                self.wavelength == other.wavelength and self.n_core == other.n_core and
                self.n_cladding == other.n_cladding)    
        

    def plot_amplitude(self, grid_size, max_plot_radius):
        """
        Returns 'plot' of amplitude of mode as a 2D numpy array. This function
        returns the 'cosine' version.

        Arguments:       
              grid_size       : int
                                size of square plot in pixels
              max_plot_radius : float
                                radial distance at edge of plot, microns                  

        Returns
              np.array        : 2D array of floats, plot of mode, cosine version  
        """

        # Find centre of grid
        gridCentre = grid_size // 2

        # Initialise output array
        mode_cos = np.zeros((grid_size, grid_size))

        # Calculate grid points
        gridPoints = np.arange(0, grid_size, 1) - gridCentre
        xMesh, yMesh = np.meshgrid(gridPoints, gridPoints)

        # Convert  to polar co-ordinates
        angle = np.arctan2(yMesh, xMesh)
        rad = np.sqrt(xMesh**2 + yMesh**2)
        rad = rad / (grid_size / 2) * max_plot_radius
     
        cosTerm = np.cos(self.l * angle)

        # Calculate the core and cladding 2D functions
        coreBessel = jv(self.l, self.u / self.core_radius *
                        rad) / jv(self.l, self.u)
        claddingBessel = kn(self.l, self.w / self.core_radius * rad) / kn(self.l, self.w)

        # Work out which grid points are in core and which in cladding
        in_core = (rad <= self.core_radius) * (rad <= max_plot_radius)
        in_cladding = (rad > self.core_radius) * (rad <= max_plot_radius)

        # Calculate sin and cos versions of core field
        mode_cos[in_core] = coreBessel[in_core] * cosTerm[in_core]

        # Calculate sin and cos versions of cladding field
        mode_cos[in_cladding] = claddingBessel[in_cladding] * cosTerm[in_cladding]

        # Normalise to have a total power of 1
        cos_norm = np.sqrt(np.sum(np.abs(mode_cos) ** 2 ))
        if cos_norm != 0:
            mode_cos = mode_cos / cos_norm
       
        return mode_cos
    

    def plot_amplitude_rotated(self, grid_size, max_plot_radius):
        """
        Returns 'plot' of amplitude of mode as a 2D numpy array. This function
        returns the 'sine' version, i.e. rotated by 4pi/l w.r.t. plot_amplitude.

        Arguments:       
              grid_size       : int
                                size of square plot in pixels
              max_plot_radius : float
                                radial distance at edge of plot, microns     

        Returns
              np.array        : 2D array of floats, plot of mode, sine version  
        """

        # For l = 0 the rotated version is the same so just return that
        if self.l == 0:
            return self.plot_amplitude(grid_size, max_plot_radius)
       
        # Find centre of grid
        gridCentre = grid_size // 2

        # Initialise output array
        mode_sin = np.zeros((grid_size, grid_size))

        # Calculate grid points
        gridPoints = np.arange(0, grid_size, 1) - gridCentre
        xMesh, yMesh = np.meshgrid(gridPoints, gridPoints)

        # Convert  to polar co-ordinates
        angle = np.arctan2(yMesh, xMesh)
        rad = np.sqrt(xMesh**2 + yMesh**2)
        rad = rad / (grid_size / 2) * max_plot_radius

        sinTerm = np.sin(self.l * angle)

        # Calculate the core and cladding 2D functions
        coreBessel = jv(self.l, self.u / self.core_radius *
                        rad) / jv(self.l, self.u)
        claddingBessel = kn(
            self.l, self.w / self.core_radius * rad) / kn(self.l, self.w)

        # Work out which grid points are in core and which in cladding
        in_core = (rad <= self.core_radius) * (rad <= max_plot_radius)
        in_cladding = (rad > self.core_radius) * (rad <= max_plot_radius)

        # Calculate sin and cos versions of core field
        mode_sin[in_core] = coreBessel[in_core] * sinTerm[in_core]

        # Calculate sin and cos versions of cladding field
        mode_sin[in_cladding] = claddingBessel[in_cladding] * sinTerm[in_cladding]

        # Normalise to total of 1
        sin_norm = np.sqrt(np.sum(np.abs(mode_sin)**2))

        if sin_norm != 0:
            mode_sin = mode_sin / sin_norm
            
        return mode_sin


    def plot_intensity(self, grid_size, max_plot_radius):
        """
        Returns 'plot' of intensity of mode as a 2D numpy array. This function
        returns the 'cosine' version.

        Arguments:       
              grid_size       : int
                                size of square plot in pixels
              max_plot_radius : float
                                radial distance at edge of plot, microns     

        Returns
             np.array        : 2D array of floats, plot of mode, cosine version  
        """
       
        return np.abs(self.plot_amplitude(grid_size, max_plot_radius))**2
   
    
    def plot_intensity_rotated(self, grid_size, max_plot_radius):
        """
        Returns 'plot' of intensity of mode as a 2D numpy array. This function
        returns the 'sine' version.
    
        Arguments:       
              grid_size       : int
                                size of square plot in pixels
              max_plot_radius : float
                                radial distance at edge of plot, microns     
    
        Returns
              np.array        : 2D array of floats, plot of mode, sine version  
        """
       
        return np.abs(self.plot_amplitude_rotated(grid_size, max_plot_radius))**2   


class Solution:
    """ Class used to store a set of lp modes for a fibre and to perform operations
    that are applied to the whole set of modes.
    
    Arguments:
        modes           : list
                          list of instances of Mode
                          as generated by find_modes
                          
    Keyword Arguments:                      
        grid_size       : int, optional
                          The size of the grid for mode plotting. Defaults to 100.
        max_plot_radius : float, optional
                          The maximum  physical radius for plotting. If 
                          None, it is set to 1.5 times the core radius.
                          Defaults to None.

    Attributes:
        modes          : list
                         A list of Mode instances, as generated by find_modes 
        num_modes      : int
                         The number of modes in the 'modes' list.
        mode_l         : np.array
                         A 1D array storing the 'l' value for each mode.
        mode_m         : np.array
                         A 1D array storing the 'm' value for each mode.
        mode_cos       : np.array
                         A 3D array storing a plot of the amplitude of the 
                         cosine component of each mode (mode, x, y)
        mode_sin       : np.array
                         A 3D array storing a plot of the amplitude of the 
                         sine component of each mode (mode, x, y)
        mode_cos_intensity : np.array
                             A 3D array storing a plot of the intensity of the 
                             cosine component of each mode (mode, x, y)
        mode_sin_intensity : np.array
                             A 3D array storing a plot of the amplitude of the 
                             sine component of each mode (mode, x, y)
        cos_amp        : np.array
                         A 1D array storing the complex amplitude coupled into 
                         the cosine version of each mode.
        sin_amp        : np.array
                         A 1D array storing the complex amplitude coupled into 
                         the sine version of each mode.
        mode_power     : np.array
                         1 1D array storing the power coupled into each mode
                         (total power in both orientations)
        power_coupled  : float
                         Total power coupled into fibre
        prop_cos_amp   : np.array
                         A 1D array storing the complex amplitude coupled into 
                         the cosine version of each mode after propagation to
                         end of fibre.
        prop_sin_amp   : np.array
                         A 1D array storing the complex amplitude coupled into 
                         the sine version of each mode after propagation to
                         end of fibre.
        prop_field     : np.array
                         2D complex array containing plot of field after
                         propagation to end of fibre
        prop_intensity : np.array
                         2D array containing plot of intensity after
                         propagation to end of fibre
    """
    
    def __init__(self, modes, grid_size = 100, max_plot_radius = None):
        """
        Initializes the Solution object with a set of modes, plots and
        stores the modes.

        Arguments:
            modes           : list
                              list of instances of Mode
                              as generated by find_modes
                              
        Keyword Arguments:                      
            grid_size       : int, optional
                              The size of the grid for mode plotting. Defaults to 100.
            max_plot_radius : float, optional
                              The maximum  physical radius for plotting. If 
                              None, it is set to 1.5 times the core radius.
                              Defaults to None.
        """
        
        # If max_plot_radius no supplied, use a sensible default
        if max_plot_radius is None:
            max_plot_radius = modes[0].core_radius * 1.5
         
        self.modes = modes
        self.num_modes = len(modes)
       
        # Store this in a more convenient way
        for mode in modes:
            self.mode_l = mode.l
            self.mode_m = mode.m

        # Plot the modes
        self.mode_cos, self.mode_sin = plot_modes_amplitude(modes, grid_size, max_plot_radius)
        self.mode_cos_intensity = np.abs(self.mode_cos)**2
        self.mode_sin_intensity = np.abs(self.mode_sin)**2
        
        # Initially set all the coupled amplitudes to 0, this will changes
        # later if we couple in a beam
        self.cos_amp = np.zeros(self.num_modes, dtype = 'complex')
        self.sin_amp = np.zeros(self.num_modes, dtype = 'complex')
        self.mode_power = np.zeros(self.num_modes)
        self.power_coupled = 0
        self.prop_field = 0
        self.prop_intensity = 0
        self.prop_cos_amp = np.zeros(self.num_modes, dtype = 'complex')
        self.prop_sin_amp = np.zeros(self.num_modes, dtype = 'complex')
        
         
    def set_amplitudes(self, cos_amp, sin_amp):
        """ Sets the amplitude coupled into each mode.
        
        Arguments:
            cos_amp     : np.array
                          1D array of floats, amplitude for cosine version of mode
            sin_amp     : np.array
                          1D array of floats, amplitude for sine (i.e. rotated) version of mode
        """
        self.cos_amp = cos_amp
        self.sin_amp = sin_amp
        
    
    def set_random_amplitudes(self):
        """Sets random amplitudes for the cosine and sine components of each mode."""
        
        mode_amplitude = np.random.rand(len(self.modes)) 
        mode_angle = np.random.rand(len(self.modes)) * np.pi * 2
        mode_phase = np.random.rand(len(self.modes)) * np.pi * 2
        self.cos_amp = np.real(mode_amplitude * np.exp(1j * mode_angle)) * np.exp(1j * mode_phase)
        self.sin_amp = np.imag(mode_amplitude * np.exp(1j * mode_angle)) * np.exp(1j * mode_phase)
    
        
    def couple_field(self, field):
        """Couples an input field into the modes and calculates the resulting mode amplitudes.

        Arguments:
           field      : np.array
                        The input field to be coupled into the modes, as a 2D
                        complex numpy array. This must be on a grid that is the same
                        dimensions and with the same pixel size as defined by
                        the grid_size and max_plot_radius argument to __init__.
        """
        self.cos_amp, self.sin_amp, self.mode_power = couple_beam(field, self.modes, self.mode_cos, self.mode_sin)
        self.power_coupled = np.sum(self.mode_power)
    
    def propagate(self, distance, rotations = False): 
        """Propagates the modes through a fibre over a given distance, storing
        the new complex amplitudes.

        Arguments:
            distance  : float
                        The distance over which the modes are propagated.
        
        Keyword Arguments:
            rotations : bool, optional
                        If True, applies random rotations of each mode during 
                        propagation. Defaults to False.
        """
        self.prop_field, self.prop_intensity, self.prop_cos_amp, self.prop_sin_amp = propagate_through_fibre(self.modes, self.mode_cos, self.mode_sin, self.cos_amp, self.cos_amp, distance, rotations = rotations)


    def get_in_amplitudes(self):
        """Returns the input (i.e. before propagation) coupled amplitudes for 
        the cosine and sine components of each mode.
        
        Returns:
            tuple : (np.array, np.array)
                    A tuple containing the two 1D numpy arrays containng the
                    input cosine amplitudes and sine amplitudes.
        """
        return self.cos_amp, self.sin_amp
    
   
    def get_out_amplitudes(self):
        """Returns the output amplitudes for the cosine and sine components 
        after propagation.

        Returns:
            tuple : (np.array, np.array)
                    A tuple containing the two 1D numpy arrays containng the
                    input cosine amplitudes and sine amplitudes.
        """
        return self.prop_cos_amp, self.prop_sin_amp
    
     
    def get_in_field(self):
        """Calculates the input field (i.e. before propagation) as represented 
        as a superposition of modes with their respective input amplitudes.

        Returns:
            np.array : The input field as a 2D complex numpy array.
        """
        return np.sum(np.expand_dims(self.cos_amp, (1,2)) * self.mode_cos + np.expand_dims(self.sin_amp, (1,2)) * self.mode_sin,0)
    
    
    def get_out_field(self):
        """Calculates the input power (intensity) distribution 
        (i.e. before propagation) as represented as a superposition of modes 
        with their respective input amplitudes, all squared.

        Returns:
            np.array : The output power (intensity map) as a 2D array.
        """
        return np.sum(np.expand_dims(self.prop_cos_amp, (1,2)) * self.mode_cos + np.expand_dims(self.prop_sin_amp, (1,2)) * self.mode_sin,0)
        
    
    def get_in_intensity(self):
        """Calculates the input power (intensity) distribution 
        (i.e. before propagation) as represented as a superposition of modes 
        with their respective input amplitudes, all squared.

        Returns:
            np.array : The input power (intensity map) as a 2D array.
        """
        return np.abs(self.get_in_field())**2
    
    
    def get_out_intensity(self):
        """Calculates the output power (intensity) distribution 
        (i.e. after propagation) as represented as a superposition of modes 
        with their respective output amplitudes, all squared.

        Returns:
            np.array : The output power (intensity map) as a 2D array.
        """
        
        return self.prop_intensity
    
    
    def plot_amplitude(self, l, m):
        """ Returns an amplitude plot of the sine and cosines components of the 
        mode with the specified l and m. Returns None, None if mode is not present.
        
        Arguments:
            l     : int
                    azimuthal mode number
            m     : int
                    radial mode number
        
        Returns:
            np.array : 2D numpy array containing plot of mode amplitude. Will be
                       None if mode is not present.
        """       
        idx = find_mode_idx(self.modes, l, m)
        if idx is not None:
            return self.mode_cos[idx], self.mode_sin[idx]
        else:
            return None, None 
        
        
    def plot_intensity(self, l, m):
        """ Returns an intensity plot of the sine and cosines components of the 
        mode with the specified l and m. Returns None, None if mode is not present.
        
        Arguments:
            l     : int
                    azimuthal mode number
            m     : int
                    radial mode number
        
        Returns:
            np.array : 2D numpy array containing plot of mode intensity. Will be
                       None if mode is not present.
        """       
        idx = find_mode_idx(self.modes, l, m)
        if idx is not None:
            return self.mode_cos_intensity[idx], self.mode_sin_intensity[idx]
        else:
            return None, None     
        
        
    def plot_mode_coupling(self):
        """ Returns a plot, as a 2D numpy array, showing the power coupled
        into each mode, with the modes organised by l and m.
        """
        return mode_coupling_plot(self.modes, self.cos_amp, self.sin_amp)


def v_number(core_radius, n_core, n_cladding, wavelength):
    """ 
    Estimates number of modes supported by a fibre using V number. Includes
    polarisation modes, so divide by 2 to get close to number of modes that
    will be found by find_modes (also need to include rotations for modes 
    with m > 1).

    Arguments: 
       core_radius    : float
                        radius of the core in microns
       n_core         : float
                        core refractive index
       n_cladding     : float
                        cladding refractive index
       wavelength     : float
                        wavelength of light in microns

    Returns:
       int            : estimated number of modes
    """

    # Numerical Aperture
    NA = fibre_na(n_core, n_cladding)

    # Fibre V Number
    V = 2 * math.pi / wavelength * core_radius * NA

    return V         
        

def est_num_modes(core_radius, n_core, n_cladding, wavelength):
    """ 
    Estimates number of modes supported by a fibre using V number. Includes
    polarisation modes, so divide by 2 to get close to number of modes that
    will be found by find_modes (also need to include rotations for modes 
    with m > 1).

    Arguments: 
       core_radius    : float
                        radius of the core in microns
       n_core         : float
                        core refractive index
       n_cladding     : float
                        cladding refractive index
       wavelength     : float
                        wavelength of light in microns

    Returns:
       int            : estimated number of modes
    """

    num_modes = (v_number(core_radius, n_core, n_cladding, wavelength)**2)/2

    return np.round(num_modes)


def find_modes(core_radius, n_core, n_cladding, wavelength):
    """
    Finds all LP modes for a step index fibre by finding solutions to the
    equation:

          besselj(l, u) / (u * besselj(l - 1, u))
               = besselk(l, w) / (w * besselk(l - 1, w))

          where w = sqrt(v^2-u^2) and 
                v = (2*pi*core_radius/wavelength)*sqrt(core_n^2-cladding_n^2)  

    Arguments: 
       core_radius    : float
                        radius of the core in microns
       n_core         : float
                        core refractive index
       n_cladding     : float
                        cladding refractive index
       wavelength     : float
                        wavelength of light in microns

    Returns:
       Mode           : instance of Mode class

    """

    modes = []

    # Calculate fibre V number
    v = (2 * math.pi * core_radius / wavelength) * \
        math.sqrt(n_core**2 - n_cladding ** 2)

    # Calculate wavenumber
    k = 2 * math.pi / wavelength

    # Sets the coarse search parameters
    # May need to be larger for very large V fibres.
    fineNPoints = 1000
    fineURange = np.linspace(0, v - .1, fineNPoints)

    # Initialise parameters
    iL = 0
    residual = np.zeros(len(fineURange))
    signChange = np.zeros(len(fineURange))

    while True:

        # Search for approximate solutions where there is a sign change and
        # the second derivative is positive
        for idx, u in enumerate(fineURange):

            residual[idx] = __calculate_LP_mismatch(
                core_radius, n_core, n_cladding, wavelength, iL, u)
            if idx > 1:
                signChange[idx] = (np.sign(residual[idx]) != np.sign(
                    residual[idx-1])) and (residual[idx] > residual[idx-1])

        # Pull out all the sign changes we found
        coarseSolutions = np.argwhere(signChange)

        # If there are no solutions for this l then we are done
        if len(coarseSolutions) == 0:
            return modes

        # Search for exact solutions around each sign change
        for idx, solution in enumerate(coarseSolutions):

            # Search where we know there is a change of sign
            minRange = fineURange[solution - 1]
            maxRange = fineURange[solution]

            # Find the exact point of the change of sign
            def wrapper(arg): return __calculate_LP_mismatch(
                core_radius, n_core, n_cladding, wavelength, iL, arg)

            u = minimize_scalar(wrapper, bounds=(minRange, maxRange)).x[0]
            beta = math.sqrt(k**2 * n_core**2 - (u/core_radius)**2)
            w = math.sqrt(v**2-u**2)
            m = idx + 1
            l = iL
            n_eff = beta / k

            mode = Mode(u, w, l, m, beta, n_eff, core_radius,
                        wavelength, n_core, n_cladding)
            modes.append(mode)

        # Move to the next l index
        iL = iL + 1


def __calculate_LP_mismatch(core_radius, core_n, cladding_n, wavelength, order, u):
    """
    Calculates:

        besselj(order, u) / (u * besselj(order - 1, u)) - 
                        besselk(order, w) / (w * besselk(order - 1, w))

        for the specified value of u, where :
            w = sqrt(v^2-u^2) and 
            v = (2*pi*core_radius/wavelength)*sqrt(core_n^2-cladding_n^2)  

    Arguments:
       core_radius    : float
                        radius of the core in microns
       core_n         : float
                        core refractive index
       cladding_n     : float
                        cladding refractive index
       wavelength     : float
                        wavlength of light in microns
       order          : int
                        mode l
       u              : int
                        mode trial u 

    Returns:
       float          : difference between core and cladding terms, when
                        this is zero a mode has been found
    """

    # Wavenumber
    k = 2 * math.pi / wavelength

    # Beta value
    beta = math.sqrt(k**2 * core_n**2 - (np.squeeze(u)/core_radius)**2)

    # w value (cladding)
    w = core_radius * math.sqrt(beta**2 - k**2 * cladding_n**2)

    # core function
    den = (u * jv(order - 1, u))

    if den != 0:
        coreTerm = jv(order, u) / (u * jv(order - 1, u))
    else:
        coreTerm = 0
        
    # cladding function
    claddingTerm = kn(order, w) / (w * kn(order - 1, w))

    # Calculate mismatch
    residual = coreTerm + claddingTerm

    return residual


def num_rotated_modes(modes):
    """ Returns number of modes when rotations are counted as separate modes.
    
    Arguments:
        modes      : list 
                     list of instances of Mode
                     as returned by find_modes
    
    Returns: 
        int        : number of modes                 
    """

    return len(modes) + sum(1 for mode in modes if mode.l > 0)


def ampcol():
    """ Generates a colour map that is useful for plotting amplitude maps.
    The colour map runs from red to blue, with white in the middle.
    
    Returns:
        LinearSegmentedColormap : colormap 
    
        
    """
    from matplotlib.colors import LinearSegmentedColormap

    colors = [
        (0, 0, 1),  # Red
        (1, 1, 1),  # White
        (1, 0, 0),  # Blue
    ]

    # Define the colormap
    return LinearSegmentedColormap.from_list("red_white_blue", colors)


def find_mode(modes, l, m):
    """
    In a list of modes (instances of Mode) returns the mode with
    the specified l and m values.

    Arguments:
        modes       : list of instances of Mode class
                      modes as returned by find_modes
        l           : int
                      azimuthal mode number
        m           : int
                      radial mode number

    Returns:
        Mode        : instance of Mode class, mode with required l and m 
                      (or None if not present)
    """
    for mode in modes:
        if mode.l == l and mode.m == m:
            return mode
    return None


def find_mode_idx(modes, l, m):
    """
    In a list of modes (instances of Mode) find the index of the mode with
    the specififed l and m values.

    Arguments:
        modes       : list of Mode
                      modes as returned by find_modes
        l           : int
                      azimuthal mode number
        m           : int
                      radial mode number

    Rerturns:
        idx         : int
                      index of required mode (or None if not present)
    """
    for idx, mode in enumerate(modes):
        if mode.l == l and mode.m == m:
            return idx
    
    return None


def plot_modes_amplitude(modes, grid_size, max_plot_radius):
    """
    Returns 'plots' of modes in a list of modes as numpy arrays. Both the 
    cosine and sine versions are returned, as two separate 3D numpy arrays.

    Arguments:
        modes           : list of instances of Mode
                          list of modes, as returned by find_modes       
        grid_size       : int
                          size of plot in pixels (square)
        max_plot_radius : float
                          radial distance at edge of plot, microns
       
    Returns:
        np.array        : 3D array, plots of sine versions of modes, 
                          ordered (mode, x, y) 
        np.array       :  3D array, plots of sine versions of modes, 
                          ordered (mode, x, y)  
    """

    # Initialise output arrays
    mode_sin = np.zeros((len(modes), grid_size, grid_size))
    mode_cos = np.zeros_like(mode_sin)

    for idx, mode in enumerate(modes):

        mode_sin[idx, :, :] = mode.plot_amplitude_rotated(grid_size, max_plot_radius)
        mode_cos[idx, :, :] = mode.plot_amplitude(grid_size, max_plot_radius)

    return mode_cos, mode_sin


def plot_modes_intensity(modes, grid_size, max_plot_radius):
    """
    Returns 'plots' of intensity of modes in a list of modes as numpy arrays. Both the 
    cosine and sine versions are returned, as two separate 3D numpy arrays.

    Arguments:
        modes            : list of instances of Mode class
                           list of modes, as returned by find_modes       
        grid_size        : int
                           size of plot in pixels (square)
        max_plot_radius  : float
                           radial distance at edge of plot, microns
     
    Returns:
        np.array         : 3D float, plots of intensity of modes, orderd (mode, x, y). 
                           Cosine version  
        np.array         : 3D float, plots of intensity of modes, orderd (mode, x, y). 
                           Sine version  
    """

    return plot_modes_amplitude(modes, grid_size, max_plot_radius)[0]**2, plot_modes_amplitude(modes, grid_size, max_plot_radius)[1]**2



def plot_field(mode_cos, mode_sin, mode_cos_coupling, mode_sin_coupling):
    """
    Computes the combined optical field from cosine and sine mode components 
    and their respective couplings.

    Arguments:
        mode_cos           : np.array
                             3D array of plot of cos components of modes 
                             (mode, x, y)
        mode_sin           : list of np.array
                             3D array of plot of sin components of modes 
                             (mode, x, y)
        mode_cos_coupling  : np.array
                             1D Array of coupling coefficients for the cos
                             components of the modes. 
        mode_sin_coupling  : np.array
                             1D Array of coupling coefficients for the sine
                             components of the modes.

    Returns:
        np.array           : The combined optical field as a 2D complex numpy 
                             array (x,y)
    """
    
    mode_cos_coupling = np.expand_dims(mode_cos_coupling, (1,2))
    mode_sin_coupling = np.expand_dims(mode_sin_coupling, (1,2))

    num_modes = np.shape(mode_cos)[0]
    
    field = np.zeros((np.shape(mode_cos)[1], np.shape(mode_cos)[2]), dtype = 'complex')
    for mode in range(num_modes):
        field += (mode_cos[mode] * mode_cos_coupling[mode] + mode_sin[mode] * mode_sin_coupling[mode] )
         
    return field
    

def fibre_na(n_core, n_cladding):
    """ Returns NA of fibre with specified core and cladding refractive index

    Arguments:
        n_core     : float 
                     core refractive index
        n_cladding : float
                     cladding refractive index

    Returns:
        float      : numerical aperture
    """
    return np.sqrt(n_core**2 - n_cladding**2)


def acceptance_angle(n_core, n_cladding):
    """ Returns acceptance angle in degrees of fibre with spcecified core and cladding
    refractive index.

    Arguments:
        n_core     : float
                     core refractive index
        n_cladding : float
                     cladding refractive index

    Returns:
        float      : angle in degrees
    """
    
    return np.arcsin(fibre_na(n_core, n_cladding)) * 180 / np.pi


def tilted_field(grid_size, field_size, wavelength, tilt):
    """ Generates a (square) complex field with a tilt.

    Arguments:
        grid_size     : int
                        number of pixels in each direction 
        field_size    : float
                        real size of field in microns
        wavelength    : float
                        wavelength of light in microns
        tilt          : float
                        tilt of field in radians

    Returns:
        np.array      : 2D complex array containing field                   

    """

    pixel_size = field_size/grid_size

    phaseIncrement = 2 * np.pi * pixel_size * np.tan(tilt) / wavelength

    xM, yM = np.meshgrid(range(1, grid_size+1), range(1, grid_size+1))

    phase_map = xM * phaseIncrement

    field = np.exp(1j * phase_map)

    return field



def gaussian_field(grid_size, field_size, center, sigma):
    """Generates a (square) complex Gaussian field.

    Arguments:
        grid_size     : int
                        number of pixels in each direction 
        field_size    : float
                        real size of field in microns
        center        : tuple of floats
                        (x, y) coordinates of the Gaussian center in microns
        sigma         : float
                        standard deviation of the Gaussian in microns

    Returns:
        np.array      : 2D complex array containing field                   
    """

    # Generate grid of coordinates
    x = np.linspace(-field_size/2, field_size/2, grid_size)
    y = np.linspace(-field_size/2, field_size/2, grid_size)
    xM, yM = np.meshgrid(x, y)

    # Calculate the distance from the center for each point in the grid
    x_center, y_center = center
    distance_squared = ((xM - x_center)**2 + (yM - y_center)**2)

    # Create Gaussian 
    amplitude_map = np.exp(-distance_squared / (2 * sigma**2))    
    intensity_map = amplitude_map**2
    
    # Normalise for total power of 1
    amplitude_map = amplitude_map / np.sqrt(np.sum(intensity_map))    
   
    # Create complex field (assuming zero phase)
    field = amplitude_map * np.exp(1j * 0)

    return field



def couple_beam(field, modes, modes_cos, modes_sin):
    """ Determines power from complex field that will be coupled into each mode.

    Arguments:
        field:      np.array
                    2D complex array representing E-field to be coupled
        modes_cos : np.array
                    plot of cosine version of mode 
        modes_sin : np.array
                    plot of sine version of mode 

    Returns:
        np.array  : 1D array, amplitude in each cosine version
        np.array  : 1D array, amplitude in each sine version
        np.array  : 1D array, power in each mode

    """

    nModes = np.shape(modes_sin)[0]

    modeCouplingSin = np.zeros(nModes, dtype='complex')
    modeCouplingCos = np.zeros(nModes, dtype='complex')
    modeCouplingIntensity = np.zeros(nModes)

    ms = modes_sin.astype('complex')
    mc = modes_cos.astype('complex')

    for idx, mode in enumerate(modes):

        # Calculate amplitude coupled into each mode by overlap integral
        # Need to do it for both sin and cos orientations 
        # and add then square to get total coupled intensity
        modeCouplingSin[idx] = np.sum(np.conj(field)*ms[idx]) 
        modeCouplingCos[idx] = np.sum(np.conj(field)*mc[idx]) 

        if mode.l == 0:
            modeCouplingSin[idx] = modeCouplingSin[idx] / np.sqrt(2)
            modeCouplingCos[idx] = modeCouplingCos[idx] / np.sqrt(2)            

        # Calculate power in mode
        modeCouplingIntensity[idx] = np.abs(
            modeCouplingSin[idx])**2 + np.abs(modeCouplingCos[idx])**2

    return modeCouplingCos, modeCouplingSin, modeCouplingIntensity


def propagate_through_fibre(modes, mode_cos, mode_sin, mode_coupling_cos, mode_coupling_sin, distance, rotations=False):
    """
    Propagates the given modes through a fibre over a specified distance, optionally applying mode rotations.

    Arguments:
        modes            : list
                           A list of mode objects, where each mode has a beta attribute.
        mode_cos         : np.array
                           The cosine components of the mode fields.
        mode_sin         : np.array
                           The sine components of the mode fields.
        mode_coupling_cos: np.array
                           The initial amplitude coupling into the cosine version of each mode.
        mode_coupling_sin: np.array
                           The initial amplitude coupling into the sine version of each mode.
        distance         : float
                           The distance over which the modes are propagated.
    Keyword Arguments:
        rotations        : bool, optional
                           If True, random power redistributions between cosine and sine components 
                           are applied during propagation. Defaults to False.

    Returns:
        tuple : (np.array, np.array, np.array, np.array)
                - end_field (np.array): The field at the end of the fibre after propagation.
                - end_intensity (np.array): The intensity at the end of the fibre after propagation.
                - end_mode_coupling_cos (np.array): The final amplitude coupling into the cosine version of each mode.
                - end_mode_coupling_sin (np.array): The final amplitude coupling into the sine version of each mode.
    """
    
    num_modes = len(mode_coupling_sin)
    end_field = np.zeros((np.shape(mode_sin)[1], np.shape(mode_sin)[2]), dtype='complex')

    end_mode_coupling_sin = np.zeros(num_modes, dtype='complex')
    end_mode_coupling_cos = np.zeros_like(end_mode_coupling_sin)

    for idx, mode in enumerate(modes):

        if rotations:
            weight = np.random.rand()
            total_power = mode_coupling_sin[idx]**2 + mode_coupling_cos[idx]**2
            power_cos = weight * total_power
            power_sin = (1 - weight) * total_power
            end_mode_coupling_sin[idx] = np.sqrt(power_sin) * np.sign(mode_coupling_sin[idx])
            end_mode_coupling_cos[idx] = np.sqrt(power_cos) * np.sign(mode_coupling_cos[idx])
        else:
            end_mode_coupling_sin[idx] = mode_coupling_sin[idx]
            end_mode_coupling_cos[idx] = mode_coupling_cos[idx]

        # Propagate phase
        end_mode_coupling_sin[idx] *= np.exp(1j * mode.beta * distance)
        end_mode_coupling_cos[idx] *= np.exp(1j * mode.beta * distance)

        # Calculate the coherent sum of modes to give the field at the end of the fibre
        end_field += end_mode_coupling_sin[idx] * mode_sin[idx] + end_mode_coupling_cos[idx] * mode_cos[idx]

    # Intensity is the square of the field
    end_intensity = np.abs(end_field)**2

    return end_field, end_intensity, end_mode_coupling_cos, end_mode_coupling_sin


def mode_coupling_plot(modes, mode_coupling_cos, mode_coupling_sin):
    """
    Generates a 2D numpy array to show amplitude in each mode. Array is
    organied by l and m.
    
    Arguments:
        modes             : list of instances of Mode
                            as created by find_modes
        mode_coupling_cos : np.array
                            1D numpy array giving amplitude in the cosine
                            component of each mode 
        mode_coupling_sin : np.array
                            1D numpy array giving amplitude in the sine
                            component of each mode                     
                        
    Returns:
        mode_coupling     : np.array
                            2D array (l,m) with ampitude in each mode                
    
    """

    # Need to know largest values of l and m to choose size of array
    max_l = -1
    max_m = -1
    for mode in modes:
        max_l = max(max_l, mode.l)
        max_m = max(max_m, mode.m)

    out_im = np.zeros((max_l + 1, max_m))

    for mode, coupling_cos, coupling_sin in zip(modes, mode_coupling_cos, mode_coupling_sin):
        if mode.l > 0: 
            out_im[mode.l, mode.m - 1] = np.abs(coupling_cos)**2 + np.abs(coupling_sin)**2
        else:
            out_im[mode.l, mode.m - 1] = np.abs(coupling_cos)**2 

            

    return out_im


def vis_all_modes(modes, intensity = False, rotated = False):
    """ Displays matplotlib plot of all modes in a single figure, arranged 
    by l and m.
    
    Arguments:
        modes     : list of instances of Mode
                    as created by find_modes
       
    Optional Keyword Arguments:
        intensity : boolean
                    If True, intensity of mode will be plotted, otherwise
                    default is to plot complex amplitude
        rotated   : boolean
                    If True, the sine version will be plotted, otherwise
                    default it to plot cosine
    
    Returns:
        fig       : handle to figure
    """

    # Need to know largest values of m and l so size of plot is known
    max_l = -1
    max_m = -1
    for mode in modes:
        max_l = max(max_l, mode.l)
        max_m = max(max_m, mode.m)

    if intensity:
        mode_cos, mode_sin = plot_modes_intensity(modes, 50, mode.core_radius * 2)
    else:
        mode_cos, mode_sin = plot_modes_amplitude(modes, 50, mode.core_radius * 2)

    if rotated:
        m = mode_sin
    else:
        m = mode_cos
    
    fig, axs = plt.subplots(max_l + 1, max_m, figsize=(2 * max_m, 2 * max_l))
    fig.tight_layout()
    for idx, mode in enumerate(modes):
        maxVal = np.max(np.abs(m[idx]))
        if intensity:
            axs[mode.l, mode.m - 1].imshow(mode_cos[idx], cmap='gray_r')
        else:                         
            
            axs[mode.l, mode.m - 1].imshow(mode_cos[idx],
                                       cmap=ampcol(), vmin=- maxVal, vmax=maxVal)
        axs[mode.l, mode.m - 1].title.set_text(f"LP {mode.l}, {mode.m}")

        axs[mode.l, mode.m - 1].axis('off')

    for ii in range(max_l + 1):
        for jj in range(max_m + 1):
            axs[ii, jj - 1].axis('off')
            axs[ii, jj - 1].set_aspect('equal')

    return fig


def power_in_core(mode, grid_size = 100, max_plot_radius = None):
    """ Calculates fraction of mode power that lies within core radius. The
    results will be more accurate if a larger grid_size is used. The 
    max_plot_radius should be significantly larger than the core radius (e.g. 
    1.5 times as large).
    
    Arguments:
        mode            : instance of Mode
                          as created by find_modes
    
    Keyword Arguments:
        grid_size       : int
                          size of plot in pixels (square). Defaults to 100
        max_plot_radius : float, optional
                          radial distance at edge of plot, microns. Defaults
                          to 1.5 times the core radius.                          
    Returns:        
        float           : fraction of power in core
    """
        
    if max_plot_radius is None:
        max_plot_radius = mode.core_radius * 1.5              
    
    # We integrate intensity within core and devide by integral over whole 
    # plot to get fraction in core
    im = np.abs(mode.plot_amplitude(grid_size, max_plot_radius))**2
    
    xMesh, yMesh = np.meshgrid(range(grid_size), range(grid_size))
    
    c = grid_size / 2
    rad = np.sqrt((xMesh - c)**2 + (yMesh - c)**2)
    rad = rad / (grid_size / 2) * max_plot_radius
    
    power_in_core = np.sum(im[rad < mode.core_radius]) / np.sum(im)

    return power_in_core
    


def modes_to_csv(modes, filename, header = True):
    """
    Exports the properties of a list of modes to a CSV file. Columns are
    "u", "w", "l", "m", "beta", "n_eff", "core_radius", "wavelength", "n_core", "n_cladding"

    Arguments:
        modes    : list
                   A list of instances of Mode class, as returned by find_modes
        filename : str
                   The path to the CSV file to which the mode data will be written.
    
    Keyword Arguments:
        header   : bool, optional
                   If True, writes a header row in the CSV file. Defaults to True.
    """
    
    
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)  
        
        if header:
            header = ["u", "w", "l", "m", "beta", "n_eff", "core_radius", "wavelength", "n_core", "n_cladding"]
            writer.writerow(header)
        
        for mode in modes:
            writer.writerow([
                mode.u,
                mode.w,
                mode.l,
                mode.m,
                mode.beta,
                mode.n_eff,
                mode.core_radius,
                mode.wavelength,
                mode.n_core,
                mode.n_cladding
        ])
            
    
def modes_from_csv(filename, header = True):
    """
    Loads the properties of modes from a CSV file and returns a list of mode objects.

    Arguments:
        filename   :  str
                     The path to the CSV file from which the mode data will be read.
      
    Keyword Arguments:
        header     : bool, optional
                     If True, skips the header row in the CSV file. Defaults to True.

    Returns:
        list       : A list of instances of Mode class, each corresponding to 
                     a row in the CSV file.
    """

    modes = []
    
    with open(filename, mode='r', newline='') as file:
        reader = csv.reader(file)
        
        if header:
            next(reader)  # Skip the header row
        
        for row in reader:
            mode = Mode(
                u = float(row[0]),
                w = float(row[1]),
                l = int(row[2]),
                m = int(row[3]),
                beta = float(row[4]),
                n_eff = float(row[5]),
                core_radius = float(row[6]),
                wavelength = float(row[7]),
                n_core = float(row[8]),
                n_cladding = float(row[9])
            )
            modes.append(mode)
    
    return modes    
    
