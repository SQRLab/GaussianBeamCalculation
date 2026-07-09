import numpy as np
import os
from gaussian_beam import GaussianBeam
from fit_gaussian_beam_data import *

# please ensure to enforce the CORRECT sign convention for your setup otherwise the results may be misleading
# for example, instead of using the convention I used in this file, if I used myBeam.set_waist_loc() to 
# define the waist as the negative of the value from the fit, I would be able to use negative distances in reconstructing
# the beam and positive distances in propagating it. Both ways should lead to identical results.

delta = -0.026914794220628382 # because we cannot measure where the collimator lens is, we must adjust for this fact be iterating until the beam waist loc of the reconstructed beam is zero
d_col = 3.5*10**-2 + delta # distance from fiberoptic pen beam to end of collimator
d_lens_col = 3*10**-2 - delta # distance from end of collimator to lens
d_lens = 15.9*10**-2 # distance from lens to zero point of data set
f_col = 8*10**-3
f_lens = 200*10**-3
d_err = 0

laser_wavelength = 650*10**-9 # update me
filepath = "C:/Users/rlhaa/Desktop/UW REU/650nm_data/lens_data_06302026.csv" # update me

# fit experimental data to determine beam parameters
beam_params = process_data(filepath, laser_wavelength)

# set beam using fitted beam parameters
myBeam = GaussianBeam(*beam_params, laser_wavelength)

# reconstruct original beam (use negatives because it is backwards but with my sign convention those become positive)
myBeam.free_space(d_lens, d_err)
myBeam.lens(f_lens)
myBeam.free_space(d_lens_col, d_err)
myBeam.lens(f_col)
myBeam.free_space(d_col, d_err)
print("\n \nOriginal Beam Parameters:")
myBeam.print()

# pass the beam back through to ensure you get the fitted beam back (with my sign convention these distances have to be negative)
myBeam.free_space(-d_col, d_err)
myBeam.lens(-f_col)
print("\n \nCollimated Beam Parameters:")
myBeam.print()
myBeam.free_space(-d_lens_col, d_err)
myBeam.lens(-f_lens)
#myBeam.free_space(d_lens, d_err)
print("\n \nFitted Beam Parameters:")
myBeam.print()
