import numpy as np
import os
from gaussian_beam import GaussianBeam
from fit_gaussian_beam_data import *

# please ensure to enforce the CORRECT sign convention for your setup otherwise the results may be misleading
# for example, instead of using the convention I used in this file, if I used myBeam.set_waist_loc() to 
# define the waist as the negative of the value from the fit, I would be able to use negative distances in reconstructing
# the beam and positive distances in propagating it. Both ways should lead to identical results.

delta = -0.026908155250368385 - 2.8722247547626617e-06 # because we cannot measure where the collimator lens is, we must adjust for this fact be iterating until the beam waist loc of the reconstructed beam is zero
d_col = 3.5*10**-2 + delta # distance from fiberoptic pen beam to end of collimator
d_lens_col = 29.6*10**-2 - delta # distance from end of collimator to lens
d_lens = 0 # distance from lens to measurement mount
f_col = 8*10**-3
f_lens = 200*10**-3
d_err = 0.0005

laser_wavelength = 397*10**-9 # update me
filepath = os.path.expanduser("~/Desktop/REU/Data/397--data_20cm.csv") # update me

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
# myBeam.free_space(-d_col, d_err)
# myBeam.lens(-f_col)
# print("\n \nCollimated Beam Parameters:")
# myBeam.print()
# myBeam.free_space(-d_lens_col, d_err)
# myBeam.lens(-f_lens)
# #myBeam.free_space(d_lens, d_err)
# print("\n \nFitted Beam Parameters:")
# myBeam.print()


f_col_new  = 25.2e-3     # AL3026-A, 26 mm is about ~25.2 mm at 397 nm bc/ ThorLabs tests at like 500 nm
f_lens_new = 194.7e-3    # LB1945-A, 200 mm is aout ~194.7 mm at 397 nm bc/ ThorLabs tests at like 500 nm
d_col_new = f_col_new      # source waist to collimator = one focal length 
d_lens_col_new = 10.5*0.0254   # Distance bebween lens and columnator

myBeam.free_space(-d_col_new, d_err)
myBeam.lens(-f_col_new)
print("\n\nCollimated Beam (proposed):")
myBeam.print()
myBeam.free_space(-d_lens_col_new, d_err)
myBeam.lens(-f_lens_new)
myBeam.free_space(d_lens, d_err)
print("\n\nFinal Beam (proposed):")
myBeam.print()