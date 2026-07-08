## GaussianBeamCalculation
This repository contains files used to determine beam parameters for experimental data and propagate a gaussian beam through optical components.

## fit_gaussian_beam
- When using the file, update the filepath to the data file and the laser wavelength
    - Data files are be .csv files with columns z (cm), x (cm), and Power (mW)
- Calculates the beam radius, beam waist location, and beam quality factor for a Gaussian Beam
- Determines the reduced chi squared and probability to exceed for the Gaussian Beam fit
- Plots the beam profile at each location z and over the full range of z values

## gaussian_beam
- The GaussianBeam class
  - Allows construction of a GaussianBeam object
  - Propagate beam through optical components using provided methods
 
## propogate_beam
- An example file that demonstrates how to use key functionalities of gaussian_beam and fit_gaussian_beam
