import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import erf 
from functools import partial
from scipy.stats import chi2

def modified_erf(x, power_min, power_max, beam_radius, beam_center):
    # function describing beam intensity
    return power_min + power_max/2 * (1-erf(np.sqrt(2) * (x - beam_center)/beam_radius))

def gaussian_beam_fn(z, beam_radius, beam_waist_loc, m, laser_wavelength):
    # gaussian beam equation
    return np.sqrt((beam_radius ** 2) * (1 + ((z - beam_waist_loc) ** 2) * (((m ** 2) * (laser_wavelength) / np.pi / (beam_radius ** 2)) ** 2)))

def fit_quality(fit_func, input_vals, output_vals, output_err, fit_params):
    # given a fit function and an experimental data set input_vals, output_vals, and output_err,
    # determines the reduced chi squared and pte of the fit_params
    yfit = fit_func(input_vals, *fit_params)
    residuals = output_vals-yfit # calculate residuals;

    # calculate and plot normalized residuals
    norm_residuals = residuals/output_err

    # Calculated chi_squared, reduced chi_squared, and PTE
    chi_squared = np.sum(norm_residuals**2)

    M = len(output_vals)
    df = M-len(fit_params)
    reduced_chi_squared = chi_squared/df

    PTE = 1 - chi2.cdf(chi_squared, df=df)
    return PTE, reduced_chi_squared

def extract_beam_parameters(z_list, x_list, power_list, laser_wavelength, plot_xs=False):
    # given a list of z values as well as a list of lists for x values and powers at each z,
    # calculates and prints the beam radius, beam waist location, and M for a gaussian beam.
    beam_radius_list = list()
    beam_radius_list_err = list()

    for i in range(len(z_list)):
        # extract variables at this z value
        z = z_list[i]
        x_vals = x_list[i]
        power_vals = power_list[i]

        # Initial guesses for parameters
        init_guess = [min(power_vals), max(power_vals), 1, np.mean(x_vals)]

        # Fit the data
        popt, pcov = opt.curve_fit(modified_erf, x_vals, power_vals, p0=init_guess)

        # Extract the fitted parameters
        power_min, power_max, beam_radius, beam_center = popt
        perr = np.sqrt(np.diag(pcov))
        power_min_err, power_max_err, beam_radius_err, beam_center_err = perr

        displacement = np.arange(np.mean(x_vals) - 5 * beam_radius, np.mean(x_vals) + 5 * beam_radius, 0.01) - beam_center
        intensity = 2 * max(power_vals) / np.pi / (beam_radius ** 2) * np.exp(-2 * ((displacement) / beam_radius) ** 2)


        beam_radius_list.append(np.abs(beam_radius*10**-2)) #convert to meters
        beam_radius_list_err.append(np.abs(beam_radius_err*10**-2)) #convert to meters
        if plot_xs == True:
            #Compare the data and the fit
            plt.plot(x_vals, power_vals, 'o', label='Data')
            x_for_plot = np.linspace(np.min(x_vals), np.max(x_vals), 500)
            plt.plot(x_for_plot, modified_erf(x_for_plot, power_min, power_max , beam_radius, beam_center), '-', label='Fit')
            plt.title(f"Knife x Position vs. Power for z={z} cm")
            plt.xlabel("Knife x Position (cm)")
            plt.ylabel("Power (mW)")
            plt.legend()
            plt.show()

            plt.plot(displacement, intensity)
            plt.title("Displacement vs. Intensity")
            plt.xlabel("Displacement (m)")
            plt.ylabel("Intensity (mW/cm^2)")
            plt.show()


    #print("Beam Radius List:", beam_radius_list)
    z_list = z_list*10**-2
    #print("Z Values: ", z_list)

    # define fit function with set laser_wavelength (have to do this way bc of how curve_fit works)
    fixed_wavelength_gaussian = partial(gaussian_beam_fn, laser_wavelength=laser_wavelength)
    
    # Compare the data and the fit
    init_guess = [np.mean(beam_radius_list), np.mean(z_list), 1]
    popt, pcov = opt.curve_fit(fixed_wavelength_gaussian, z_list, beam_radius_list, p0=init_guess)

    # Extract the fitted parameters
    beam_waist, beam_waist_loc, m = popt
    perr = np.sqrt(np.diag(pcov))
    beam_waist_err, beam_waist_loc_err, m_err = perr
    print("Fitted Beam Parameters:")
    print("Beam Waist:", beam_waist, "±", beam_waist_err, "m")
    print("Beam Waist Location:", beam_waist_loc, "±", beam_waist_loc_err, "m")
    print("M^2:", m ** 2, "±", 2*m*m_err)
    print("Rayleigh Range:", np.pi * (beam_waist ** 2) / (m ** 2) / (laser_wavelength), "±", np.pi/(laser_wavelength)*np.sqrt((2*beam_waist*beam_waist_err/m**2)**2+(beam_waist**2*2*m_err/m**2)**2),"m")
    print("Divergence Angle:", (m ** 2) * (laser_wavelength) / np.pi / beam_radius * 360 / (2 * np.pi), "±", 180/np.pi**2* (laser_wavelength)*np.sqrt((2*m*m_err/beam_waist)**2+(m**2/beam_waist**2*beam_waist_err)**2),"degrees")
    
    # Plot the data and fit
    plt.errorbar(z_list, beam_radius_list, beam_radius_list_err, fmt='o', label='Data')
    z = np.linspace(np.min(z_list), np.max(z_list), 100)
    plt.plot(z, fixed_wavelength_gaussian(z, *popt), '-', label='Fit')
    plt.title("z-coordinate vs. Beam Radius")
    plt.xlabel("z-coordinate (m)")
    plt.ylabel("Beam Radius (m)")
    plt.legend()
    plt.show()
    plt.savefig("Beam Fitted")

    # determine quality of fit
    PTE, reduced_chi2 = fit_quality(fixed_wavelength_gaussian, z_list, beam_radius_list, beam_radius_list_err, popt)
    print("PTE", PTE)
    print("Reduced Chi Squared", reduced_chi2)

    return beam_waist, beam_waist_loc, m**2, beam_waist_err, beam_waist_loc_err, 2*m*m_err


def process_data(filepath, laser_wavelength, plot_xs=False):
    # read in data using the data frame
    df = pd.DataFrame(pd.read_csv(f"{filepath}"))

    # compress table so each entry corresponds to list of values for a given z
    df = df.groupby("z (cm)")[["x (cm)", "Intensity (mW)"]].agg(list).reset_index()
    
    z_list = df["z (cm)"].apply(np.array).to_numpy()
    
    x_list = df["x (cm)"].apply(np.array).to_numpy()
    
    power_list = df["Intensity (mW)"].apply(np.array).to_numpy()
    
    return extract_beam_parameters(z_list, x_list, power_list, laser_wavelength, plot_xs=False)
    

if __name__ == '__main__':
    laser_wavelength = 650*10**-9 # update me
    filepath = "C:/Users/rlhaa/Desktop/School/UW REU/lens_data_06302026.csv" # update me

    process_data(filepath, laser_wavelength, plot_xs=False)