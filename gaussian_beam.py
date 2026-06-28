import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import erf 

def modified_erf(x, power_min, power_max, beam_radius, beam_center):
    # function describing beam intensity
    return power_min + power_max/2 * (1-erf(np.sqrt(2) * (x - beam_center)/beam_radius))

def gaussian_beam_fn(z, beam_radius, beam_waist_loc, m):
    # gaussian beam equation
    return np.sqrt((beam_radius ** 2) * (1 + ((z - beam_waist_loc) ** 2) * (((m ** 2) * (397e-6) / np.pi / (beam_radius ** 2)) ** 2)))

def extract_beam_parameters(z_list, x_list, power_list):
    # given a list of z values as well as a list of lists for x values and powers at each z,
    # calculates and prints the beam radius, beam waist location, and M for a gaussian beam.
    beam_radius_list = list()

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

        # Compare the data and the fit
        plt.plot(x_vals, power_vals, 'o', label='Data')
        plt.plot(x_vals, modified_erf(x_vals, power_min, power_max , beam_radius, beam_center), '-', label='Fit')
        plt.title(f"Knife x Position vs. Power for z={z} cm")
        plt.xlabel("Knife x Position (mm)")
        plt.ylabel("Power (mW)")
        plt.legend()
        plt.show()

        displacement = np.arange(np.mean(x_vals) - 5 * beam_radius, np.mean(x_vals) + 5 * beam_radius, 0.01) - beam_center
        intensity = 2 * max(power_vals) / np.pi / (beam_radius ** 2) * np.exp(-2 * ((displacement) / beam_radius) ** 2)
        plt.plot(displacement, intensity)
        plt.title("Displacement vs. Intensity")
        plt.xlabel("Displacement (mm)")
        plt.ylabel("Intensity (mW/mm^2)")
        plt.show()

        beam_radius_list.append(beam_radius)

    print("Beam Radius List:", beam_radius_list)

    # Compare the data and the fit
    init_guess = [np.mean(beam_radius_list), np.mean(z_list), 1]
    popt, pcov = opt.curve_fit(gaussian_beam_fn, z_list, beam_radius_list, p0=init_guess)

    # Extract the fitted parameters
    beam_radius, beam_waist_loc, m = popt
    print("Beam Radius:", beam_radius, "mm")
    print("Beam Waist Location:", beam_waist_loc, "mm")
    print("M^2:", m ** 2)
    print("Rayleigh Range:", np.pi * (beam_radius ** 2) / (m ** 2) / (397e-6), "mm")
    print("Divergence Angle:", (m ** 2) * (397e-6) / np.pi / beam_radius * 360 / (2 * np.pi), "degrees")

    # Plot the data and fit
    plt.plot(z_list, beam_radius_list, 'o', label='Data')
    z = np.arange(np.min(z_list), np.max(z_list), 1)
    plt.plot(z, gaussian_beam_fn(z, beam_radius, beam_waist_loc, m), '-', label='Fit')
    plt.title("z-coordinate vs. Beam Radius")
    plt.xlabel("z-coordinate (mm)")
    plt.ylabel("Beam Radius (mm)")
    plt.show()
    return beam_radius, beam_waist_loc, m


if __name__ == '__main__':
    filepath = "C:/Users/rlhaa/Desktop/School/UW REU/lens_data_06262026.csv"
    # read in data using the data frame
    df = pd.DataFrame(pd.read_csv(f"{filepath}"))
    # compress table so each entry corresponds to list of values for a given z
    df = df.groupby("z (cm)")[["x (cm)", "Intensity (microW)"]].agg(list).reset_index()
    z_list = df["z (cm)"].apply(np.array).to_numpy()
    print(z_list)
    x_list = df["x (cm)"].apply(np.array).to_numpy()
    print(x_list)
    power_list = df["Intensity (microW)"].apply(np.array).to_numpy()
    print(power_list)
    
    beam_radius, beam_waist_loc, m = extract_beam_parameters(z_list, x_list, power_list)
    