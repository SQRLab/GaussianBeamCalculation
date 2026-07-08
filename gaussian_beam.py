import numpy as np

# define the gaussian beam class
# a Gaussian Beam object contains a laser wavelength, waist location, rayleigh range, beam quality factor, and errors for all but the laser wavelength
class GaussianBeam:
    def __init__(self, waist, waist_loc, m2, waist_err, waist_loc_err, m2_err, laser_wavelength):
        """ This class is used to describe a Gaussian beam. 
            
            Parameters:
                waist: the smallest radius (or waist) of the beam
                waist_loc: location of the beam waist
                m2: beam quality factor
                waist_err: beam waist size uncertainity
                waist_loc_err: beam waist location uncertainty
                m2_err: beam quality factor uncertainty
                laser_wavelength: the wavelength of the given laser

        """
        self.waist_loc = waist_loc
        self.m2 = m2
        self.wavelength = laser_wavelength

        self.waist_loc_err = waist_loc_err
        self.m2_err = m2_err

        self.rayleigh_range = np.pi * (waist ** 2) / (self.m2) / (self.wavelength)
        self.rayleigh_range_err = np.pi/(self.wavelength)*np.sqrt((2*waist*waist_err/m2)**2+(waist**2*m2_err/m2**(3/2))**2)

    ## GETTERS ##    
    def get_waist_loc(self):
        return self.waist_loc

    def get_waist_loc_err(self):
        return self.waist_loc_err

    def get_rayleigh_range(self):
        return self.rayleigh_range

    def get_rayleigh_range_err(self):
        return self.rayleigh_range_err

    def get_m2(self):
        return self.m2

    def get_m2_err(self):
        return self.m2_err
    
    def get_waist(self):
        # calculates the waist size
        return np.sqrt(self.rayleigh_range*self.wavelength*self.m2/np.pi)
    
    def get_waist_err(self):
        # calculates the uncertainty in waist size
        return np.sqrt(self.wavelength/np.pi)*np.sqrt((0.5*self.rayleigh_range_err*np.sqrt(self.m2/self.rayleigh_range))**2+(0.5*self.m2_err*np.sqrt(self.rayleigh_range/self.m2))**2)


    ## SETTERS ##    
    def set_waist_loc(self, loc):
        self.waist_loc = loc

    def set_waist_loc_err(self, loc_err):
        self.waist_loc_err = loc_err
    
    ## METHODS ##

    def free_space(self, d, d_err):
        # propogates beam through free space a distance d where distance is known within d_err
        # for backwards propagation, use -d
        self.waist_loc += d
        self.waist_loc_err = np.sqrt((self.waist_loc_err)**2+(d_err)**2)

    def lens(self, f):
        # propogates beam through a lens with focal length f
        # for backwards propagation, use -f
        q = self.waist_loc + 1j*self.rayleigh_range
        q_err = self.waist_loc_err + 1j*self.rayleigh_range_err
        q_update = q/(-1/f*q+1)
        q_update_err = q_err/(-1/f*q+1)**2

        self.waist_loc = np.real(q_update)
        self.waist_loc_err = np.real(q_update_err)
        self.rayleigh_range = np.imag(q_update)
        self.rayleigh_range_err = np.imag(q_update_err)

    def print(self):
        print("Beam waist", self.get_waist(), "±", self.get_waist_err())
        print("Beam waist loc", self.get_waist_loc(), "±", self.get_waist_loc_err())
        print("Rayleigh Range", self.get_rayleigh_range(), "±", self.get_rayleigh_range_err())
        print("M^2", self.get_m2(), "±", self.get_m2_err())


