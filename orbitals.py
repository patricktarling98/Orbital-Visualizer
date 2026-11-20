import numpy as np
import physics
import json
from scipy.special import genlaguerre, sph_harm_y, factorial

gridLen = 100
rangeLim = 10
percent_of_points = 0.9

class Orbital:
    def __init__(self, n, l, m, Z=1):
        self.n = n
        self.l = l
        self.m = m
        self.Z = Z

    def probability_density(self, atom, r, theta, phi):
        psi = self.wave_function(atom, r, theta, phi) #REAL PHYSICS

        pdf = np.abs(psi)**2
        return pdf

    """CALCULATE WAVE FUNCTION FOR THIS ORBITAL, returns sqrt(probability) at all points"""
    def wave_function(self, atom, r, theta, phi):
        # Calculate wave function for this orbital
        R = self.radial_component(r)
        Y = self.angular_component(theta, phi)
        return R * Y    
    "---------------------------------------------------------STR8 VIBE CODED--------------------------------------------------------------------"
    
    def radial_component(self, r):
        n, l, Z = self.n, self.l, self.Z

        rho = 2 * Z * r / n

        prefactor = np.sqrt(
            (2*Z/n)**3 *
            factorial(n-l-1) /
            (2*n * factorial(n+l))
        )

        Lag = genlaguerre(n-l-1, 2*l+1)

        return prefactor * np.exp(-rho/2) * rho**l * Lag(rho)  

    def angular_component(self, theta,phi):
        return sph_harm_y(self.m, self.l, phi, theta)

    "-------------------------------------------------BREAD N BUTTER---------------------------------------------------------------"
    def generate_electron_cloud(self, num_points=100):
        #MAKE GRID N ASSIGN EACH POINT R,THETA,PHI VALUES
        xv, yv, zv = np.linspace(-rangeLim, rangeLim, gridLen), np.linspace(-rangeLim, rangeLim, gridLen), np.linspace(-rangeLim, rangeLim, gridLen)
        x, y, z = np.meshgrid(xv, yv, zv, indexing='ij')

        r = np.sqrt(x**2 + y**2 + z**2)

        #Scale r now?
        scale_factor = self.find_scale_factor(r)
        r_scaled = r * scale_factor

        #Vibe coded tbh
        theta = np.zeros_like(r)
        nonzero_mask = r > 0
        cos_theta = z[nonzero_mask] / r[nonzero_mask]
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta[nonzero_mask] = np.arccos(cos_theta)        
        phi = np.arctan2(y, x)

        total_pdf = np.zeros_like(r)

        
        orbital_pdf = self.probability_density(self, r, theta, phi) 
        total_pdf += orbital_pdf
        total_pdf /= total_pdf.sum()

        #RETURNS TO DISPLAY AFTER DOING WEIGHTED SAMPLING
        return physics.sample_points_from_pdf(total_pdf, x, y, z, numPoints=num_points)
    
    """Finds the scale factor such that r*scale results in the largest r being the r containing 90% of points"""
    def find_scale_factor(self, r): #eventually change this so that it takes x%
        #Hydrogenic orbitals scale as 4*n^2/Z for 90% so we will use this for now
        Z = self.Z
        return 4*(self.n**2)/(Z*r.max())

    