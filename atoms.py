import numpy as np
import physics
import json
from scipy.special import genlaguerre, sph_harm_y, factorial

gridLen = 100
rangeLim = 10
percent_of_points = 0.9


"""Atom class: Has a name, atomic number and a list of all orbital objects which contain orbital data"""
class Atom:
    def __init__(self, name):
        self.name = name
        self.Z = self.load_atomic_number()
        self.orbitals = self.generate_orbitals()
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

        #Gets pdf of each orbital, adds them up, and normalizes AFTER all have been summed
        for orbital in self.orbitals:
            orbital_pdf = orbital.probability_density(self, r, theta, phi) 
            total_pdf += orbital_pdf
        total_pdf /= total_pdf.sum()

        #RETURNS TO DISPLAY AFTER DOING WEIGHTED SAMPLING
        return physics.sample_points_from_pdf(total_pdf, x, y, z, numPoints=num_points)
    "--------------------------------------------------------------LOSER HELPERS BELOW--------------------------------------------------"

    def load_atomic_number(self): #Load atomic number "Z"
        atom_data = json.load(open('orbitals.json'))
        atom = atom_data[self.name]
        atomic_number = atom['Z']
        return atomic_number

    def generate_orbitals(self): #Load list of orbital objects containing n,m,l
        atom_data = json.load(open('orbitals.json'))
        self_atom = atom_data[self.name]
        self_orbitals = self_atom["orbitals"]
        orbital_objects = []
        for orbital in self_orbitals:
            n = orbital['n']
            l = orbital['l']
            m = orbital['m']
            occupancy = orbital['occupancy']
            orbital_objects.append(Orbital(n, l, m, self.Z, occupancy))

        return orbital_objects
    
    """Finds the scale factor such that r*scale results in the largest r being the r containing 90% of points"""
    def find_scale_factor(self, r): #eventually change this so that it takes x%
        #Hydrogenic orbitals scale as 4*n^2/Z for 90% so we will use this for now
        n_max = max(orbital.n for orbital in self.orbitals)
        Z = self.Z
        return 4*(n_max**2)/(Z*r.max())
    
"----------------------------------------------------------------------------------------------------------------------------------------"




class Orbital:
    def __init__(self, n, l, m, Z, occupancy):
        self.n = n
        self.l = l
        self.m = m
        self.Z = Z
        self.occupancy = occupancy

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
    