import numpy as np
from scipy.special import genlaguerre, sph_harm_y, factorial
import visualizer

gridLen = visualizer.gridLen
rangeLim = visualizer.rangeLim
percent_of_points = visualizer.percent_of_points

class Orbital:
    def __init__(self, n, l, m, Z=1):
        self.n = n
        self.l = l
        self.m = m
        self.Z = Z

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
        return sph_harm_y(self.m, self.l, theta, phi)

    def wave_function(self, atom, r, theta, phi):
        """Calculate the wave function at all points; atom is unused."""
        # Calculate wave function for this orbital
        R = self.radial_component(r)
        Y = self.angular_component(theta, phi)
        return R * Y    
    
    def probability_density(self, atom, r, theta, phi):
        psi = self.wave_function(atom, r, theta, phi) #REAL PHYSICS

        pdf = np.abs(psi)**2 #TODO is this correct? Should we multiply by r^2 sin(theta) for the volume element in spherical coordinates?
        return pdf    
    

    "-------------------------------------------------BREAD N BUTTER---------------------------------------------------------------"
    def generate_electron_cloud(self, num_points=100):
        #MAKE GRID N ASSIGN EACH POINT R,THETA,PHI VALUES
        xv, yv, zv = np.linspace(-rangeLim, rangeLim, gridLen), np.linspace(-rangeLim, rangeLim, gridLen), np.linspace(-rangeLim, rangeLim, gridLen)
        x, y, z = np.meshgrid(xv, yv, zv, indexing='ij')

        r = np.sqrt(x**2 + y**2 + z**2)

        #TODO scale r
        scale_factor = self.find_scale_factor(r)
        r_scaled = r * scale_factor

        #just prevents division by zero
        theta = np.zeros_like(r)
        nonzero_mask = r > 0
        cos_theta = z[nonzero_mask] / r[nonzero_mask]
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta[nonzero_mask] = np.arccos(cos_theta)  

        phi = np.arctan2(y, x)

        orbital_pdf = self.probability_density(self, r, theta, phi) 
        orbital_pdf /= orbital_pdf.sum()  # Normalize the PDF

        return self.sample_points_from_pdf(orbital_pdf, x, y, z, num_points)
    
    """Samples points from a given pdf, returns an array of points"""
    def sample_points_from_pdf(self, pdf, x, y, z, numPoints):
        #Flatten arrays
        pdfFlat = pdf.ravel()
        xFlat = x.ravel()
        yFlat = y.ravel()
        zFlat = z.ravel()

        rng = np.random.default_rng(seed=0)
        indexs = rng.choice(len(pdfFlat), size=numPoints, p=pdfFlat) #Select indices based on the pdf probabilities
        points = np.column_stack((xFlat[indexs], yFlat[indexs], zFlat[indexs]))

        return points
    
    """Finds the scale factor such that r*scale results in the largest r being the r containing 90% of points"""
    def find_scale_factor(self, r): #eventually change this so that it takes x%
        #TODO make this correct
        Z = self.Z
        return 4*(self.n**2)/(Z*r.max())

    def find_cutoff_radius(self, r, percentage=percent_of_points):
        # Flatten the array and sort it
        r_flat = r.ravel()
        r_sorted = np.sort(r_flat)

        # Find the index corresponding to the desired percentage
        index = int(len(r_sorted) * percentage)

        # Return the cutoff radius
        return r_sorted[index]
    