import numpy as np
from scipy.special import genlaguerre, sph_harm_y, factorial
from scipy.integrate import quad
import visualizer

gridLen = visualizer.gridLen
percent_of_points = visualizer.percent_of_points

class Orbital:
    def __init__(self, n, l, m, Z=1):
        self.n = n
        self.l = l
        self.m = m
        self.Z = Z
        self.rangeLim = self.find_cutoff_radius(percentage=percent_of_points)  # Set rangeLim based on the cutoff radius

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
        return sph_harm_y(self.l, self.m, theta, phi)

    def wave_function(self, atom, r, theta, phi):
        """Calculate the wave function at all points; atom is unused."""
        # Calculate wave function for this orbital
        R = self.radial_component(r)
        Y = self.angular_component(theta, phi)
        print("Number of NaN values in R:", np.isnan(R).sum())  # Debugging line to check for NaN values
        print("Number of NaN values in Y:", np.isnan(Y).sum())  # Debugging line to check for NaN values
        print("R sum:", np.sum(R))  # Debugging line to check the sum of R values
        print("Y sum:", np.sum(Y))  # Debugging line to check the sum of Y values
        return R * Y    
    
    def probability_density(self, atom, r, theta, phi):
        psi = self.wave_function(atom, r, theta, phi) #REAL PHYSICS

        pdf = np.abs(psi)**2
        pdf[r > self.rangeLim] = 0  # Set probability density to zero outside the cutoff radius
        print("Number of NaN values in pdf:", np.isnan(pdf).sum())  # Debugging line to check for NaN values
        print("pdf sum:", np.sum(pdf))  # Debugging line to check the sum of pdf values
        return pdf    
    

    "-------------------------------------------------BREAD N BUTTER---------------------------------------------------------------"
    def generate_electron_cloud(self, num_points=100):
        #MAKE GRID N ASSIGN EACH POINT R,THETA,PHI VALUES
        xv, yv, zv = np.linspace(-self.rangeLim, self.rangeLim, gridLen), np.linspace(-self.rangeLim, self.rangeLim, gridLen), np.linspace(-self.rangeLim, self.rangeLim, gridLen)
        x, y, z = np.meshgrid(xv, yv, zv, indexing='ij')

        r = np.sqrt(x**2 + y**2 + z**2)

        #just prevents division by zero
        theta = np.zeros_like(r)
        nonzero_mask = r > 0
        cos_theta = z[nonzero_mask] / r[nonzero_mask]
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta[nonzero_mask] = np.arccos(cos_theta)  

        phi = np.arctan2(y, x)

        orbital_pdf = self.probability_density(self, r, theta, phi) 
        orbital_pdf /= np.sum(orbital_pdf)  # Normalize the PDF to ensure it sums to 1

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
    

    def find_cutoff_radius(self, percentage=0.9):
        """Finds the radius that contains a given percentage of the electron cloud, with units of bhor radius"""
        #start at this radius, and increase till desired r
        r_cutoff = 4*(self.n**2)/self.Z

        net_prob, _ = quad(self.radial_pdf, 0, r_cutoff)

        while(net_prob < percentage):
            r_cutoff *= 2
            net_prob, _ = quad(self.radial_pdf, 0, r_cutoff)

        return r_cutoff

    def radial_pdf(self, r):
        """Calculate the radial probability density function for this orbital."""
        R = self.radial_component(r)
        pdf = (r**2) * np.abs(R)**2
        return pdf

    
    