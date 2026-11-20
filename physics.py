import numpy as np
bohr_radius = 1 # in meters, but we can use it as a unit of distance for our model
gridLen = 100
rangeLim = 10 # limit for the grid in each direction, can be adjusted to capture more of the electron cloud

"""OLD FILE FOR PHYSICS, THESE ARE ALL UNDER ATOM/ORBITAL CLASSES NOW EXCEPT THE SAMPLE POINTS FUNCTION"""

#Bread n Butter
"------------------------------------------------------------------------------------------------------------------------------------------------------------- "
"""Generates a cloud of points"""
def generate_electron_cloud(numPoints = 100, gridSize = (gridLen, gridLen, gridLen)):
    xv, yv, zv = np.linspace(-rangeLim, rangeLim, gridSize[0]), np.linspace(-rangeLim, rangeLim, gridSize[1]), np.linspace(-rangeLim, rangeLim, gridSize[2])
    x, y, z = np.meshgrid(xv, yv, zv, indexing='ij')

    pdf = calculate_probability_distribution(x, y, z, numPoints=numPoints)#Get pdf, assigne each point a probability value
    points = sample_points_from_pdf(pdf, x, y, z, numPoints=numPoints)#Sample points from the pdf to create the electron cloud

    return points


#Boring Helper Methods
"--------------------------------------------------------------------------------------------------------------------------------------------------------------"
"""Calculates the pdf, and populates an array with the pdf values for all points in space"""
def calculate_probability_distribution(x, y, z,numPoints = 100):
    r = np.sqrt(x**2 + y**2 + z**2)#Find r for all points
    theta = np.zeros_like(r)
    nonzero_mask = r > 0
    theta[nonzero_mask] = np.arccos(z[nonzero_mask]/r[nonzero_mask])#find theta for all points, avoid division by zero
    phi = np.arctan2(y, x)#find phi for all points

    #TODO:Need to downscale r before calculating pdf - radial distribution?
    downscale = r.max()/(5*bohr_radius) #Scale down distances such that rmax is bohr radius
    r_scaled = r / downscale

    #TODO: Need to standardize this for all atoms, currently only works for hydrogen
    pdf = 1/(np.pi*bohr_radius**3) * np.exp(-2*r_scaled/bohr_radius)

    pdf/= pdf.sum()
    return pdf

"""Samples points from a given pdf, returns an array of points"""
def sample_points_from_pdf(pdf, x, y, z, numPoints):
    pdfFlat = pdf.ravel()
    xFlat = x.ravel()
    yFlat = y.ravel()
    zFlat = z.ravel()

    rng = np.random.default_rng(seed=0)
    indexs = rng.choice(len(pdfFlat), size=numPoints, p=pdfFlat)
    points = np.column_stack((xFlat[indexs], yFlat[indexs], zFlat[indexs]))

    return points

def find_scale_factor(points):
    pass

