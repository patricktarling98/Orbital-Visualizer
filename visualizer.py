import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
import time

WIDTH = 600
HEIGHT = 600
percent_of_points = 0.9
gridLen = 100
rangeLim = 10
numPoints = 10000

#General plot settings
pv.set_plot_theme('dark')
pv.global_theme.window_size = [WIDTH, HEIGHT]
pv.global_theme.axes.show = False
pv.global_theme.anti_aliasing = 'fxaa'
pv.global_theme.show_scalar_bar = False


"""Render image"""
def render(orbital, numPoints = 100):
    start_time = time.time()
    plotter = pv.Plotter()

    electrons = orbital.generate_electron_cloud(num_points=numPoints) #REAL WORK

    electronsPV = pv.PolyData(electrons)
    color = compute_density(electrons) #Compute density for coloring
    color_scaled = np.log(color + 1e-9)
    color_scaled = (color_scaled - color_scaled.min()) / (color_scaled.max() - color_scaled.min())
    plotter.add_mesh(electronsPV, scalars=color_scaled, point_size=5, render_points_as_spheres=True, opacity=0.5, cmap='magma')
    nucleous = pv.Sphere(radius=rangeLim/20, center = (0,0,0))
    plotter.add_mesh(nucleous, color='blue', opacity=0.9)

    #Add density coloring
    end_time = time.time()
    print(f"Generated {numPoints} points with gridsize {gridLen} and computed density in {end_time - start_time:.2f} seconds.")
    plotter.show()

"""Deals with assigning points scalar color values"""
def compute_density(points, k=100):
    tree = cKDTree(points)

    distances, _ = tree.query(points, k=k)
    furthest_distance = distances[:,-1]
    volume = (4/3) * np.pi * furthest_distance**3
    density = k / volume

    return density