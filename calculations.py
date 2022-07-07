import numpy as np
import MDAnalysis
from rdkit import Chem

def sum_grid_volume(positions,radius=1.2, volume_grid_size=0.2):
    """
    Calculate the volume of a cavity
    positions: position
    radius:
    volume_grid_size:
    return (float): volume of the cavity
    """
    # Find min, max for the grid box, and create x,y,z of grid
    x_space = np.arange(np.min(positions.T[0])-radius, np.max(positions.T[0])+radius, volume_grid_size)
    y_space = np.arange(np.min(positions.T[1])-radius, np.max(positions.T[1])+radius, volume_grid_size)
    z_space = np.arange(np.min(positions.T[2])-radius, np.max(positions.T[2])+radius, volume_grid_size)
    # convert x,y,z axis to 3d vector
    grid_points = np.vstack(np.meshgrid(x_space,y_space,z_space)).reshape(3,-1).T
    # calcualte distance between the points
    dist_matrix = MDAnalysis.lib.distances.distance_array(positions, grid_points)
    # if grid point within radius of position then we add volume of the grid
    return np.sum(np.sum(dist_matrix<radius,axis=0)>0)*(volume_grid_size**3)
    
