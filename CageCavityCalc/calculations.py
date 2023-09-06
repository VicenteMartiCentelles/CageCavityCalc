import numpy as np
from rdkit import Chem

from scipy.spatial import distance_matrix
from CageCavityCalc.log import logger


def sum_grid_volume(positions,radius=1.2, volume_grid_size=0.2, max_memory=1e9):
    """
    Calculate the volume of a cavity
    positions: position
    radius:
    volume_grid_size:
    max_memory: # we use 1000 MB as maximal size (it is still a lot)
    return (float): volume of the cavity
    """
    # Find min, max for the grid box, and create x,y,z of grid
    x_space = np.arange(np.min(positions.T[0])-radius, np.max(positions.T[0])+radius, volume_grid_size)
    y_space = np.arange(np.min(positions.T[1])-radius, np.max(positions.T[1])+radius, volume_grid_size)
    z_space = np.arange(np.min(positions.T[2])-radius, np.max(positions.T[2])+radius, volume_grid_size)
    # convert x,y,z axis to 3d vector
    grid_points = np.vstack(np.meshgrid(x_space,y_space,z_space)).reshape(3,-1).T
    # calcualte distance between the points

    logger.info(f"Memory for fast volume calculations ~{len(positions)*len(grid_points)*8/1e9:.1f} GB")

    # if grid point within radius of position then we add volume of the grid
    if len(positions)*len(grid_points)*8 < max_memory:
        logger.info(f"Fast summing of the grid")
        # This is twice as fast, but it requires more memory. For large cages is not an option
        dist_matrix = distance_matrix(positions, grid_points)
        volume = np.sum(np.sum(dist_matrix<radius,axis=0)>0)*(volume_grid_size**3)
    else:
        logger.info(f"Slightly slower summing of the grid")
        # Slower, but does not create large matrix
        volume = 0
        for grid_point in grid_points:
            dist_matrix = distance_matrix([grid_point], positions)
            volume += 1*(np.sum(dist_matrix < radius) > 0)
        volume*=(volume_grid_size ** 3)

    return volume
    
