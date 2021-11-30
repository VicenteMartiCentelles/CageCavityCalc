import numpy as np
import math

class GridPoint:
    def __init__(self, i, j, k, points, center, grid_spacing):
        self.i = i
        self.j = j
        self.k = k     
        
        self.i_from_center = i - points
        self.j_from_center = j - points
        self.k_from_center = k - points
        
        self.x = self.i_from_center * grid_spacing + center[0]
        self.y = self.j_from_center * grid_spacing + center[1]
        self.z = self.k_from_center * grid_spacing + center[2]
        
        self.pos = [self.x, self.y, self.z]
        
        self.d_from_center = np.linalg.norm( np.array(self.pos) -  np.array(center) )
        
        self.inside_cavity = 0
        
        self.overlapping_with_cage = 0
        
        self.number_of_neighbors = 0
        
        self.is_window = 0
        
        self.vector_angle = 0
        
        self.neighbors = []

class CageGrid:
    def __init__(self, center, radius, delta, grid_spacing):
        self.center = center
        self.size = radius + delta
        self.grid_spacing = grid_spacing
        self.points = math.ceil(self.size / self.grid_spacing)
        self.n = 2*self.points + 1  
        self.grid = [GridPoint(i,j,k, self.points, self.center, self.grid_spacing) for k in range(self.n) for j in range(self.n) for i in range(self.n)]
        self.grid = np.array(self.grid)
        self.gridPosList = []
        for i in self.grid:
            self.gridPosList.append(i.pos)

