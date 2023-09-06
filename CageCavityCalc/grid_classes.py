import numpy as np
import math

class GridPoint:
    def __init__(self, i, j, k, points, center, grid_spacing):
        if isinstance(points, float): ## Not used, in this previous version we used a square box
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
        
        if isinstance(points, list):
            self.i = i
            self.j = j
            self.k = k     
            
            self.i_from_center = i - ((points[0]-1)/2)
            self.j_from_center = j - ((points[1]-1)/2)
            self.k_from_center = k - ((points[2]-1)/2)
            
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
    def __init__(self, center, box_size, delta, grid_spacing):
        
        if isinstance(box_size, float): ## Not used, in this previous version we used a square box
            self.center = center
            self.size = box_size + delta
            self.grid_spacing = grid_spacing
            self.points = math.ceil(self.size / self.grid_spacing)
            self.n = 2*self.points + 1
            self.grid = [GridPoint(i,j,k, self.points, self.center, self.grid_spacing) for k in range(self.n) for j in range(self.n) for i in range(self.n)]
            self.grid = np.array(self.grid)
            self.gridPosList = []
            for i in self.grid:
                self.gridPosList.append(i.pos)
                
        if isinstance(box_size, list):
            self.center = center
            self.size_x_min = box_size[0]
            self.size_x_max = box_size[1]
            self.size_y_min = box_size[2]
            self.size_y_max = box_size[3]
            self.size_z_min = box_size[4]
            self.size_z_max = box_size[5]
            self.grid_spacing = grid_spacing
            self.points_x = math.ceil((self.size_x_max-self.size_x_min + delta) / self.grid_spacing)+1
            self.points_y = math.ceil((self.size_y_max-self.size_y_min + delta) / self.grid_spacing)+1
            self.points_z = math.ceil((self.size_z_max-self.size_z_min + delta) / self.grid_spacing)+1




            self.points = [self.points_x,self.points_y,self.points_z]
            self.grid = [GridPoint(i,j,k, self.points, self.center, self.grid_spacing) for k in range(self.points_z) for j in range(self.points_y) for i in range(self.points_x)]
            self.grid = np.array(self.grid)
            self.gridPosList = []
            for i in self.grid:
                self.gridPosList.append(i.pos)

