"""
Computational domain related classes
Created on Sat Jun 20 19:45:03 2020

@author: mirsandiharyo
"""

import numpy as np
import math

class Domain:
    def __init__(self, lx, ly, nx, ny, gravx, gravy):
        """ 
        Initialize the domain parameters.
        """
        self.lx = lx
        self.ly = ly
        self.nx = nx
        self.ny = ny
        self.gravx = gravx
        self.gravy = gravy
        self.dx = self.lx/self.nx;
        self.dy = self.ly/self.ny;

    def get_cell_index(self, x, y, axis):
        """ 
        Fetch the indices of the eulerian cell located on the left of a 
        given point.
        """
        if (axis == 1):      # x-dir
            index_x = math.floor(x/self.dx);
            index_y = math.floor((y+0.5*self.dy)/self.dy)      
        else:                # y-dir
            index_x = math.floor((x+0.5*self.dx)/self.dx)
            index_y = math.floor(y/self.dy)
        return index_x, index_y
        
    def get_weight_coeff(self, x, y, index_x, index_y, axis):
        """ 
        Calculate the weight coefficients of a point with respect to its 
        location inside the eulerian cell.
        """
        if (axis == 1):      # x-dir
            coeff_x = x/self.dx-index_x
            coeff_y = (y+0.5*self.dy)/self.dy-index_y
        else:                # y-dir
            coeff_x = (x+0.5*self.dx)/self.dx-index_x
            coeff_y = y/self.dy-index_y
        return coeff_x, coeff_y

class Face:
    def __init__(self, domain):
        """ 
        Initialize variables (liquid is at rest at the beginning).
        """
        # velocity in x-direction
        self.u = np.zeros((domain.nx+1, domain.ny+2))
        self.u_old = np.zeros((domain.nx+1, domain.ny+2))
        self.u_temp = np.zeros((domain.nx+1, domain.ny+2))
        # velocity in y-direction
        self.v = np.zeros((domain.nx+2, domain.ny+1))
        self.v_old = np.zeros((domain.nx+2, domain.ny+1))
        self.v_temp = np.zeros((domain.nx+2, domain.ny+1))
        # forces
        self.force_x = np.zeros((domain.nx+2, domain.ny+2))
        self.force_y = np.zeros((domain.nx+2, domain.ny+2))
    
    def initialize_force(self, domain):
        """ 
        Set the forces to zero.
        """
        self.force_x = np.zeros((domain.nx+2, domain.ny+2))
        self.force_y = np.zeros((domain.nx+2, domain.ny+2))
        
    def store_old_variables(self):
        """ 
        Store old variables for second order scheme.
        """
        self.u_old = self.u.copy()
        self.v_old = self.v.copy()
        
    def store_2nd_order_variables(self):
        """ 
        Store second order variables.
        """
        self.u = 0.5*(self.u+self.u_old)
        self.v = 0.5*(self.v+self.v_old)
        
class Center:
    def __init__(self, domain):
        """
        Initialize variables stored at cell center.
        """
        # set the grid
        self.x = np.linspace(-0.5, domain.nx+2-1.5, domain.nx+2)*domain.dx
        self.y = np.linspace(-0.5, domain.ny+2-1.5, domain.ny+2)*domain.dy;
        # pressure
        self.pres = np.zeros((domain.nx+2, domain.ny+2))  
    
  