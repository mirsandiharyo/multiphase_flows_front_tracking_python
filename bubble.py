"""
Bubble class
Created on Sat Jun 20 19:06:17 2020

@author: mirsandiharyo
"""

import numpy as np
import math

class Bubble:
    total = 0
    
    def __init__(self, center_x, center_y, radius, point): 
        self.center_x = center_x
        self.center_y = center_y
        self.radius = radius
        self.point = point
        self.x = np.zeros(point+2)
        self.y = np.zeros(point+2)
        self.x_old = np.zeros(point+2)
        self.y_old = np.zeros(point+2)        
        Bubble.total += 1
        
    def initialize_front(self):
        # determine the location of the initial spherical bubble
        for i in range(self.point+2):
            self.x[i] = self.center_x-self.radius*math.sin(2.0*math.pi*i/self.point)
            self.y[i] = self.center_y-self.radius*math.cos(2.0*math.pi*i/self.point)