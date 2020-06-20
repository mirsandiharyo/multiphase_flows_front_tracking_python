"""
Bubble class
Created on Sat Jun 20 19:06:17 2020

@author: mirsandiharyo
"""

import numpy as np

class Bubble:
    total = 0
    def __init__(self, center_x, center_y, radius, point): 
        self.center_x = center_x
        self.center_y = center_y
        self.radius = radius
        self.point = point
        self.x = np.zeros(point)
        self.y = np.zeros(point)
        self.x_old = np.zeros(point)
        self.y_old = np.zeros(point)        
        Bubble.total += 1
