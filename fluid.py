"""
Fluid related classes
Created on Sat Jun 20 19:55:00 2020

@author: mirsandiharyo
"""

import numpy as np

class FluidProp:
    def __init__(self, cont_rho, cont_mu, disp_rho, disp_mu, sigma):
        # continuous phase
        self.cont_rho = cont_rho
        self.cont_mu = cont_mu
        # dispersed phase
        self.disp_rho = disp_rho
        self.disp_mu = disp_mu
        self.sigma = sigma
        
class Fluid:
    def __init__(self, domain, fluid_prop):
        # initialize the density and viscosity using values from the 
        # continuous phase       
        self.rho = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_rho
        self.rho_old = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_rho
        self.mu = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_mu
        self.mu_old = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_mu
        
    def initialize_domain(self, domain, center, bubble_list, fluid, fluid_prop):
        # set the fluid properties inside the discrete phase with an initial
        # spherical shape
        for i in range(2,domain.nx+1):
            for j in range(2,domain.ny+1):
                for bub in bubble_list:
                    if ((center.x[i]-bub.center_x)**2+
                        (center.y[j]-bub.center_y)**2 < bub.radius**2):
                        fluid.rho[i][j] = fluid_prop.disp_rho
                        fluid.mu[i][j]  = fluid_prop.disp_mu
