"""
Fluid related classes
Created on Sat Jun 20 19:55:00 2020

@author: mirsandiharyo
"""

import numpy as np

class FluidProp:
    def __init__(self, cont_rho, cont_mu, disp_rho, disp_mu, sigma):
        """ 
        Initialize the fluid properties of the continuous and dispersed phases.
        """
        # continuous phase
        self.cont_rho = cont_rho
        self.cont_mu = cont_mu
        # dispersed phase
        self.disp_rho = disp_rho
        self.disp_mu = disp_mu
        self.sigma = sigma
        
class Fluid:
    def __init__(self, domain, fluid_prop):
        """ 
        Initialize the density and viscosity fields using the properties from
        continuous phase.          
        """
        self.rho = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_rho
        self.rho_old = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_rho
        self.mu = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_mu
        self.mu_old = np.zeros((domain.nx+2, domain.ny+2))+fluid_prop.cont_mu
        
    def initialize_domain(self, domain, center, bubble_list, fluid_prop):
        """ 
        Set the fluid properties inside the discrete phase with an initial
        spherical shape. 
        """
        for i in range(1,domain.nx+1):
            for j in range(1,domain.ny+1):
                for bub in bubble_list:
                    if ((center.x[i]-bub.center_x)**2+
                        (center.y[j]-bub.center_y)**2 < bub.radius**2):
                        self.rho[i,j] = fluid_prop.disp_rho
                        self.mu[i,j]  = fluid_prop.disp_mu
                           
    def store_old_variables(self):
        """ 
        Store old variables for second order scheme.
        """
        self.rho_old = self.rho
        self.mu_old = self.mu
 
    def store_2nd_order_variables(self):
        """ 
        Store second order variables.
        """
        self.rho = 0.5*(self.rho+self.rho_old)
        self.mu = 0.5*(self.mu+self.mu_old)
        
    def update_density(self):
        """ 
        Update the density field using the density jump at the lagrangian interface.
        Linear averaging is used to get the value at each cell.
        """
        pass
    
    def update_viscosity(self, fluid_prop):
        """
        Update the viscosity field using harmonic averaging.
        """
        self.mu = self.rho-fluid_prop.cont_rho
        self.mu = self.mu*(fluid_prop.disp_mu -fluid_prop.cont_mu )/ \
                          (fluid_prop.disp_rho-fluid_prop.cont_rho)
        self.mu = self.mu+fluid_prop.cont_mu