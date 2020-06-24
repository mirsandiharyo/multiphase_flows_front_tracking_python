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
        self.rho_old = self.rho.copy()
        self.mu_old = self.mu.copy()
 
    def store_2nd_order_variables(self):
        """ 
        Store second order variables.
        """
        self.rho = 0.5*(self.rho+self.rho_old)
        self.mu = 0.5*(self.mu+self.mu_old)
        
    def update_density(self, param, domain, bubble_list, fluid_prop):
        """ 
        Update the density field using the density jump at the lagrangian 
        interface.
        Linear averaging is used to get the value at each cell.
        """
        # initialize the variables to store the density jump
        face_x = np.zeros((domain.nx+2, domain.ny+2))
        face_y = np.zeros((domain.nx+2, domain.ny+2))
        # distribute the density jump to the eulerian grid
        for bub in bubble_list:
            for i in range(1, bub.point+1):
                # density jump in x-direction
                force_x = -0.5*(bub.y[i+1]-bub.y[i-1])* \
                    (fluid_prop.disp_rho-fluid_prop.cont_rho)
                bub.distribute_lagrangian_to_eulerian(domain, face_x, bub.x[i],
                                                      bub.y[i], force_x, 1)
                # density jump in y-direction
                force_y = 0.5*(bub.x[i+1]-bub.x[i-1])* \
                    (fluid_prop.disp_rho-fluid_prop.cont_rho); 
                bub.distribute_lagrangian_to_eulerian(domain, face_y, bub.x[i], 
                                                      bub.y[i], force_y, 2)
            
        # construct the density field using SOR
        # TODO: create SOR function 
        for it in range(param.max_iter):
            old_rho = self.rho.copy()
            for i in range(1, domain.nx+1):
                for j in range(1, domain.ny+1):
                    self.rho[i,j] = (1.0-param.beta)* \
                         self.rho[i  ,j  ]+ param.beta*0.25* \
                        (self.rho[i+1,j  ]+self.rho[i-1,j  ]+
                         self.rho[i  ,j+1]+self.rho[i  ,j-1]+
                         domain.dx*face_x[i-1,j  ]-
                         domain.dx*face_x[i  ,j  ]+
                         domain.dy*face_y[i  ,j-1]-
                         domain.dy*face_y[i  ,j  ])
            if (np.abs(old_rho-self.rho).max() < param.max_err):
                break
            
    def update_viscosity(self, fluid_prop):
        """
        Update the viscosity field using harmonic averaging.
        """
        self.mu = self.rho-fluid_prop.cont_rho
        self.mu = self.mu*(fluid_prop.disp_mu -fluid_prop.cont_mu )/ \
                          (fluid_prop.disp_rho-fluid_prop.cont_rho)
        self.mu = self.mu+fluid_prop.cont_mu