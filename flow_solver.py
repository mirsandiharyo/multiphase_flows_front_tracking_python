"""
Flow solver
Created on Mon Jun 22 18:48:00 2020

@author: mirsandiharyo
"""

import numpy as np

class FlowSolver:
    @staticmethod
    def update_wall_velocity(domain, face):
        """
        Update the wall velocity (the domain is currently assumed as a box 
        with no-slip boundary condition).
        """ 
        u_south = 0;
        u_north = 0;
        v_west = 0;
        v_east = 0;
        face.u[:,0] = 2*u_south-face.u[:,1]
        face.v[0,:] = 2*v_west -face.v[1,:]
        face.u[:,domain.ny+1] = 2*u_north-face.u[:,domain.ny]
        face.v[domain.nx+1,:] = 2*v_east -face.v[domain.nx,:]
      
    @staticmethod
    def calculate_temporary_velocity(param, domain, fluid_prop, fluid, face):
        """
        Calculate the temporary velocities without accounting for the pressure
        (first step of the second order projection method).
        """ 
        # temporary u velocity (advection term)
        face.u_temp[1:domain.nx,1:domain.ny+1] = \
        	   face.u[1:domain.nx  ,1:domain.ny+1]+param.dt*(-0.25*
            (((face.u[2:domain.nx+1,1:domain.ny+1]+
               face.u[1:domain.nx  ,1:domain.ny+1])**2-
              (face.u[1:domain.nx  ,1:domain.ny+1]+
               face.u[0:domain.nx-1,1:domain.ny+1])**2)/domain.dx+
             ((face.u[1:domain.nx  ,2:domain.ny+2]+
               face.u[1:domain.nx  ,1:domain.ny+1])*
              (face.v[2:domain.nx+1,1:domain.ny+1]+
               face.v[1:domain.nx  ,1:domain.ny+1])-
              (face.u[1:domain.nx  ,1:domain.ny+1]+
               face.u[1:domain.nx  ,0:domain.ny  ])*
              (face.v[2:domain.nx+1,0:domain.ny  ]+
               face.v[1:domain.nx  ,0:domain.ny  ]))/domain.dy)+
               face.force_x[1:domain.nx,1:domain.ny+1]/
              (0.5*(fluid.rho[2:domain.nx+1,1:domain.ny+1]+
                    fluid.rho[1:domain.nx  ,1:domain.ny+1]))-
              (1.0 -fluid_prop.cont_rho/
              (0.5*(fluid.rho[2:domain.nx+1,1:domain.ny+1]+
                    fluid.rho[1:domain.nx  ,1:domain.ny+1])))*domain.gravx)
                       
        # temporary v velocity (advection term)
        face.v_temp[1:domain.nx+1,1:domain.ny] = \
               face.v[1:domain.nx+1,1:domain.ny  ]+param.dt*(-0.25*
            (((face.u[1:domain.nx+1,2:domain.ny+1]+
               face.u[1:domain.nx+1,1:domain.ny  ])*
              (face.v[2:domain.nx+2,1:domain.ny  ]+
               face.v[1:domain.nx+1,1:domain.ny  ])-
              (face.u[0:domain.nx  ,2:domain.ny+1]+
               face.u[0:domain.nx  ,1:domain.ny  ])*
              (face.v[1:domain.nx+1,1:domain.ny  ]+
               face.v[0:domain.nx  ,1:domain.ny  ]))/domain.dx+
             ((face.v[1:domain.nx+1,2:domain.ny+1]+
               face.v[1:domain.nx+1,1:domain.ny  ])**2-
              (face.v[1:domain.nx+1,1:domain.ny]+
               face.v[1:domain.nx+1,0:domain.ny-1])**2)/domain.dy)+
               face.force_y[1:domain.nx+1,1:domain.ny]/
              (0.5*(fluid.rho[1:domain.nx+1,2:domain.ny+1]+
                    fluid.rho[1:domain.nx+1,1:domain.ny  ]))-
              (1.0 -fluid_prop.cont_rho/
              (0.5*(fluid.rho[1:domain.nx+1,2:domain.ny+1]+
                    fluid.rho[1:domain.nx+1,1:domain.ny  ])))*domain.gravy)
                       
        # temporary u velocity (diffusion term)
        face.u_temp[1:domain.nx  ,1:domain.ny+1] = \
        face.u_temp[1:domain.nx  ,1:domain.ny+1]+param.dt*((1./domain.dx)*2.*
          (fluid.mu[2:domain.nx+1,1:domain.ny+1]*(1./domain.dx)*
            (face.u[2:domain.nx+1,1:domain.ny+1]-
             face.u[1:domain.nx  ,1:domain.ny+1])-
           fluid.mu[1:domain.nx  ,1:domain.ny+1]*(1./domain.dx)*
            (face.u[1:domain.nx  ,1:domain.ny+1]-
             face.u[0:domain.nx-1,1:domain.ny+1]))+(1./domain.dy)*(0.25*
          (fluid.mu[1:domain.nx  ,1:domain.ny+1]+
           fluid.mu[2:domain.nx+1,1:domain.ny+1]+         
           fluid.mu[2:domain.nx+1,2:domain.ny+2]+
           fluid.mu[1:domain.nx  ,2:domain.ny+2])*((1./domain.dy)*
            (face.u[1:domain.nx  ,2:domain.ny+2]-
             face.u[1:domain.nx  ,1:domain.ny+1])+(1./domain.dx)*
            (face.v[2:domain.nx+1,1:domain.ny+1]-
             face.v[1:domain.nx  ,1:domain.ny+1]))-0.25*
          (fluid.mu[1:domain.nx  ,1:domain.ny+1]+
           fluid.mu[2:domain.nx+1,1:domain.ny+1]+
           fluid.mu[2:domain.nx+1,0:domain.ny  ]+
           fluid.mu[1:domain.nx  ,0:domain.ny  ])*((1./domain.dy)*
            (face.u[1:domain.nx  ,1:domain.ny+1]-
             face.u[1:domain.nx  ,0:domain.ny  ])+(1./domain.dx)*
            (face.v[2:domain.nx+1,0:domain.ny  ]-
             face.v[1:domain.nx  ,0:domain.ny  ]))))/(0.5*
         (fluid.rho[2:domain.nx+1,1:domain.ny+1]+
          fluid.rho[1:domain.nx  ,1:domain.ny+1]))     
                                                       
        # temporary v velocity (diffusion term)
        face.v_temp[1:domain.nx+1,1:domain.ny  ] = \
        face.v_temp[1:domain.nx+1,1:domain.ny  ]+param.dt*((1./domain.dx)*(0.25*
          (fluid.mu[1:domain.nx+1,1:domain.ny  ]+
           fluid.mu[2:domain.nx+2,1:domain.ny  ]+
           fluid.mu[2:domain.nx+2,2:domain.ny+1]+
           fluid.mu[1:domain.nx+1,2:domain.ny+1])*((1./domain.dy)*
            (face.u[1:domain.nx+1,2:domain.ny+1]-
             face.u[1:domain.nx+1,1:domain.ny  ])+(1./domain.dx)*
            (face.v[2:domain.nx+2,1:domain.ny  ]-
             face.v[1:domain.nx+1,1:domain.ny  ]))-0.25*      
          (fluid.mu[1:domain.nx+1,1:domain.ny  ]+
           fluid.mu[1:domain.nx+1,2:domain.ny+1]+                  
           fluid.mu[0:domain.nx  ,2:domain.ny+1]+
           fluid.mu[0:domain.nx  ,1:domain.ny  ])*((1./domain.dy)*
            (face.u[0:domain.nx  ,2:domain.ny+1]-
             face.u[0:domain.nx  ,1:domain.ny  ])+(1./domain.dx)*
            (face.v[1:domain.nx+1,1:domain.ny  ]-
             face.v[0:domain.nx  ,1:domain.ny  ])))+(1./domain.dy)*2.*
          (fluid.mu[1:domain.nx+1,2:domain.ny+1]*(1./domain.dy)*
            (face.v[1:domain.nx+1,2:domain.ny+1]-
             face.v[1:domain.nx+1,1:domain.ny  ])-
           fluid.mu[1:domain.nx+1,1:domain.ny  ]*(1./domain.dy)*
            (face.v[1:domain.nx+1,1:domain.ny  ]-
             face.v[1:domain.nx+1,0:domain.ny-1])))/(0.5*
         (fluid.rho[1:domain.nx+1,2:domain.ny+1]+
          fluid.rho[1:domain.nx+1,1:domain.ny  ])) 
                                                              
    @staticmethod
    def solve_pressure(param, domain, fluid, face, center):
        """
        Calculate the pressure field.
        """
        # initialize variables
        temp1 = np.zeros((domain.nx+2, domain.ny+2))
        temp2 = np.zeros((domain.nx+2, domain.ny+2))
        
        # calculate source term and the coefficient for pressure
        rho_temp = fluid.rho.copy()
        large_num = 1000;
        rho_temp[:,0] = large_num;
        rho_temp[:,domain.ny+1] = large_num;
        rho_temp[0,:] = large_num;
        rho_temp[domain.nx+1,:] = large_num;
        
        temp1[1:domain.nx+1,1:domain.ny+1] = (0.5/param.dt)* \
            ((face.u_temp[1:domain.nx+1,1:domain.ny+1]
             -face.u_temp[0:domain.nx  ,1:domain.ny+1])/domain.dx+
             (face.v_temp[1:domain.nx+1,1:domain.ny+1]
             -face.v_temp[1:domain.nx+1,0:domain.ny  ])/domain.dy)
        
        temp2[1:domain.nx+1,1:domain.ny+1] = 1.0/((1./domain.dx)*
             (1./(domain.dx*
             (rho_temp[2:domain.nx+2,1:domain.ny+1]+
              rho_temp[1:domain.nx+1,1:domain.ny+1]))+
              1./(domain.dx*
             (rho_temp[0:domain.nx  ,1:domain.ny+1]+
              rho_temp[1:domain.nx+1,1:domain.ny+1])))+(1./domain.dy)*
             (1./(domain.dy*
             (rho_temp[1:domain.nx+1,2:domain.ny+2]+
              rho_temp[1:domain.nx+1,1:domain.ny+1]))+
              1./(domain.dy*
             (rho_temp[1:domain.nx+1,0:domain.ny  ]+
              rho_temp[1:domain.nx+1,1:domain.ny+1]))))  
        
        # construct the pressure field using SOR
        # TODO: create SOR function
        for it in range(param.max_iter):
            old_pres = center.pres.copy()
            for iskip in range(2):
                rb = iskip
                center.pres[1+rb:domain.nx+1:2,1:domain.ny+1:2] = \
                ((1.0-param.beta)*
                center.pres[1+rb:domain.nx+1:2,1:domain.ny+1:2]+param.beta*
                      temp2[1+rb:domain.nx+1:2,1:domain.ny+1:2]*
               ((1.0/domain.dx/domain.dx)* 
               (center.pres[2+rb:domain.nx+2:2,1:domain.ny+1:2]/
                  (rho_temp[2+rb:domain.nx+2:2,1:domain.ny+1:2]+
                   rho_temp[1+rb:domain.nx+1:2,1:domain.ny+1:2])+
                center.pres[  rb:domain.nx  :2,1:domain.ny+1:2]/
                  (rho_temp[  rb:domain.nx  :2,1:domain.ny+1:2]+
                   rho_temp[1+rb:domain.nx+1:2,1:domain.ny+1:2]))+
               (1.0/domain.dy/domain.dy)*
               (center.pres[1+rb:domain.nx+1:2,2:domain.ny+2:2]/
                  (rho_temp[1+rb:domain.nx+1:2,2:domain.ny+2:2]+
                   rho_temp[1+rb:domain.nx+1:2,1:domain.ny+1:2])+
                center.pres[1+rb:domain.nx+1:2,0:domain.ny  :2]/
                  (rho_temp[1+rb:domain.nx+1:2,0:domain.ny  :2]+
                   rho_temp[1+rb:domain.nx+1:2,1:domain.ny+1:2]))-
                      temp1[1+rb:domain.nx+1:2,1:domain.ny+1:2]))
    
                rb=1-iskip
                center.pres[1+rb:domain.nx+1:2,2:domain.ny+1:2] = \
                ((1.0-param.beta)*
                center.pres[1+rb:domain.nx+1:2,2:domain.ny+1:2]+param.beta*
                      temp2[1+rb:domain.nx+1:2,2:domain.ny+1:2]*
               ((1.0/domain.dx/domain.dx)*
               (center.pres[2+rb:domain.nx+2:2,2:domain.ny+1:2]/
                  (rho_temp[2+rb:domain.nx+2:2,2:domain.ny+1:2]+
                   rho_temp[1+rb:domain.nx+1:2,2:domain.ny+1:2])+
                center.pres[  rb:domain.nx:2  ,2:domain.ny+1:2]/
                  (rho_temp[  rb:domain.nx:2  ,2:domain.ny+1:2]+
                   rho_temp[1+rb:domain.nx+1:2,2:domain.ny+1:2]))+
               (1.0/domain.dy/domain.dy)*
               (center.pres[1+rb:domain.nx+1:2,3:domain.ny+2:2]/
                  (rho_temp[1+rb:domain.nx+1:2,3:domain.ny+2:2]+
                   rho_temp[1+rb:domain.nx+1:2,2:domain.ny+1:2])+
                center.pres[1+rb:domain.nx+1:2,1:domain.ny  :2]/
                  (rho_temp[1+rb:domain.nx+1:2,1:domain.ny  :2]+
                   rho_temp[1+rb:domain.nx+1:2,2:domain.ny+1:2]))-
                      temp1[1+rb:domain.nx+1:2,2:domain.ny+1:2]))
    
            center.pres[0,:] = center.pres[1,:]; 
            center.pres[domain.nx+1,:] = center.pres[domain.nx,:]
            center.pres[:,0] = center.pres[:,1]; 
            center.pres[:,domain.ny+1] = center.pres[:,domain.ny]
                                  
            if (np.abs(old_pres-center.pres).max() < param.max_err):
                break

    @staticmethod
    def correct_velocity(param, domain, fluid, face, center):
        """
        Correct the velocity by adding the pressure gradient.
        """
        # correct velocity in x-direction
        face.u[1:domain.nx,1:domain.ny+1] = \
            face.u_temp[1:domain.nx  ,1:domain.ny+1]-param.dt*(2.0/domain.dx)* \
           (center.pres[2:domain.nx+1,1:domain.ny+1]-
            center.pres[1:domain.nx  ,1:domain.ny+1])/ \
             (fluid.rho[2:domain.nx+1,1:domain.ny+1]+
              fluid.rho[1:domain.nx  ,1:domain.ny+1])
                
        # correct velocity in y-direction
        face.v[1:domain.nx+1,1:domain.ny] = \
            face.v_temp[1:domain.nx+1,1:domain.ny  ]-param.dt*(2.0/domain.dy)* \
           (center.pres[1:domain.nx+1,2:domain.ny+1]-
            center.pres[1:domain.nx+1,1:domain.ny  ])/ \
             (fluid.rho[1:domain.nx+1,2:domain.ny+1]+
              fluid.rho[1:domain.nx+1,1:domain.ny  ])