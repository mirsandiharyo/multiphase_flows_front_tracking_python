"""
Flow solver
Created on Mon Jun 22 18:48:00 2020

@author: mirsandiharyo
"""

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
        face.u[:, 0] = 2*u_south-face.u[:, 1]
        face.u[:, domain.ny+1] = 2*u_north-face.u[:, domain.ny+1];
        face.v[0, :] = 2*v_west -face.v[1, :];
        face.v[domain.nx+1, :] = 2*v_east -face.v[domain.nx+1, :];
        
    @staticmethod
    def calculate_temporary_velocity(param, domain, fluid_prop, fluid, face):
        """
        Calculate the temporary velocities without accounting for the pressure
        (first step of the second order projection method).
        """ 
        # temporary u velocity (advection term)
        for i in range(1, domain.nx):
            for j in range(1, domain.ny+1):
                face.u_temp[i,j] = face.u[i,j]+param.dt*(-0.25*
                    (((face.u[i+1,j  ]+face.u[i  ,j  ])**2-
                      (face.u[i  ,j  ]+face.u[i-1,j  ])**2)/domain.dx+
                     ((face.u[i  ,j+1]+face.u[i  ,j  ])*
                      (face.v[i+1,j  ]+face.v[i  ,j  ])-
                      (face.u[i  ,j  ]+face.u[i  ,j-1])*
                      (face.v[i+1,j-1]+face.v[i  ,j-1]))/domain.dy)+
                       face.force_x[i,j]/
                      (0.5*(fluid.rho[i+1,j]+fluid.rho[i,j]))-
                      (1.0 -fluid_prop.cont_rho/
                      (0.5*(fluid.rho[i+1,j]+fluid.rho[i,j])))*domain.gravx)
        # temporary v velocity (advection term)
        for i in range(1, domain.nx+1):
            for j in range(1, domain.ny):
                face.v_temp[i,j] = face.v[i,j]+param.dt*(-0.25*
                    (((face.u[i  ,j+1]+face.u[i  ,j  ])*
                      (face.v[i+1,j  ]+face.v[i  ,j  ])-
                      (face.u[i-1,j+1]+face.u[i-1,j  ])*
                      (face.v[i  ,j  ]+face.v[i-1,j  ]))/domain.dx+
                     ((face.v[i  ,j+1]+face.v[i  ,j  ])**2-
                      (face.v[i  ,j  ]+face.v[i  ,j-1])**2)/domain.dy)+
                       face.force_y[i,j]/
                      (0.5*(fluid.rho[i,j+1]+fluid.rho[i,j]))-
                      (1.0 -fluid_prop.cont_rho/
                      (0.5*(fluid.rho[i,j+1]+fluid.rho[i,j])))*domain.gravy)
        # temporary u velocity (diffusion term)
        for i in range(1, domain.nx):
            for j in range (1, domain.ny+1):
                face.u_temp[i,j] = face.u_temp[i,j]+param.dt*((1./domain.dx)*2.*
                (fluid.mu[i+1,j  ]*(1./domain.dx)*(face.u[i+1,j  ]-face.u[i  ,j  ]) -
                 fluid.mu[i  ,j  ]*(1./domain.dx)*(face.u[i  ,j  ]-face.u[i-1,j  ]))+
                 (1./domain.dy)*(0.25*
                (fluid.mu[i  ,j  ]+fluid.mu[i+1,j  ]+         
                 fluid.mu[i+1,j+1]+fluid.mu[i  ,j+1])*                            
                ((1./domain.dy)*(face.u[i  ,j+1]-face.u[i  ,j  ])+
                 (1./domain.dx)*(face.v[i+1,j  ]-face.v[i  ,j  ]))-0.25*
                (fluid.mu[i  ,j  ]+fluid.mu[i+1,j  ]+
                 fluid.mu[i+1,j-1]+fluid.mu[i  ,j-1])*
                ((1./domain.dy)*(face.u[i  ,j  ]-face.u[i  ,j-1])+
                 (1./domain.dx)*(face.v[i+1,j-1]-face.v[i  ,j-1]))))/ \
                 (0.5*(fluid.rho[i+1,j]+fluid.rho[i,j]))
        # temporary v velocity (diffusion term)
        for i in range(1, domain.nx+1):
            for j in range(1, domain.ny):
                face.v_temp[i,j] = face.v_temp[i,j]+param.dt*((1./domain.dx)*(0.25*
                (fluid.mu[i  ,j  ]+fluid.mu[i+1,j  ]+
                 fluid.mu[i+1,j+1]+fluid.mu[i,j+1  ])*
                ((1./domain.dy)*(face.u[i  ,j+1]-face.u[i  ,j  ])+            
                 (1./domain.dx)*(face.v[i+1,j  ]-face.v[i  ,j  ]))-0.25*      
                (fluid.mu[i  ,j  ]+fluid.mu[i  ,j+1]+                  
                 fluid.mu[i-1,j+1]+fluid.mu[i-1,j  ])*          
                ((1./domain.dy)*(face.u[i-1,j+1]-face.u[i-1,j  ])+
                 (1./domain.dx)*(face.v[i  ,j  ]-face.v[i-1,j  ])))+ \
                 (1./domain.dy)*2.* \
                (fluid.mu[i  ,j+1]*(1./domain.dy)*(face.v[i  ,j+1]-face.v[i  ,j  ])-
                 fluid.mu[i  ,j  ]*(1./domain.dy)*(face.v[i  ,j  ]-face.v[i  ,j-1])))/ \
                 (0.5*(fluid.rho[i,j+1]+fluid.rho[i,j]))
                 
    @staticmethod
    def solve_pressure():
        """
        Calculate the pressure field.
        """
        pass
    
    @staticmethod
    def correct_velocity(param, domain, fluid, face, center):
        """
        Correct the velocity by adding the pressure gradient.
        """
        # correct velocity in x-direction
        face.u[1:domain.nx, 1:domain.ny+1] = \
            face.u_temp[1:domain.nx, 1:domain.ny+1]-param.dt*(2.0/domain.dy)* \
            (center.pres[2:domain.nx+1, 1:domain.ny+1]-center.pres[1:domain.nx, 1:domain.ny+1])/ \
            (fluid.rho[2:domain.nx+1, 1:domain.ny+1]  +fluid.rho[1:domain.nx, 1:domain.ny+1])
                    
        # correct velocity in y-direction
        face.v[1:domain.nx+1, 1:domain.ny] = \
            face.v_temp[1:domain.nx+1, 1:domain.ny]-param.dt*(2.0/domain.dy)* \
            (center.pres[1:domain.nx+1, 2:domain.ny+1]-center.pres[1:domain.nx+1, 1:domain.ny])/ \
            (fluid.rho[1:domain.nx+1, 2:domain.ny+1]  +fluid.rho[1:domain.nx+1, 1:domain.ny])           



