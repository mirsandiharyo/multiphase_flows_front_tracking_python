"""
Flow solver
Created on Mon Jun 22 18:48:00 2020

@author: mirsandiharyo
"""

class FlowSolver:
    # def __init__(self):
    #     pass

    @staticmethod
    def update_wall_velocity(domain, face):
        # the domain is currently assumed as a box with no-slip boundary condition
        u_south = 0;
        u_north = 0;
        v_west = 0;
        v_east = 0;
        face.u[:, 0] = 2*u_south-face.u[:, 1]
        face.u[:, domain.ny+1] = 2*u_north-face.u[:, domain.ny+1];
        face.v[0, :] = 2*v_west -face.v[1, :];
        face.v[domain.nx+1, :] = 2*v_east -face.v[domain.nx+1, :];        