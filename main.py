"""
A two-dimensional gas-liquid multiphase flows using a front-tracking type
method. A set of Navier-Stokes equation is solved on a eulerian grid 
using a second order projection method. The fluid properties are advected 
by lagrangian marker points. The time marching is second order by using 
predictor-corrector method. The code can be used to simulate a bubble 
rising in a rectangular box.
Created by: Haryo Mirsandi
"""

# import
from io_manager import IOManager
from domain import Face, Center
from fluid import Fluid
from flow_solver import FlowSolver

# clean output folder
io_man = IOManager()
io_man.create_dir('output')
io_man.clean_dir('output','bub*.png')
    
# read input file
filepath = 'input.txt'
[param, domain, fluid_prop, bubble_list] = io_man.read_input(filepath)

# initialize variables (grid, velocity, pressure, and force)
flow_solver = FlowSolver()
face = Face(domain)
center = Center(domain)
    
# initialize the physical properties inside the domain
fluid = Fluid(domain, fluid_prop)
fluid.initialize_domain(domain, center, bubble_list, fluid_prop)

# set the initial front (gas-liquid interface)
for bub in bubble_list:
    bub.initialize_front()
     
# start time-loop
# visualize the initial condition
io_man.visualize_results(face, domain, fluid, fluid_prop, bubble_list, 
                         param.time, 0)

for nstep in range(1, param.nstep+1):
    # store old variables
    face.store_old_variables()
    fluid.store_old_variables()
    for bub in bubble_list:
        bub.store_old_variables()
    
    for substep in range(2):  # second order loop
        # calculate the surface tension force at the front (lagrangian grid)
        # and distribute it to eulerian grid
        face.initialize_force(domain)
        
        for bub in bubble_list:
            bub.calculate_surface_tension(domain, fluid_prop, face)

        # update the tangential velocity at boundaries
        flow_solver.update_wall_velocity(domain, face)
        
        # calculate the (temporary) velocity
        flow_solver.calculate_temporary_velocity(param, domain, fluid_prop, 
                                                 fluid, face)
        
        # solve pressure
        flow_solver.solve_pressure()
        
        # correct the velocity by adding the pressure gradient
        flow_solver.correct_velocity()
        
        # update the front location 
        for bub in bubble_list:
            bub.update_front_location()
            
        # update physical properties
        fluid.update_density()
        fluid.update_viscosity()
        
    # substep end
    # store second order variables
    face.store_2nd_order_variables()
    fluid.store_2nd_order_variables()
    for bub in bubble_list:
        bub.store_2nd_order_variables()

    # restructure the front
    for bub in bubble_list:
        bub.restructure_front(domain)

    # visualize the results
    param.time = param.time+param.dt
    if (nstep % param.out_freq == 0):
        print('visualize results')
        io_man.visualize_results(face, domain, fluid, fluid_prop, bubble_list, 
                         param.time, nstep)        
# end time-loop
print('program finished')
