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

# clean output folder
io_man = IOManager()
io_man.create_dir('output')
io_man.clean_dir('output','bub*.png')
    
# read input file
filepath = 'input.txt'
[param, domain, fluid_prop, bubble_list] = io_man.read_input(filepath)

# initialize variables (grid, velocity, pressure, and force)
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

for nstep in range(param.nstep):
    # store second order variables
    face.store_old_variables()
    fluid.store_old_variables()
    for bub in bubble_list:
        bub.store_old_variables()
    
    for substep in range(2):  # second order loop
        pass
        # calculate the surface tension force at the front (lagrangian grid)
        # and distribute it to eulerian grid

        # update the tangential velocity at boundaries

        # calculate the (temporary) velocity

        # solve pressure
        
        # correct the velocity by adding the pressure gradient
        
        # update the front location 

        # update physical properties
  
    # end
    # store second order variables
    
    # restructure the front
    
    # visualize the results

 # end time-loop

