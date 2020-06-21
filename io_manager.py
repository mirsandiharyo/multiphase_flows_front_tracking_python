"""
Input output manager
Created on Thu Jun 18 12:23:28 2020

@author: mirsandiharyo
"""

import glob, os
from parameter import Parameter
from domain import Domain
from fluid import FluidProp
from bubble import Bubble

def clean_dir(dir, pattern):
    for file in glob.glob(dir+"/"+pattern):
        os.remove(file)
        
def create_dir(dir):
    os.makedirs(dir, exist_ok=True)

def read_input(filepath):
    with open(filepath) as file:
        # solver parameters
        file.readline()
        nstep = int(file.readline().split("=")[1])    
        dt = float(file.readline().split("=")[1])  
        max_iter = int(file.readline().split("=")[1]) 
        max_err = float(file.readline().split("=")[1])   
        beta = float(file.readline().split("=")[1])
        out_freq = int(file.readline().split("=")[1])
        param = Parameter(nstep, dt, max_iter, max_err, beta, out_freq)
        file.readline()    
        # numerical parameters
        file.readline()
        lx = float(file.readline().split("=")[1])    
        ly = float(file.readline().split("=")[1])  
        nx = int(file.readline().split("=")[1]) 
        ny = int(file.readline().split("=")[1])   
        gravx = float(file.readline().split("=")[1])
        gravy = float(file.readline().split("=")[1])
        domain = Domain(lx, ly, nx, ny, gravx, gravy)
        file.readline()
        # physical properties
        # dispersed phase
        file.readline()
        file.readline()
        disp_rho = float(file.readline().split("=")[1])
        disp_mu = float(file.readline().split("=")[1])
        sigma = float(file.readline().split("=")[1])
        # continuous phase
        file.readline()    
        cont_rho = float(file.readline().split("=")[1])
        cont_mu = float(file.readline().split("=")[1])
        fluid_prop = FluidProp(cont_rho, cont_mu, disp_rho, disp_mu, sigma)
        file.readline()
        # bubble size and location
        bubble_list = []
        file.readline()
        nbub = int(file.readline().split("=")[1])
        for num in range(nbub):
            radius = float(file.readline().split("=")[1])
            center_x = float(file.readline().split("=")[1])
            center_y = float(file.readline().split("=")[1])
            point = int(file.readline().split("=")[1])
            bubble_list.append(Bubble(center_x, center_y, radius, point))
    return param, domain, fluid_prop, bubble_list