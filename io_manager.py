"""
Input output manager
Created on Thu Jun 18 12:23:28 2020

@author: mirsandiharyo
"""

import glob, os
import numpy as np
import matplotlib.pyplot as plt
from parameter import Parameter
from domain import Domain
from fluid import FluidProp
from bubble import Bubble

class IOManager:
    @staticmethod
    def clean_dir(dir, pattern):
        """
        Clean output directory.
        """
        for file in glob.glob(dir+"/"+pattern):
            os.remove(file)
    
    @staticmethod
    def create_dir(dir):
        """
        Create output directory.
        """
        os.makedirs(dir, exist_ok=True)
    
    @staticmethod
    def read_input(filepath):
        """
        Read simulation parameters from input file.
        """
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
    
    @staticmethod
    def visualize_results(face, domain, fluid, fluid_prop, bubble_list, time, nstep):
        """
        Visualize the phase fraction field, velocity vector, and marker points.
        """
        # calculate phase fraction
        alpha = fluid.rho - fluid_prop.cont_rho
        alpha = alpha * 1/(fluid_prop.disp_rho-fluid_prop.cont_rho)     
        # plot the phase fraction
        plt.clf()
        plt.imshow(np.rot90(alpha[1:domain.nx+1,1:domain.ny+1]), cmap='jet', 
                          extent=[0,domain.lx,0,domain.ly], aspect=1)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        # set figure title
        caption = 'Time = %.3f s'% time
        plt.title(caption, fontsize=8)
        # set the colorbar
        cbar = plt.colorbar()
        cbar.ax.set_title('Phase fraction', rotation=0, size=8)
        cbar.ax.tick_params(labelsize=7)           
        # create grid and calculate velocity at the cell center
        grid_x = np.linspace(0, domain.lx, domain.nx+1)
        grid_y = np.linspace(0, domain.ly, domain.ny+1)
        u_center = np.zeros((domain.nx+1, domain.ny+1))
        v_center = np.zeros((domain.nx+1, domain.ny+1))
        u_center[0:domain.nx+1,0:domain.ny+1]=0.5*(face.u[0:domain.nx+1,1:domain.ny+2]+
                                                   face.u[0:domain.nx+1,0:domain.ny+1])
        v_center[0:domain.nx+1,0:domain.ny+1]=0.5*(face.v[1:domain.nx+2,0:domain.ny+1]+
                                                   face.v[0:domain.nx+1,0:domain.ny+1])
        # plot the velocity vector
        plt.quiver(grid_x, grid_y, np.rot90(u_center), np.rot90(v_center))
        # plot the marker points
        for bub in bubble_list:
            plt.plot(bub.x[0:bub.point],bub.y[0:bub.point],'k',linewidth=1)  
        # save the plot
        caption = 'output/bub_%03d.png' % nstep
        plt.savefig(caption,dpi=150)
        plt.pause(0.001)