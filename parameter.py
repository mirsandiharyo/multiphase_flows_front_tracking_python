"""
Parameter class
Created on Sat Jun 20 19:48:20 2020

@author: mirsandiharyo
"""

class Parameter:
    def __init__(self, nstep, dt, max_iter, max_err, beta, out_freq):
        self.nstep = nstep
        self.dt = dt
        self.max_iter = max_iter
        self.max_err = max_err
        self.beta = beta
        self.out_freq = out_freq