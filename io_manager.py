"""
Input output manager
Created on Thu Jun 18 12:23:28 2020

@author: mirsandiharyo
"""

import glob, os

def clean_dir(dir, pattern):
    for file in glob.glob(dir+"/"+pattern):
        os.remove(file)
        
def create_dir(dir):
    os.makedirs(dir, exist_ok=True)
