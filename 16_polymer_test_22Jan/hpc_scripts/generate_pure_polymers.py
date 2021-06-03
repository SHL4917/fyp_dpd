import numpy as np
import sys
import os
from pathlib import Path


__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%%
poly_array = np.array([10, 15, 20, 25, 30, 35])
bounds = 40
conc = 0.05

#%%
if not os.path.exists("pure_e"):
    os.makedirs("pure_e")

#%%
for i, val in np.ndenumerate(poly_array):   
    path = f"pure_e/e{val}"
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    poly_structure = np.zeros(val, dtype = "int") + 2
    poly_num = create_poly_datafile(bounds, conc, poly_structure, 
                                    filename = path + "/data.txt")
    
    np.save(path + f"/e{val}", poly_num)