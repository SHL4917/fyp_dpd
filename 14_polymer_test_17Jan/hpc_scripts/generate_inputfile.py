#%% Imports and setup
import numpy as np
import os
from pathlib import Path
from input_func import *

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% 
conc = 0.05
bounds = 40
filebase = "pure_e/"

#%%
for i in np.arange(3, 10):
    poly_str_array = np.zeros(i) + 2
    path = filebase + "e" + str(i)
    if not os.path.exists(path):
        os.makedirs(path)
    
    poly_num = create_poly_datafile(bounds, conc, poly_str_array, 
                                    filename = path + "/data.txt")
    
    np.savetxt(path + "/e" + str(i) + ".txt", poly_num, fmt = "%.d")
    