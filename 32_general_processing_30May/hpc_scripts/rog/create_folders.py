import numpy as np
import sys
import os
from pathlib import Path
import re
from shutil import copy2

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
lengths = [2, 3, 4, 5, 6]
folder_name = [f"rog_h{2*x:.0f}t{x:.0f}" for x in lengths]

bounds = 40
conc = 0.07

#%%
for i, folder in enumerate(folder_name):
    directory = f"{base}\{folder}"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    copy2("in.polymer_test", directory)
    copy2("polymer.pbs", directory)
    
    poly_structure = np.zeros(3 * lengths[i], dtype = "int") + 2
    poly_structure[2 * lengths[i]:] += 1
    
    poly_num = create_poly_datafile(bounds, conc, poly_structure, 
                                    filename = directory + "\data.txt",
                                    mass_array = [1, 1, 1])
    np.save(directory + f"\\rog_{lengths[i]:.0f}", poly_num)