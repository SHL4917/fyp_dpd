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
lengths = [5, 10, 15, 20, 25, 30, 35]
folder_name = [f"e2e_{x:.0f}" for x in lengths]

bounds = 40
conc = 0.05

#%%
for i, folder in enumerate(folder_name):
    directory = f"{base}\{folder}"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    copy2("in.polymer_test", directory)
    copy2("polymer.pbs", directory)
    
    poly_structure = np.zeros(lengths[i], dtype = "int") + 2
    poly_num = create_poly_datafile(bounds, conc, poly_structure, 
                                    filename = directory + "\data.txt",
                                    mass_array = [1, 1, 1])
    np.save(directory + f"\length_{lengths[i]:.0f}", poly_num)