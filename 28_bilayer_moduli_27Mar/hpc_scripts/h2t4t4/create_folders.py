import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import shutil

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
bounds = 36
head = np.array([2, 2])
tail1 = np.array([3, 3, 3, 3])
tail2 = np.array([3, 3, 3, 3])
dense = np.linspace(65, 85, 20)

#%%
for i in dense:
    directory = f"{base}\conc_{i * 100:.0f}"
    filename = directory + "\data.txt"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    shutil.copy2("in.bilayer_test", directory)
    shutil.copy2("polymer.pbs", directory)
    shutil.copy2("extract_data.py", directory)
    shutil.copy2("extract_data.pbs", directory)
    shutil.copy2("fyp_functions_short.py", directory)
    
    polymer_coords, conc, poly_count = create_bilayer_datafile_v2(bounds,
                                                                  head, 
                                                                  tail1, 
                                                                  tail2,
                                                                  i/100,
                                                                  filename)
    
    np.save(directory + "\polymer_coords", polymer_coords)
    np.savetxt(directory + "\concentration.txt", np.array([conc]))
    np.savetxt(directory + "\polymer_count.txt", np.array([poly_count]))