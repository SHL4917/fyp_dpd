import numpy as np
import sys
import os
from pathlib import Path


__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%% Get polymer coordinates for each folder and save in each folder
for i in range(10, 36, 5):
    polyfile = np.load(f"e{i}/e{i}.npy")
    
    np.savetxt(f"e{i}/e{i}.txt", polyfile, fmt = "%.d")
    polymer_coord = get_poly_from_dump(f"e{i}/dump.atom", f"e{i}/e{i}.txt")
    np.save("e" + str(i) + "/polymer_coords.npy", polymer_coord)
