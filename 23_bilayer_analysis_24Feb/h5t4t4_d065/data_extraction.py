import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%%
polyfile = "polymer_coords.npy"
dump_poly = "dump.polymeronly"
dump_all = "dump.allparticles"

head = np.array([2, 2, 2, 2, 2])
tail1 = np.array([3, 3, 3, 3])
tail2 = np.array([3, 3, 3, 3])
polymer_struct = np.concatenate((head, tail1, tail2))

polymer_masses = np.array([1, 1, 1])

resolution = 400
bounds = 30

#%%
poly_coords = get_poly_from_dump(dump_poly, polyfile)
poly_com = get_com_all(poly_coords, polymer_struct, polymer_masses)
np.save("poly_com.npy", poly_com)

#%%

poly_com_last = poly_com[:, :, -1]
planes = get_planes(poly_com_last, resolution, bounds, 1)

np.save("plane.npy", planes)

#%%
coord1, coord2, coord3 = get_allparticles_from_dump2(dump_all)

np.save("coord1.npy", coord1)
np.save("coord2.npy", coord2)
np.save("coord3.npy", coord3)
