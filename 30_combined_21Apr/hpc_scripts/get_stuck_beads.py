import numpy as np
import sys
import os
from pathlib import Path
# import matplotlib.pyplot as plt


__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append(base)
from fyp_functions_short import *

#%%
def get_hist(coords, num_bins, bounds, index = 2, density = True):

    axis = np.linspace(0, bounds[index], num_bins)
    count = np.zeros(num_bins)
    height_range = axis.max() - axis.min()    
    
    for i in range(len(coords)):
        height = coords[i, index]
        height_index = (height/height_range) * num_bins
        height_index = np.trunc(height_index).astype(int) - 1
        count[height_index] += 1
    
    vol = bounds[0] * bounds[1] * bounds[2] * (1/num_bins)
    
    if density:
        count /= vol
    return count, axis
    
#%%
dump = "dump.allparticles"
timestep, data = get_data_from_dump(dump, 5, 2)

np.save("coord1.npy", data[1])
np.save("coord2.npy", data[2])
np.save("coord3.npy", data[3])
np.save("coord4.npy", data[4])
np.save("coord5.npy", data[5])

#%%
coord1 = np.load("coord1.npy")
coord2 = np.load("coord2.npy")
coord3 = np.load("coord3.npy")
coord4 = np.load("coord4.npy")
coord5 = np.load("coord5.npy")

bounds = np.array([36, 36, 36 * 3])

file = os.path.basename(os.getcwd())

bilayer_struct = "$H_5 T_4 T_4$"
micelle_struct = "$H_6 T_4$"

val = file.split("_")[-1]
varied = "$x_1=" + val + "$"
filename = "x1_" + val

#%%
num_slices = np.shape(coord1)[2]
bilayer_inner_bounds = np.array([19, 80])

lower_ind = axis[axis < bilayer_inner_bounds[0]].size - 1
upper_ind = axis.size - axis[axis > bilayer_inner_bounds[1]].size

stuck4 = np.zeros(num_slices)
stuck5 = np.zeros(num_slices)

for i in range(num_slices):
    snap4 = coord4[:, :, i]
    snap5 = coord5[:, :, i]
    count4, axis = get_hist(snap4, 70, bounds, density = False)
    count5, axis = get_hist(snap5, 70, bounds, density = False)
    
    stuck4[i] = count4[0:lower_ind + 1].sum() + count4[upper_ind:].sum()
    stuck5[i] = count5[0:lower_ind + 1].sum() + count5[upper_ind:].sum()

np.save(f"stuck4_{filename}.npy", stuck4)    
np.save(f"stuck5_{filename}.npy", stuck5)
np.save(f"timestep_{filename}.npy", timestep)

#%%
num4 = len(coord4)
num5 = len(coord5)

np.savetxt(f"particle_number_4_5_{filename}.npy", np.array([num4, num5]))