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
coord1, coord2, coord3, coord4, coord5 = get_allparticles_from_dump3(dump)

np.save("coord1.npy", coord1)
np.save("coord2.npy", coord2)
np.save("coord3.npy", coord3)
np.save("coord4.npy", coord4)
np.save("coord5.npy", coord5)

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
index = -1
snap1 = coord1[:, :, index]
snap2 = coord2[:, :, index]
snap3 = coord3[:, :, index]
snap4 = coord4[:, :, index]
snap5 = coord5[:, :, index]

#%%
count1, axis = get_hist(snap1, 70, bounds)
count2, axis = get_hist(snap2, 70, bounds)
count3, axis = get_hist(snap3, 70, bounds)
count4, axis = get_hist(snap4, 70, bounds)
count5, axis = get_hist(snap5, 70, bounds)


#%%
# fig, ax = plt.subplots()
# ax.scatter(axis, count1, s = 15, label = "Solvent")
# ax.scatter(axis, count2, s = 15, label = "HB - Bilayer")
# ax.scatter(axis, count3, s = 15, label = "TB - Bilayer")
# ax.scatter(axis, count4, s = 15, label = "HB - Micelle")
# ax.scatter(axis, count5, s = 15, label = "TB - Micelle")

# plt.figtext(.5,1.02, "Number Density along the Z-axis", 
#             fontsize=14, ha='center')

# plt.figtext(.5,.96,
#             f"Bilayer - {bilayer_struct}, Micelle - {micelle_struct}", 
#             fontsize=11, ha='center')

# plt.figtext(.5,.91,
#             "Modified coefficient: " + varied, 
#             fontsize=11, ha='center')

# ax.set_xlabel("Distance $x^*$", fontsize = 12)
# ax.set_ylabel(r"Number Density $\rho^*$", fontsize = 12)
# ax.legend(bbox_to_anchor = (1, 1))

# plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

# ax.set_xlim([0, bounds[2]])
# ax.set_ylim([0, 4])    

# fig.savefig(f'{filename}_lastslice_hist.png', format='png', dpi=1200, 
#              bbox_inches='tight')

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

#%%
# fig, ax = plt.subplots()
# step = 20000
# timesteps = np.arange(0, num_slices) * step/1000

# ax.scatter(timesteps, stuck4, s = 15, label = "HB - Micelle")
# ax.scatter(timesteps, stuck5, s = 15, label = "TB - Micelle")

# plt.figtext(.5,.97, f"Number of Beads near Bilayer over Time", 
#             fontsize=14, ha='center')

# plt.figtext(.5,.91,
#             f"For {varied}",
#             fontsize=11, ha='center')

# plt.legend(bbox_to_anchor = (1.35, 0.6))
# ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
# ax.set_ylabel(r"Count", fontsize = 12)
# plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

# ax.set_xlim([0, timesteps.max() * 1.1])
# ax.set_ylim([0, max(stuck4.max(), stuck5.max()) * 1.1])    

# fig.savefig(f'{filename}_beads_near_bilayer.png', format='png', dpi=1200, 
#              bbox_inches='tight')

#%%
# num4 = len(coord4)
# num5 = len(coord5)

# fig, ax = plt.subplots()
# step = 20000
# timesteps = np.arange(0, num_slices) * step/1000

# ax.scatter(timesteps, stuck4/num4, s = 15, label = "HB - Micelle")
# ax.scatter(timesteps, stuck5/num4, s = 15, label = "TB - Micelle")

# plt.figtext(.5,.97, f"Fraction of Beads near Bilayer over Time", 
#             fontsize=14, ha='center')

# plt.figtext(.5,.91,
#             f"For {varied}",
#             fontsize=11, ha='center')

# plt.legend(bbox_to_anchor = (1.35, 0.6))
# ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
# ax.set_ylabel(r"Count", fontsize = 12)
# plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

# ax.set_xlim([0, timesteps.max() * 1.1])
# ax.set_ylim([0, max((stuck4/num4).max(), (stuck5/num5).max()) * 1.1])    

# fig.savefig(f'{filename}_frac_beads_near_bilayer.png', format='png', dpi=1200, 
#              bbox_inches='tight')