import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import open3d as o3d

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
bound = 30

poly_com = np.load("poly_com.npy")
plane = np.load("plane.npy").T
plane = plane[:, 0]

coord1 = np.load("coord1.npy")
coord2 = np.load("coord2.npy")
coord3 = np.load("coord3.npy")

#%% Get the averaged histogram values, assuming that the bilayer doesn't move!
for i in range(np.shape(coord1)[2]):
    c1 = coord1[:, :, i]
    c2 = coord2[:, :, i]
    c3 = coord3[:, :, i]

    t1, t2, t3, axis = get_bilayer_histogram_norotate(c1[:, 0:3], c2[:, 0:3], 
                                                       c3[:, 0:3], 1500, 
                                                       bound)
    if i == 0:
        hist1 = t1.copy()
        hist2 = t2.copy()
        hist3 = t3.copy()
        continue   
    hist1 = hist1 + t1
    hist2 = hist2 + t2
    hist3 = hist3 + t3
    
    print(f"Timestep: {i}")
    
hist1 = hist1 / (i + 1)
hist2 = hist2 / (i + 1)
hist3 = hist3 / (i + 1)

#%%
fig, ax = plt.subplots()
ax.scatter(axis, hist1, s = 15, color = "k", label = "Water Beads")
ax.scatter(axis, hist2, s = 15, color = "b", label = "Hydrophilic Beads")
ax.scatter(axis, hist3, s = 15, color = "r", label = "Hydrophobic Beads")

ax.set_xlim([5, 20])
ax.set_ylim([0, 4])

plt.figtext(.5,.97, "Density Profile of Lipid Bilayer", 
            fontsize=14, ha='center')

plt.figtext(.5,.91,
            f"Averaged over {np.shape(coord1)[2]} Timesteps", 
            fontsize=11, ha='center')

plt.legend(bbox_to_anchor = (0.8, 0.6))
ax.set_xlabel("Distance $x^*$", fontsize = 12)
ax.set_ylabel(r"Density $\rho^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

#fig.savefig('bilayer_density_profile.png', format='png', dpi=1200, 
#             bbox_inches='tight')

#%%
num_bins = 300
for i in range(np.shape(coord1)[2]):
    c1 = coord1[:, :, i]
    c2 = coord2[:, :, i]
    c3 = coord3[:, :, i]
    
    t_profile, axis = get_p_profile_norotate(np.concatenate([c1, c2, c3]), 
                                             num_bins, bound)
    
    if i == 0:
        p_profile = t_profile.copy()
        continue
    
    p_profile = p_profile + t_profile
    
    print(f"Timestep: {i}")

p_profile = p_profile/(i + 1)

#%%
fig, ax = plt.subplots()
ax.scatter(axis, p_profile, s = 15)

ax.set_xlim([5, 20])

plt.figtext(.5,.97, "Pressure Profile of Lipid Bilayer", 
            fontsize=14, ha='center')

plt.figtext(.5,.91,
            f"Averaged over {np.shape(coord1)[2]} Timesteps", 
            fontsize=11, ha='center')

ax.set_xlabel("Distance $x^*$", fontsize = 12)
ax.set_ylabel(r"Pressure $P^*(z)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('bilayer_pressure_profile_bad.png', format='png', dpi=1200, 
            bbox_inches='tight')
    
    
