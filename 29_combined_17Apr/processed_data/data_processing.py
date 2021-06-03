import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import re

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

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
bounds = np.array([36, 36, 36 * 3])
filenames = os.listdir()
filenames4 = [x for x in filenames if "stuck4_" in x]
filenames5 = [x for x in filenames if "stuck5_" in x]

val = [int(x.split("_")[-1].split(".")[0]) for x in filenames4]
for i in range(len(val)):
    if val[0] > val[-1]:
        val = val[1:] + [val[0]]
        filenames4 = filenames4[1:] + [filenames4[0]]
        filenames5 = filenames5[1:] + [filenames5[0]]

bilayer_struct = "$H_5 T_4 T_4$"
micelle_struct = "$H_6 T_4$"

num4 = np.load("coord4.npy")
num4 = len(num4)
num5 = np.load("coord5.npy")
num5 = len(num5)

#%%
# coord1 = np.load("coord1.npy")
# coord2 = np.load("coord2.npy")
# coord3 = np.load("coord3.npy")
# coord4 = np.load("coord4.npy")
# coord5 = np.load("coord5.npy")

# index = -1
# snap1 = coord1[:, :, index]
# snap2 = coord2[:, :, index]
# snap3 = coord3[:, :, index]
# snap4 = coord4[:, :, index]
# snap5 = coord5[:, :, index]

# count1, axis = get_hist(snap1, 70, bounds)
# count2, axis = get_hist(snap2, 70, bounds)
# count3, axis = get_hist(snap3, 70, bounds)
# count4, axis = get_hist(snap4, 70, bounds)
# count5, axis = get_hist(snap5, 70, bounds)

#%% Just to check if the bounds set (19, 80) for bilayer thresholds are valid!
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

# ax.set_xlabel("Distance $x^*$", fontsize = 12)
# ax.set_ylabel(r"Number Density $\rho^*$", fontsize = 12)
# ax.legend(bbox_to_anchor = (1, 1))

# plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

# ax.set_xlim([0, bounds[2]])
# ax.set_ylim([0, 4])    

#%%

for i in range(len(val)):
    fig, ax = plt.subplots()
    varied = f"$x_1={val[i]}$"
    filename = f"x1_{val[i]}"

    stuck4 = np.load(filenames4[i])
    stuck5 = np.load(filenames5[i])
    num_slices = len(stuck4)
    
    step = 20000
    timesteps = np.arange(0, num_slices) * step/1000

    ax.scatter(timesteps, stuck4, s = 15, label = "HB - Micelle")
    ax.scatter(timesteps, stuck5, s = 15, label = "TB - Micelle")
    
    plt.figtext(.5,.97, f"Number of Beads near Bilayer over Time", 
                fontsize=14, ha='center')
    
    plt.figtext(.5,.91,
                f"For {varied}",
                fontsize=11, ha='center')
    
    plt.legend(bbox_to_anchor = (1.35, 0.6))
    ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
    ax.set_ylabel(r"Count", fontsize = 12)
    plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
    
    ax.set_xlim([0, timesteps.max() * 1.1])
    ax.set_ylim([0, max(stuck4.max(), stuck5.max()) * 1.1])    
    
    fig.savefig(f'{filename}_beads_near_bilayer.png', format='png', dpi=1200, 
                  bbox_inches='tight')
    
    plt.clf()
    plt.close(fig)
    
    fig, ax = plt.subplots()
    ax.scatter(timesteps, stuck4/num4, s = 15, label = "HB - Micelle")
    ax.scatter(timesteps, stuck5/num5, s = 15, label = "TB - Micelle")
    
    plt.figtext(.5,.97, f"Fraction of Beads near Bilayer over Time", 
                fontsize=14, ha='center')
    
    plt.figtext(.5,.91,
                f"For {varied}",
                fontsize=11, ha='center')
    
    plt.legend(bbox_to_anchor = (1.35, 0.6))
    ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
    ax.set_ylabel(r"Count", fontsize = 12)
    plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
    
    ax.set_xlim([0, timesteps.max() * 1.1])
    ax.set_ylim([0, max((stuck4/num4).max(), (stuck5/num5).max()) * 1.1])    
    
    fig.savefig(f'{filename}_frac_beads_near_bilayer.png', format='png', dpi=1200, 
                  bbox_inches='tight')  
    
    plt.clf()
    plt.close(fig)
    
#%%
fig, ax = plt.subplots()

color_array = np.linspace(0, 255, num = len(val)).astype(int)
cmap = cm.nipy_spectral(color_array)

for i in range(len(val)):
    if i == 1:
        continue
    if i == 3:
        continue
    if i == 5:
        continue
    
    varied = f"$x_1={val[i]}$"
    filename = f"x1_{val[i]}"

    stuck4 = np.load(filenames4[i])
    num_slices = len(stuck4)
    
    step = 20000
    timesteps = np.arange(0, num_slices) * step/1000

    ax.scatter(timesteps, stuck4/num4, s = 12, label = varied, color = cmap[i])
    ax.plot(timesteps, stuck4/num4, ls = "--", color = cmap[i], lw = 1.2)
    
plt.figtext(.5,.95, f"Number of Micelle Head Beads near Bilayer over Time", 
            fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1.32, 1))
ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
ax.set_ylabel(r"Count", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 1300])
ax.set_ylim([0, 1])    

fig.savefig(f'all_head_beads_near_bilayer.png', format='png', dpi=1200, 
              bbox_inches='tight')

#%%
fig, ax = plt.subplots()

color_array = np.linspace(0, 255, num = len(val)).astype(int)
cmap = cm.nipy_spectral(color_array)

for i in range(len(val)):
    if i == 1:
        continue
    if i == 3:
        continue
    if i == 5:
        continue
    
    varied = f"$x_1={val[i]}$"
    filename = f"x1_{val[i]}"

    stuck5 = np.load(filenames5[i])
    num_slices = len(stuck5)
    
    step = 20000
    timesteps = np.arange(0, num_slices) * step/1000

    ax.scatter(timesteps, stuck5/num5, s = 12, label = varied, color = cmap[i])
    ax.plot(timesteps, stuck5/num5, ls = "--", color = cmap[i], lw = 1.2)
    
plt.figtext(.5,.95, f"Number of Micelle Tail Beads near Bilayer over Time", 
            fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1.32, 1))
ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
ax.set_ylabel(r"Count", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 1300])
ax.set_ylim([0, 1])    

fig.savefig(f'all_tail_beads_near_bilayer.png', format='png', dpi=1200, 
              bbox_inches='tight')


