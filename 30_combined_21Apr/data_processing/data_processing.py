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
filenames_num = [x for x in filenames if "particle_number" in x]

val = [int(x.split("_")[-1].split(".")[0]) for x in filenames4]
for i in range(len(val)):
    if val[0] > val[-1]:
        val = val[1:] + [val[0]]
        filenames4 = filenames4[1:] + [filenames4[0]]
        filenames5 = filenames5[1:] + [filenames5[0]]
        filenames_num = filenames_num[1:] + [filenames[0]]

bilayer_struct = "$H_5 T_4 T_4$"
micelle_struct = "$H_6 T_4$"

num4 = len(np.load("coord4.npy"))
num5 = len(np.load("coord5.npy"))

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

#color_array = np.linspace(0, 255, num = len(val)).astype(int)
#cmap = cm.nipy_spectral(color_array)

for i in range(len(val)):
    if i == 1:
        continue
    if i == 3:
        continue
    if i == 5:
        continue
    
    varied = f"$x_1={val[i]/100:.2f}$"
    filename = f"x1_{val[i]}"

    stuck4 = np.load(filenames4[i])
    num_slices = len(stuck4)
    
    step = 20000
    timesteps = np.arange(0, num_slices) * step * 0.04

    ax.scatter(timesteps, stuck4/num4, s = 8, label = varied)
    ax.plot(timesteps, stuck4/num4, ls = "--", lw = 1)
    
# plt.figtext(.5,.95, f"Number of Micelle Beads near Bilayer", fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1, 0.7))
ax.set_xlabel("Time $t^*$", fontsize = 12)
ax.set_ylabel(r"Fraction of all Micelle Particles", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 90000])
ax.set_ylim([0, 1])    

fig.savefig(f'frac_beads_bilayer_all_x1.png', format='png', dpi=1200, bbox_inches='tight')

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
    
    varied = f"$x_1={val[i]/100:.2f}$"
    filename = f"x1_{val[i]}"

    stuck5 = np.load(filenames5[i])
    num_slices = len(stuck5)
    
    step = 20000
    timesteps = np.arange(0, num_slices) * step/1000

    ax.scatter(timesteps, stuck5/num5, s = 12, label = varied, color = cmap[i])
    ax.plot(timesteps, stuck5/num5, ls = "--", color = cmap[i], lw = 1.2)
    
plt.figtext(.5,.95, f"Number of Micelle Tail Beads near Bilayer over Time", 
            fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1.1, 0.65))
ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
ax.set_ylabel(r"Count", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 1600])
ax.set_ylim([0, 1])    

fig.savefig(f'all_tail_beads_near_bilayer.png', format='png', dpi=1200, 
              bbox_inches='tight')


