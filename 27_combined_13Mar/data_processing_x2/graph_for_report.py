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
values = [15, 40, 70, 100]
filenames_4 = [f"stuck4_$x_2={x}$.npy" for x in values]
filenames_5 = [f"stuck5_$x_2={x}$.npy" for x in values]

#%%
c4 = np.load("coord4.npy")
c5 = np.load("coord5.npy")

tot = len(c4) + len(c5)

#%%
fig, ax = plt.subplots()

for i in range(len(values)):
    stuck_4 = np.load(filenames_4[i])
    stuck_5 = np.load(filenames_5[i])
    stuck_all = stuck_4 + stuck_5
    
    timestep = np.arange(len(stuck_4)) * 20000 * 0.04
    
    ax.scatter(timestep, stuck_all/tot, s = 15, 
               label = "$x_2=" + str(values[i]) + "$")
    ax.plot(timestep, stuck_all/tot, ls = "--", lw = 1)

#plt.figtext(.5,.95, f"Number of Micelle Particles near Bilayer over Time", fontsize=14, ha='center')

plt.legend(loc = "upper left")
ax.set_xlabel("Time $t^*$", fontsize = 12)
ax.set_ylabel(r"Fraction of All Particles", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 15000])
ax.set_ylim([0, 0.05])    

fig.savefig(f'all_beads_near_bilayer_x2.png', format='png', dpi=1200, bbox_inches='tight')

