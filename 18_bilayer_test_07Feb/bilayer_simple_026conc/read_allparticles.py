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
bound = 30
dump = "dump.allparticles"
coord1, coord2, coord3 = get_allparticles_from_dump(dump)

#%%
plane = np.loadtxt("plane_coeff.txt")
snap1 = coord1[:, :, 7]
snap2 = coord2[:, :, 7]
snap3 = coord3[:, :, 7]

c1norm, c2norm, c3norm, axis = get_bilayer_histogram(snap1, snap2, snap3, 500, 
                                               plane, bound)

#%% Plot shit
fig, ax = plt.subplots()
ax.scatter(axis, c1norm, linewidth = 0, s = 10, label = "Water Beads")
ax.scatter(axis, c2norm, linewidth = 0, s = 20, label = "Hydrophilic Beads", 
           color = "k")
ax.scatter(axis, c3norm, linewidth = 0, s = 20, label = "Hydrophobic Beads", 
           color = "r")

ax.set_xlim([-10, 5])
ax.set_ylim([0, 3])

plt.figtext(.5,.91,
            "Density Profile of Lipid Bilayer", 
            fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1, 0.5))
ax.set_xlabel("Distance $x^*$", fontsize = 12)
ax.set_ylabel(r"Density $\rho^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('bilayer_density_profile.png', format='png', dpi=1200, 
            bbox_inches='tight')


