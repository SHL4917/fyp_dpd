import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import shutil
from scipy.signal import savgol_filter
import scipy as sp

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
bound = 36
struct = "h5t3t3"
filenames = os.listdir()
conc = [x.split("_")[-1].split(".")[0] for x in filenames if "hist_axis" in x]
conc = [int(x) for x in conc]
conc = np.array(conc)

poly_count = np.loadtxt("poly_count.txt")

#%%
tension = np.zeros(np.size(poly_count))
preff_area = np.zeros(np.size(poly_count))

fig, ax = plt.subplots()

for i in range(len(conc)):   
    
    p_profile = np.load(f"p_profile_conc_{conc[i]}.npy") 
    p_axis = np.load(f"p_profile_axis_conc_{conc[i]}.npy")
    
    tension[i] = sp.integrate.simps(p_profile, p_axis)
    preff_area = (bound**2) / (poly_count * 0.5)

preff_area, tension = zip(*sorted(zip(preff_area, tension)))
preff_area = np.array(preff_area)
tension = np.array(tension)

ax.scatter(preff_area, tension, s = 20)
ax.plot(preff_area, tension, linestyle = "--")

plt.figtext(.5,.97, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

plt.figtext(.5,.91,
            f"Polymer Structure: $H_5 T_3 T_3$", 
            fontsize=11, ha='center')

ax.set_xlabel("Preferred Area $x^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $P^*(z)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

#fig.savefig('surface_tension_vs_area_h5t3t3.png', format='png', dpi=1200, 
            #bbox_inches='tight')

#%%
np.savetxt(f"preff_area_vs_tension_{struct}.txt", 
           np.concatenate(([preff_area], [tension]), axis = 0).T)

#%%
index = 9
fig, ax = plt.subplots()

for i in range(len(conc)):
    if i != index:
        continue
    
    hist1 = np.load(f"hist1_conc_{conc[i]}.npy")
    hist2 = np.load(f"hist2_conc_{conc[i]}.npy")
    hist3 = np.load(f"hist3_conc_{conc[i]}.npy")
    
    hist_axis = np.load(f"hist_axis_conc_{conc[i]}.npy")
    
    ax.scatter(hist_axis, hist1, s = 15, label = "Water Beads",
               color = "tab:blue")
    ax.scatter(hist_axis, hist2, s = 15, label = "Hydrophilic Beads",
               color = "tab:orange")
    ax.scatter(hist_axis, hist3, s = 15, label = "Hydrophobic Beads",
               color = "tab:green")

plt.figtext(.5,.97, "Density Profile of Lipid Bilayer", 
            fontsize=14, ha='center')

plt.figtext(.5,.91,
            f"Averaged over 125 Timesteps at Concentration of "
            f"{poly_count[index] * 11/(3 * 36**3):.3f}", 
            fontsize=11, ha='center')

plt.legend(bbox_to_anchor = (0.8, 0.6))
ax.set_xlabel("Distance $x^*$", fontsize = 12)
ax.set_ylabel(r"Density $\rho^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([5, 20])
ax.set_ylim([0, 4])    

#fig.savefig('bilayer_density_profile_hole.png', format='png', dpi=1200, 
             #bbox_inches='tight')