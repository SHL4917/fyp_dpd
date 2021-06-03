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
struct = "h5t4t4"
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
            f"Polymer Structure: $H_5 T_4 T_4$", 
            fontsize=11, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $P^*(z)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('surface_tension_vs_area_h5t4t4.png', format='png', dpi=1200, 
            bbox_inches='tight')

#%%
np.savetxt(f"preff_area_vs_tension_{struct}.txt", 
           np.concatenate(([preff_area], [tension]), axis = 0).T)

#%%
fig, ax = plt.subplots()

for i in range(len(conc)):
    if i != 6:
        continue
    
    p_profile = np.load(f"p_profile_conc_{conc[i]}.npy") 
    p_axis = np.load(f"p_profile_axis_conc_{conc[i]}.npy")
    
    hist1 = np.load(f"hist1_conc_{conc[i]}.npy")
    hist2 = np.load(f"hist2_conc_{conc[i]}.npy")
    hist3 = np.load(f"hist3_conc_{conc[i]}.npy")
    hist_axis = np.load(f"hist_axis_conc_{conc[i]}.npy")

par1 = ax.twinx()

dot1 = ax.scatter(p_axis, p_profile, s = 50, color = "k", 
           label = "Bilayer Stress Profile")
ax.plot(p_axis, p_profile, linestyle = "--", color = "k")
dot2 = par1.scatter(hist_axis, hist1, s = 15, color = "tab:blue",
             label = "Water Beads")
dot3 = par1.scatter(hist_axis, hist2, s = 15, color = "tab:orange",
             label = "Hydrophilic Beads")
dot4 = par1.scatter(hist_axis, hist3, s = 15, color = "tab:green",
             label = "Hydrophobic Beads")

ax.set_xlim([7, 20])
par1.set_ylim([0, 5])

plt.figtext(.5,.97, "Surface Tension Overlayed with Density Profile", 
            fontsize=14, ha='center')

plt.figtext(.5,.91,
            f"Polymer Structure: $H_5 T_4 T_4$", 
            fontsize=11, ha='center')

ax.set_xlabel("Distance $x^*$", fontsize = 12)
ax.set_ylabel(r"Pressure Profile $P^*(z)$", fontsize = 12)

par1.set_ylabel(r"Density $\rho^*$", fontsize = 12)

dots = [dot1, dot2, dot3, dot4]
plt.legend(dots, [d.get_label() for d in dots], bbox_to_anchor = (1.15, 0.62))

plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
fig.savefig('overlaid_surface_tension.png', format='png', dpi=1200, 
            bbox_inches='tight')
