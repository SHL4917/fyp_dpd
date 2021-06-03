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
struct = "h5t5t5"

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
            f"Polymer Structure: $H_5 T_5 T_5$", 
            fontsize=11, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('surface_tension_vs_area_h5t5t5.png', format='png', dpi=1200, 
            bbox_inches='tight')

#%%
np.savetxt(f"preff_area_vs_tension_{struct}.txt", 
           np.concatenate(([preff_area], [tension]), axis = 0).T)