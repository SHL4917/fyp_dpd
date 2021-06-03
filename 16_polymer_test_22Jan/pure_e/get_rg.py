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

#%% Iterating over the different folders
r_g_array = np.zeros(6)
mass_array = np.zeros(6)
j = 0

for i in range(10, 36, 5):
    print(i)
    folder = "e" + str(i) + "/"
    polymer_struct = np.zeros(i, dtype = int) + 2
    
    poly_data = np.load(folder + "polymer_coords.npy")
    com = get_com_all(poly_data, polymer_struct)
    np.save(folder + "com.npy", com)
    
    r_g_ave, r_g_std = get_r_g_stats(get_r_g_all(poly_data, com))
    r_g_array[j] = np.average(r_g_ave)
    mass_array[j] = i * 44
    
    j = j + 1

#%% Plot 
fig, ax = plt.subplots()
ax.scatter(np.log(mass_array), np.log(r_g_array), marker = ".", 
           label = "Points")

c = np.polyfit(np.log(mass_array[:]), np.log(r_g_array[:]), deg = 1)
x_vals = np.array(ax.get_xlim())
y_vals = c[0] * x_vals + c[1]

ax.plot(x_vals, y_vals, ls = "--", label = f"Gradient = {c[0]:.3f}")

plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
plt.figtext(.5,.95,
            "Radius of Gyration against Polymer Molar Mass", 
            fontsize=14, ha='center')
plt.figtext(.5,.90,
            "For a 3:1 CG Scheme (Xia and Zhong, 2006)", 
            fontsize=12, ha='center')
ax.set_ylabel("$ln(R_g)$")
ax.set_xlabel("$ln(M)$")
plt.legend(loc = "lower right")

fig.savefig('r_g_3to1.png', format='png', dpi=1200, bbox_inches='tight')