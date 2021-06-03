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
bound = 36
filenames = os.listdir()

struct_filename = [x for x in filenames if "preff_area_vs" in x]
struct_filename = struct_filename[1:] + [struct_filename[0]]

struct = [x.split("_")[-1].split(".")[0] for x in struct_filename]

for i in range(len(struct)):
    if len(struct[0]) > 6:
        struct = struct[1:] + [struct[0]]
        struct_filename = struct_filename[1:] + [struct_filename[0]]
        
struct_latex = ["$H_{" + re.split('(\d+)', x)[1] + \
                "} T_4 T_4$" for x in struct]
    
bead_num = [int(re.split('(\d+)', x)[1]) for x in struct]

bead_num = np.array(bead_num)

#%%
fig, ax = plt.subplots()

b1 = np.array([14, 11, 10, 5, 0, 1, 2, 0, 0, 0])
b2 = np.array([25, 17, 25, 10, 0, 10, 10, 9, 0, 3])

for i, file in enumerate(struct_filename):
    data = np.loadtxt(file)
    preff_area = data[:, 0]
    tension = data[:, 1] * -1
    if b1[i] == 0 and b2[i] == 0:
        continue
    
    ax.scatter(preff_area[b1[i]:b2[i]], tension[b1[i]:b2[i]], 
               s = 25, label = struct_latex[i])
    ax.plot(preff_area[b1[i]:b2[i]], tension[b1[i]:b2[i]], 
            linestyle = "--", linewidth = 1.2)

plt.figtext(.5,.92, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Projected Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.05, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('surface_tension_vs_area_trimmed.png', format='png', dpi=1200, 
            bbox_inches='tight')

#%%
proj_area = np.zeros(len(struct))
elastic_mod = np.zeros(len(struct))
fig, ax = plt.subplots()

c1 = np.array([1.35, 1.37, 1.38, 1.39, 0, 1.39, 1.39, 1.38, 0, 1.46])
c2 = np.array([1.47, 1.53, 1.57, 1.58, 0, 1.61, 1.61, 1.62, 0, 1.65])

for i, file in enumerate(struct_filename):
    data = np.loadtxt(file)
    preff_area = data[:, 0]
    tension = data[:, 1] * -1
    
    if b1[i] == 0 and b2[i] == 0:
        continue
    
    ax.scatter(preff_area[b1[i]:b2[i]], tension[b1[i]:b2[i]], 
               s = 25, label = struct_latex[i])
    coeff = np.polyfit(preff_area[b1[i]:b2[i]], tension[b1[i]:b2[i]], 1)
    
    x = np.linspace(c1[i], c2[i], 100)
    y = np.polyval(coeff, x)
    
    ax.plot(x, y, linestyle = "--", linewidth = 1.2)
    
    proj_area[i] = -coeff[1]/coeff[0]
    elastic_mod[i] = proj_area[i] * coeff[0]

ax.axhline(0, color = "k", linestyle = "--", linewidth = 0.5)

plt.figtext(.5,.93, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.1, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('surface_tension_vs_area_best_fit.png', format='png', dpi=1200, 
            bbox_inches='tight')
            
#%%
fig, ax = plt.subplots()

ax.scatter(bead_num[elastic_mod != 0], elastic_mod[elastic_mod != 0], s = 25)
ax.plot(bead_num[elastic_mod != 0], elastic_mod[elastic_mod != 0], 
        linestyle = "--", linewidth = 1.2)

ax.set_xlim([1, 12])
#ax.set_ylim([5, 12])

plt.figtext(.5,.97, "Membrane Elastic Modulus against Head Bead Number", 
            fontsize=14, ha='center')

plt.figtext(.5,.91, "Structure - $H_x T_4 T_4$", 
            fontsize=12, ha='center')

ax.set_xlabel("Number of Head Beads", fontsize = 12)
ax.set_ylabel(r"Elastic Modulus $E^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('elastic_mod_vs_tail_number.png', format='png', dpi=1200, 
            bbox_inches='tight')

