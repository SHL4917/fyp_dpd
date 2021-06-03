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
bound = 36
filenames = os.listdir()

struct_filename = [x for x in filenames if "preff_area_vs" in x]
struct = [x.split("_")[-1].split(".")[0] for x in struct_filename]

#%%
fig, ax = plt.subplots()

for i, file in enumerate(struct_filename):
    data = np.loadtxt(file)
    preff_area = data[:, 0]
    tension = data[:, 1] * -1
    ax.scatter(preff_area, tension, s = 15, label = struct[i])
    ax.plot(preff_area, tension, linestyle = "--", linewidth = 1)

plt.figtext(.5,.92, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.1, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

#fig.savefig('surface_tension_vs_area_all.png', format='png', dpi=1200, 
            #bbox_inches='tight')

#%%
fig, ax = plt.subplots()

data = np.loadtxt(struct_filename[0])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 11
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[0])
ax.plot(preff_area[0:ind], tension[0:ind], linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[1])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 12
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[1])
ax.plot(preff_area[0:ind], tension[0:ind], linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[2])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 20
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[2])
ax.plot(preff_area[0:ind], tension[0:ind], linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[3])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 19
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[3])
ax.plot(preff_area[0:ind], tension[0:ind], linestyle = "--", linewidth = 1)

ax.axhline(0, color = "k", linestyle = "--", linewidth = 0.5)

plt.figtext(.5,.92, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.1, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

#fig.savefig('surface_tension_vs_area_trimmed.png', format='png', dpi=1200, 
            #bbox_inches='tight')

#%%
proj_area = np.zeros(4)
elastic_mod = np.zeros(4)
fig, ax = plt.subplots()

data = np.loadtxt(struct_filename[0])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 11
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[0])
coeff = np.polyfit(preff_area[0:ind], tension[0:ind], 1)
x = np.linspace(1.35, 1.55, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[0] = -coeff[1]/coeff[0]
elastic_mod[0] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[1])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 12
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[1])
coeff = np.polyfit(preff_area[0:ind], tension[0:ind], 1)
x = np.linspace(1.37, 1.68, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[1] = -coeff[1]/coeff[0]
elastic_mod[1] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[2])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 20
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[2])
coeff = np.polyfit(preff_area[0:ind], tension[0:ind], 1)
x = np.linspace(1.37, 1.68, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[2] = -coeff[1]/coeff[0]
elastic_mod[2] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[3])
preff_area = data[:, 0]
tension = data[:, 1] * -1

ind = 19
ax.scatter(preff_area[0:ind], tension[0:ind], s = 15, label = struct[3])
coeff = np.polyfit(preff_area[0:ind], tension[0:ind], 1)
x = np.linspace(1.37, 1.68, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[3] = -coeff[1]/coeff[0]
elastic_mod[3] = proj_area[0] * coeff[0]

ax.axhline(0, color = "k", linestyle = "--", linewidth = 0.5)

plt.figtext(.5,.93, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.1, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

#fig.savefig('surface_tension_vs_area_best_fit.png', format='png', dpi=1200, 
            #bbox_inches='tight')
            
#%%
bead_num = np.array([3, 4, 5, 6])

fig, ax = plt.subplots()

ax.scatter(bead_num, elastic_mod/2, s = 20)
ax.plot(bead_num, elastic_mod/2, linestyle = "--", linewidth = 1)

ax.set_xlim([2, 7])
ax.set_ylim([4, 9])

plt.figtext(.5,.97, "Membrane Elastic Modulus against Tail Bead Number", 
            fontsize=14, ha='center')
plt.figtext(.5,.92, "Surfactant Structure of $H_5 T_x T_x$", 
            fontsize=12, ha='center')

ax.set_xlabel("Number of Tail Beads", fontsize = 12)
ax.set_ylabel(r"Elastic Modulus $E^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('elastic_mod_vs_tail_number.png', format='png', dpi=1200, 
            bbox_inches='tight')

