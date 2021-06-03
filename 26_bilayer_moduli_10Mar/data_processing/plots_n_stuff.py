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
struct_filename = struct_filename[1:] + [struct_filename[0]]

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

fig.savefig('surface_tension_vs_area_all.png', format='png', dpi=1200, 
            bbox_inches='tight')

#%%
fig, ax = plt.subplots()

data = np.loadtxt(struct_filename[0])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 0
ind = 15
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[0])
ax.plot(preff_area[beg:ind], tension[beg:ind], 
        linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[1])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 2
ind = 12
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[1])
ax.plot(preff_area[beg:ind], tension[beg:ind], 
        linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[2])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 1
ind = 15
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[2])
ax.plot(preff_area[beg:ind], tension[beg:ind], 
        linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[3])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 1
ind = 11
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[3])
ax.plot(preff_area[beg:ind], tension[beg:ind], 
        linestyle = "--", linewidth = 1)

data = np.loadtxt(struct_filename[4])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 0
ind = 11
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[4])
ax.plot(preff_area[beg:ind], tension[beg:ind], 
        linestyle = "--", linewidth = 1)

ax.axhline(0, color = "k", linestyle = "--", linewidth = 0.5)

plt.figtext(.5,.92, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

ax.set_xlabel("Preferred Area $A^*$", fontsize = 12)
ax.set_ylabel(r"Surface Tension $\sigma_z^*$", fontsize = 12)
ax.legend(bbox_to_anchor = (1.1, 1))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig('surface_tension_vs_area_trimmed.png', format='png', dpi=1200, 
            bbox_inches='tight')

#%%
proj_area = np.zeros(len(struct))
elastic_mod = np.zeros(len(struct))
fig, ax = plt.subplots()

data = np.loadtxt(struct_filename[0])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 0
ind = 15
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[0])
coeff = np.polyfit(preff_area[beg:ind], tension[beg:ind], 1)
x = np.linspace(1.4, 1.6, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[0] = -coeff[1]/coeff[0]
elastic_mod[0] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[1])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 2
ind = 12
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[1])
coeff = np.polyfit(preff_area[beg:ind], tension[beg:ind], 1)
x = np.linspace(1.4, 1.68, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[1] = -coeff[1]/coeff[0]
elastic_mod[1] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[2])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 1
ind = 15
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[2])
coeff = np.polyfit(preff_area[beg:ind], tension[beg:ind], 1)
x = np.linspace(1.4, 1.68, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[2] = -coeff[1]/coeff[0]
elastic_mod[2] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[3])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 1
ind = 11
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[3])
coeff = np.polyfit(preff_area[beg:ind], tension[beg:ind], 1)
x = np.linspace(1.4, 1.6, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[3] = -coeff[1]/coeff[0]
elastic_mod[3] = proj_area[0] * coeff[0]

data = np.loadtxt(struct_filename[4])
preff_area = data[:, 0]
tension = data[:, 1] * -1

beg = 0
ind = 11
ax.scatter(preff_area[beg:ind], tension[beg:ind], s = 15, label = struct[4])
coeff = np.polyfit(preff_area[beg:ind], tension[beg:ind], 1)
x = np.linspace(1.45, 1.55, 100)
y = np.polyval(coeff, x)
ax.plot(x, y, linestyle = "--", linewidth = 1)

proj_area[4] = -coeff[1]/coeff[0]
elastic_mod[4] = proj_area[0] * coeff[0]

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
bead_num = np.array([3, 5, 6, 9, 11])

fig, ax = plt.subplots()

ax.scatter(bead_num, elastic_mod, s = 15)
ax.plot(bead_num, elastic_mod, linestyle = "--", linewidth = 1)

ax.set_xlim([2, 12])
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

