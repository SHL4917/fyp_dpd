import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import re

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
dump = "dump.poly"
polyfile = "length_15.npy"
polymer_coord = get_poly_from_dump(dump, polyfile)

#%%
e2e = np.zeros(1)

for i in range(np.shape(polymer_coord)[3]):
    e2e_slice = get_e2e(polymer_coord[:, :, :, i])
    e2e = np.concatenate((e2e, e2e_slice))

e2e = e2e[1:]

#%%
e2e_rms = (np.average(e2e**2))**0.5
e2e_normed = e2e/e2e_rms
num_bins = 100

x = np.linspace(0, 3, num_bins)
y =  4 * np.pi * x**2 * ((1/(2 * np.pi))**1.5) * np.exp(-1.5 * x**2)

#%%
fig, ax = plt.subplots()

hist, bins = np.histogram(e2e_normed, bins = num_bins)
bins = bins - (0.5 * (bins[1] - bins[0]))
bins = bins[1:]


ax.scatter(bins, hist/(np.sum(hist * bins)), s = 15, label = "Simulation Data")
ax.plot(x, y/(np.sum(y * x)), lw = 1, label = "Ideal Distribution")

plt.figtext(.5,.93, f"Plot of Polymer End-to-End Distance Distribution", 
            fontsize=14, ha='center')
plt.legend()
ax.set_xlabel(r"$R / \sqrt{\langle R^2 \rangle}$", fontsize = 12)
ax.set_ylabel(r"$P(x)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1) 

ax.set_xlim([0, 2.5])
ax.set_ylim([0, 0.035])    

fig.savefig(f'e2e_dist.png', format='png', dpi=1200, bbox_inches='tight')
    