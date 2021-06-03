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
filenames = os.listdir()
filenames = [x for x in filenames if "e2e" in x]

poly_lengths = [int(x.split("_")[-1]) for x in filenames]
poly_id_filename = [f"length_{x}.npy" for x in poly_lengths]

#%%
e2e = {}
for i in range(len(filenames)):
    dump = base + f"\{filenames[i]}\dump.poly"
    polyfile = base + f"\{filenames[i]}\{poly_id_filename[i]}"
    polymer_coord = get_poly_from_dump(dump, polyfile)
    
    e2e_slice = np.zeros(np.shape(polymer_coord)[3])
    for j in range(np.shape(polymer_coord)[3]):
        e2e_slice[j] = np.average(get_e2e(polymer_coord[:, :, :,j]))
    
    e2e[poly_lengths[i]] = np.average(e2e_slice)

#%%
fig, ax = plt.subplots()

x = np.zeros(len(e2e))
y = np.zeros(len(e2e))

i = 0
for key, item in e2e.items():
    x[i] = key
    y[i] = item
    i = i + 1

x = np.sort(x)
y = np.sort(y)

coeffs = np.polyfit(np.log(x[2:]), np.log(y[2:]), deg = 1)
linex = np.linspace(np.log(x[2]), np.log(x[-1]))
liney = np.polyval(coeffs, linex)

ax.scatter(np.log(x[2:]), np.log(y[2:]), s = 20, label = "Points")
ax.plot(linex, liney, ls = "--", lw = 1, 
        label = f"Best-fit Line of Gradient = {coeffs[0]:.3f}")

plt.figtext(.5,.93, f"Plot of Polymer End-to-End Distance against Number of Beads", 
            fontsize=14, ha='center')
plt.legend()
ax.set_xlabel("$ln(N)$", fontsize = 12)
ax.set_ylabel(r"$ln\left( \sqrt{\langle R^2 \rangle}\right)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1) 

fig.savefig(f'e2e_vs_num.png', format='png', dpi=1200, bbox_inches='tight')