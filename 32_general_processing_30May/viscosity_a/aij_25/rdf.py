#%% Imports and setup
import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import math

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%%
# dump = "dump.atom"

# timesteps, coord = get_data_from_dump(dump, num_particles = 1, 
#                                       num_variables = 0)
# np.save("coords.npy", coord[1])

#%%
rho = 3
coords = np.load("coords.npy")[:, 0:3, :]
N = len(coords)
bounds = 10

#%%
to_add = coords.copy()   
for i in range(3):
    for j in range(3):
        for k in range(3):
            copied = to_add.copy()
            copied[:, 0] = copied[:, 0] + i * bounds
            copied[:, 1] = copied[:, 1] + j * bounds
            copied[:, 2] = copied[:, 2] + k * bounds
            coords = np.concatenate((coords, copied), axis = 0)

coords = coords[np.shape(to_add)[0]:, :]
coords[:, 0] = coords[:, 0] - bounds
coords[:, 1] = coords[:, 1] - bounds
coords[:, 2] = coords[:, 2] - bounds
    
#%%
num_bins = 150
limit = 3
histogram = np.linspace(0.0000001, limit, num_bins)
denom = N * rho * 4 * np.pi * histogram**2 * (histogram[2] - histogram[1])

#%%
rdf_ave = np.zeros(num_bins)

for i in range(40):
    rdf = np.zeros(num_bins)
    coord_slice = coords[:, :, i * 5]
    
    for j in range(len(coord_slice)):
        rdf_slice = np.zeros(num_bins)
        r_j = coord_slice[j, :]
    
        if np.any(r_j < 0) or np.any(r_j > 10):
            continue
        
        dist = np.linalg.norm(coord_slice[np.arange(len(coord_slice))!=j , :] - r_j, axis = 1)
        dist = dist[dist < limit]
        
        int_dist = ((dist/limit) * num_bins).astype(int)
        
        for k in range(len(int_dist)):
            rdf_slice[int_dist[k]] = rdf_slice[int_dist[k]] + 1
        
        rdf += rdf_slice
    
    rdf = rdf/denom
    rdf_ave += rdf
    
    print(i)

rdf_ave = rdf_ave/40

#%%
fig, ax = plt.subplots()
ax.scatter(histogram, rdf_ave, s = 10)
ax.plot(histogram, rdf_ave, ls = "--", lw = "0.5")

plt.figtext(.5,.93, f"Radial Distribution Function", 
            fontsize=14, ha='center')

ax.set_xlabel("$r^*$", fontsize = 12)
ax.set_ylabel("$g(r)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, limit])
ax.set_ylim([0, 1.5])    

fig.savefig(f'rdf_ideal.png', format='png', dpi=1200, 
              bbox_inches='tight')