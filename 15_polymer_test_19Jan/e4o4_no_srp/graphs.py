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

#%% Load data
poly_data = np.load("poly_coords.npy")
com = np.load("com.npy")

#%% Set Slicing Parameters
bounds = 40
step = 550

#%% Get cluster dictionaries
coords = com[:, :, step]
poly_coords = poly_data[:, :, :, step]

cluster_dict, poly_dict = get_clusters(coords, poly_coords, bounds)

#%% Plot number of clusters per time step
timesteps = 600
nclusters = np.zeros(timesteps)
steps = np.arange(0, timesteps, 1)
for i in range(timesteps):
    coords = com[:, :, i]
    poly_coords = poly_data[:, :, :, i]
    
    cluster_dict, poly_dict = get_clusters(coords, poly_coords, bounds)
    print("Timestep " + str(i))
    
    nclusters[i] = len(cluster_dict)
    
#%% Plotting
fig, ax = plt.subplots()
ax.scatter(steps * 1000, nclusters, s = 5, marker = "o", linewidth = 0)

plt.figtext(.5,.95,
            "Number of Clusters Detected Against Timestep", 
            fontsize=14, ha='center')
plt.figtext(.5,.90,
            "Snapshot every 1000 Timesteps", 
            fontsize=12, ha='center')

ax.set_xlabel("Timestep ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("Clusters", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
ax.set_xlim([0, 600000])
ax.set_ylim([0, 60])

fig.savefig('clusters_600steps.png', format='png', dpi=1200, 
            bbox_inches='tight')