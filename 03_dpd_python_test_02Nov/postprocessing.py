#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% Load datafiles - Average quantities
data = np.loadtxt("density.txt")
timestep = data[:, 0]
avdensity = data[:, 1]

#%% Load datafiles - msd
data = np.loadtxt("msd_vec.txt")
msd = data[:, 1]

#%% Plot density w.r.t. timestep
fig, ax = plt.subplots()
ax.scatter(timestep, avdensity, marker = ".")
ax.set_ylim([0, 1])

plt.title("Reduced Density against Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel(r"$\rho^* $", fontsize = 12)
plt.show()

#%% Plot msd w.r.t. timestep
fig, ax = plt.subplots()
ax.scatter(timestep, msd, marker = ".")
ax.plot(np.unique(timestep), 
        np.poly1d(np.polyfit(timestep, msd, 1))(np.unique(timestep)), 
        color = "k")

plt.title("MSD against Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$msd^*$", fontsize = 12)

fig.legend(handles = (ax.get_children() ), 
           labels = ('Raw Data', 'Linear Fit'), loc = "center right")

plt.show()

#%% Get diffusion coefficient
diffusion_star = np.polyfit(timestep, msd, 1)[0]

