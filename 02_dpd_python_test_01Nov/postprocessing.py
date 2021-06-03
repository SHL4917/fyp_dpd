#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% Load datafiles - Average quantities
data = np.loadtxt("avpressure.txt")
timestep_av_pres_and_temp = data[:, 0]
avpressure = data[:, 1]

data = np.loadtxt("avtemp.txt")
avtemp = data[:, 1]

data = np.loadtxt("press_xz.txt")
timestep_press = data[:, 0]
press_xz = data[:, 1]

#%% Load datafiles - Velocity profiles
data = np.genfromtxt("velo_profile_vx.txt", invalid_raise = False)
timestep_velo_profile = data[:, 0]
velo_profile_vx = np.genfromtxt("velo_profile_vx.txt", invalid_raise = False, 
                      skip_header = 4)

num_bins = data[0, 1].astype(int)
box_length = 10
num_entries = (len(velo_profile_vx)/num_bins).astype(int)
velo_profile_vx = np.reshape(velo_profile_vx, (num_entries, num_bins, 4))
velo_profile_vx[:, :, 1] = velo_profile_vx[:, :, 1] * box_length

#%% Plotting pressure and temperature
fig, ax = plt.subplots(1, 2)
ax[0].scatter(timestep_av_pres_and_temp, avpressure, marker = ".")
ax[0].set_ylim([0, 1.2])
ax[1].scatter(timestep_av_pres_and_temp, avtemp, marker = ".")
ax[1].set_ylim([0, 1.2])

ax[0].set_title("Reduced Pressure", fontsize = 14)
ax[1].set_title("Reduced Temperature", fontsize = 14)

ax[0].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[0].set_ylabel("P*", fontsize = 12)

ax[1].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[1].set_ylabel("T*", fontsize = 12)
fig.tight_layout()

plt.show()

#%% Plotting velocity profile
fig, ax = plt.subplots()
x = velo_profile_vx[66, :, 1]
y = velo_profile_vx[66, :, 3]
ax.scatter(x, y)
ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), 
        color = "k")

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("V*", fontsize = 12)

fig.legend(handles = (ax.get_children() ), 
           labels = ('Raw Data', 'Linear Fit'), loc = "center right")

plt.show()

#%% Get viscosities, assuming couette flow in x-direction
shearrate_star = np.zeros(len(timestep_velo_profile))

for i in range(len(shearrate_star)):
    x = velo_profile_vx[i, :, 1]
    y = velo_profile_vx[i, :, 3]
    shearrate_star[i] = np.polyfit(x, y, 1)[0]

#%% Plot shear rate w.r.t timestep
fig, ax = plt.subplots()
ax.scatter(timestep_velo_profile, shearrate_star, marker = ".")
ax.set_ylim([0, 0.25])
ax.hlines(np.mean(shearrate_star), 
          timestep_velo_profile[0], timestep_velo_profile[-1])

plt.title("Shear Rate w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$\gamma^*$", fontsize = 12)

plt.show()

# Can see that exceeds original set shear rate of 1.5

#%% Plot pressure in xz direction
fig, ax = plt.subplots()
ax.scatter(timestep_press, press_xz, marker = ".")
ax.set_ylim([-0.05, 0])
ax.hlines(np.mean(press_xz), 
          timestep_press[0], timestep_press[-1])

plt.title("Pressure (XZ component) w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$P_{xz}^*$", fontsize = 12)
plt.show()

#%%
viscosity_star = -np.mean(press_xz)/np.mean(shearrate_star)



