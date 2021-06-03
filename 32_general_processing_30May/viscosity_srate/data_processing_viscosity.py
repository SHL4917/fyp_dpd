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
def movingaverage (x, y, lump = 1):
    entries = len(x)
    num_aves = math.floor(entries/lump)
    
    y_ave = np.zeros(num_aves)
    x_ave = np.zeros(num_aves)
    
    for i in range(num_aves):
        ave = np.sum(y[i * lump: (i+1) * lump]) / lump
        y_ave[i] = (np.sum(y_ave) + ave)/(i + 1)
        x_ave[i] = x[i * lump]
    
    return (x_ave, y_ave)

#%%
base = os.getcwd()
subfolders = [name for name in os.listdir(".") if 
              os.path.isdir(os.path.join(base, name))]

srate = [float(x.split("_")[-1])/100 for x in subfolders]

srate_values = {}
for i in range(len(srate)):
    srate_values[srate[i]] = subfolders[i]
    
#%% Load datafiles - Average quantities
avpressure = {}
avtemp = {}
press_xz = {}
density = {}
timestep_av_values = {}
timestep_press = {}

for key, entry in srate_values.items():
    data = np.loadtxt(entry + "/" + "avpressure.txt")
    timestep_av_values[key] = data[:, 0]
    avpressure[key] = data[:, 1]
    
    data = np.loadtxt(entry + "/" + "avtemp.txt")
    avtemp[key] = data[:, 1]
    
    data = np.loadtxt(entry + "/" + "press_xz.txt")
    timestep_press[key] = data[:, 0]
    press_xz[key] = data[:, 1]

#%% Load datafiles - Velocity profiles
timestep_velo_profile = {}
vprofile_vx = {}
tprofile = {}
rhoprofile = {}

box_length = 10

for key, entry in srate_values.items():
    file_name = entry + "/" + "velo_profile_vx.txt"
    data = np.genfromtxt(file_name, invalid_raise = False)
    timestep_velo_profile[key] = data[:, 0]
    v_vx = np.genfromtxt(file_name, invalid_raise = False, skip_header = 4)
    
    num_bins = data[0, 1].astype(int)
    num_entries = (len(v_vx)/num_bins).astype(int)
    
    v_vx = np.reshape(v_vx, (num_entries, num_bins, 4))
    v_vx[:, :, 1] = v_vx[:, :, 1] * box_length
    vprofile_vx[key] = v_vx
    
    file_name = entry + "/" + "tstar_profile.txt"
    t_profile = np.genfromtxt(file_name, invalid_raise = False, 
                              skip_header = 4)
    t_profile = np.reshape(t_profile, (num_entries, num_bins, 4))
    t_profile[:, :, 1] = t_profile[:, :, 1] * box_length
    tprofile[key] = t_profile
    
    file_name = entry + "/" + "rhostar_profile.txt"
    rho = np.genfromtxt(file_name, invalid_raise = False, 
                              skip_header = 4)
    rho = np.reshape(rho, (num_entries, num_bins, 4))
    rho[:, :, 1] = rho[:, :, 1] * box_length
    rhoprofile[key] = rho

#%% Plotting velocity profile -vx
fig, ax = plt.subplots()
for key, entry in srate_values.items():
    x = vprofile_vx[key][45, :, 1]
    y = np.average(vprofile_vx[key][:, :, 3], axis = 0)
    ax.scatter(x, y, s = 15, label = "$\dot{\gamma}_{shear}^* = " + str(key) + "$")
    ax.plot(x, y, linestyle = "--")
    # ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("$z^*$", fontsize = 12)
ax.set_ylabel("$v_x^*$", fontsize = 12)
fig.legend(bbox_to_anchor = (1.15, 0.7))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 10])  
    
fig.savefig("velo_profile_srate_vary.png", format='png', dpi=1200, bbox_inches='tight')  

plt.show()

#%%
fig, ax = plt.subplots()
for key, entry in srate_values.items():
    x = tprofile[key][45, :, 1]
    y = tprofile[key][45, :, 3]
    ax.scatter(x, y, s = 15, label = "$\dot{\gamma}_{shear}^* = " + str(key) + "$")
    ax.plot(x, y, ls = "--", lw = 1)
    
plt.title("Temperature Profile", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("$T^*$", fontsize = 12)
ax.set_ylim([0.5, 3])
ax.set_xlim([0, 10])
fig.legend(bbox_to_anchor = (1.15, 0.7))
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

fig.savefig("temp_profile_srate_vary.png", format='png', dpi=1200, bbox_inches='tight')  
plt.show()






