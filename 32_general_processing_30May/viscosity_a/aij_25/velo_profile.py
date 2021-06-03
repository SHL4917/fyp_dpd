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

#%%
file_name = "velo_profile_vx.txt"
box_length = 10
data = np.genfromtxt(file_name, invalid_raise = False)
v_vx = np.genfromtxt(file_name, invalid_raise = False, skip_header = 4)
    
num_bins = data[0, 1].astype(int)
num_entries = (len(v_vx)/num_bins).astype(int)
    
v_vx = np.reshape(v_vx, (num_entries, num_bins, 4))
v_vx[:, :, 1] = v_vx[:, :, 1] * box_length

#%% Plot velocity profile

fig, ax = plt.subplots()
x = v_vx[45, :, 1]
y = np.average(v_vx[:, :, 3], axis = 0)
ax.scatter(x, y, s = 15)
ax.plot(x, y, linestyle = "--")

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("$z^*$", fontsize = 12)
ax.set_ylabel("$v_x^*$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 10])  
    
fig.savefig("velo_profile_ideal.png", format='png', dpi=1200, bbox_inches='tight')  

plt.show()