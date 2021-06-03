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

a_ij = [float(x.split("_")[-1]) for x in subfolders]

a_ij_values = {}
for i in range(len(a_ij)):
    a_ij_values[a_ij[i]] = subfolders[i]
    
#%% Load datafiles - Average quantities
avpressure = {}
avtemp = {}
press_xz = {}
density = {}
timestep_av_values = {}
timestep_press = {}

for key, entry in a_ij_values.items():
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

for key, entry in a_ij_values.items():
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
for key, entry in a_ij_values.items():
    x = vprofile_vx[key][45, :, 1]
    y = np.average(vprofile_vx[key][:, :, 3], axis = 0)
    ax.scatter(x, y, s = 15, label = "$a_{solvent} = " + str(key) + "$")
    ax.plot(x, y, linestyle = "--")
    # ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("$z^*$", fontsize = 12)
ax.set_ylabel("$v_x^*$", fontsize = 12)
fig.legend(loc = "center right")
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 10])  
    
fig.savefig("velo_profile_aij_vary.png", format='png', dpi=1200, bbox_inches='tight')  

plt.show()

#%% Get viscosities, assuming couette flow in x-direction
shearrate_star = {}
viscosity_star = {}
for key, entry in a_ij_values.items():
    shearrate_star[key] = np.zeros(len(vprofile_vx[key]))
    
    x = vprofile_vx[key][i, :, 1]
    y = np.average(vprofile_vx[key][:, :, 3], axis = 0)
    shearrate_star[key] = np.polyfit(x, y, 1)[0]
    
    y_ave = np.average(press_xz[key])
    
    viscosity_star[key] = -y_ave/shearrate_star[key]

#%%
fig, ax = plt.subplots()
aij_val = []
stress_xz = []

for key, entry in a_ij_values.items():
    aij_val += [key]
    stress_xz += [-np.average(press_xz[key])]

plt.title("Bulk Stress $\sigma_{xz}$ against Repulsion Coefficient", fontsize = 14)
ax.scatter(aij_val, stress_xz, s = 20)
ax.plot(aij_val, stress_xz, ls = "--", linewidth = 1)
ax.set_xlabel("$a_{solvent}$", fontsize = 12)
ax.set_ylabel("$\sigma_{xz}^*$", fontsize = 12)

#ax.set_xlim([0.5, 3]) 
#ax.set_ylim([0, 3])
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
fig.savefig("stress_xz_aij_vary.png", format='png', dpi=1200, bbox_inches='tight')
    
#%% Performing conversion to real units
N_A = 6.02214 * 10**(23)
k_b = 1.38064852 * 10**(-23)

# Volume Properties of water
M_H2O = 18.01528 * 10**(-3)
V_H2O = 30 * 10**(-30)

# Simulation properties
rho_starr = {}
for key, entry in a_ij_values.items():
    rho_starr[key] = 3
    
dt_star = 0.04
r_c_star = 1
eta_star = 1
m_star = 1
temp_star = 1
temp = 298 # Predetermined mapping of 1 : 298K for temperature

# CG factor of 3 H2O molecules to 1 bead
N_m = 20

# Cutoff r_c_star was defined as 1
epsilon = temp * k_b/temp_star

r_c = {}
rho = {}
dt = {}
eta = {}
P = {}
m = m_star * (N_m * M_H2O)/N_A
for key, rho_star in rho_starr.items():
    r_c[key] = r_c_star * (rho_star * N_m * V_H2O)**(1/3)
    rho[key] = rho_star * (m/r_c[key]**3) # Results in grams!
    dt[key] = 1 * (m**0.5 * r_c[key])/((temp * k_b)**0.5)
    eta[key] = eta_star * (dt[key] * temp * k_b)/(r_c[key]**3 * dt_star)
    P[key] = epsilon/r_c[key]**3

#%% Get actual viscosities and plot
data = np.genfromtxt("temp_viscosity_data.txt")

fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 2] * 1000, c = "k", linewidth = 1.5,
        label = "Experimental Viscosity")


t_simulation  = {}

viscosity = {}
for key, entry in a_ij_values.items():
    viscosity[key] = eta[key] * viscosity_star[key]
    x_ave, y_ave = movingaverage(timestep_av_values[key], avtemp[key])
    t_simulation[key] = y_ave[-1] * (epsilon/k_b) - 273.15
    ax.scatter(t_simulation[key], viscosity[key] * 1000, s = 30,
               label = "$a_{solvent} = " + str(key) + "$") 


plt.title("Dynamic Viscosity against Temperature", fontsize = 14)
ax.set_xlabel("T ($\circ C$)", fontsize = 12)
ax.set_ylabel("$\eta\ (mPa\ s)$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)


ax.set_xlim([40, 50])  
ax.set_ylim([0.4, 0.9])  
fig.legend(bbox_to_anchor = (1.3, 0.7))

fig.savefig("viscosity_aij_vary.png", format='png', dpi=1200, bbox_inches='tight')  








