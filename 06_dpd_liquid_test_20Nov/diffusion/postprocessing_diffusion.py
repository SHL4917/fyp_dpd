#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% Get folder names
base = os.getcwd()
subfolders = [name for name in os.listdir(".") if 
              os.path.isdir(os.path.join(base, name))]

a_ij_values = [float(x.split("_")[-1]) for x in subfolders]

#%% Load datafiles - MSD
msd = {}
timestep = {}
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    data = np.loadtxt(subfolders[i] + "/" + "msd_vec.txt")
    msd[entry] = data[:, 1]
    timestep[entry] = data[:, 0]
    
#%% Load datafiles - Density
density = {}
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    data = np.loadtxt(subfolders[i] + "/" + "density.txt")
    density[entry] = data[:, 1]
    
#%% Plot density w.r.t. timestep
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    ax.scatter(timestep[entry], density[entry], marker = ".")

plt.title("Density against Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel(r"$\rho^*$", fontsize = 12)
ax.set_xlim([0, 300000])

fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")
    
#%% Plot msd w.r.t. timestep
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    ax.scatter(timestep[entry] * 0.04, msd[entry], marker = ".")

plt.title("MSD against Timestep", fontsize = 14)
ax.set_xlabel("Dimensionless Time $t^*$", fontsize = 12)
ax.set_ylabel("$msd^*$", fontsize = 12)
ax.set_xlim([0, 10000])
ax.set_ylim([0, 1500])

fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")

plt.show()

#%% Get diffusion coefficient
diffusion_star = np.zeros(len(subfolders))
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    diffusion_star[i] = np.polyfit(timestep[entry] * 0.04, msd[entry], 1)[0]
    
#%% Performing conversion to real units
N_A = 6.02214 * 10**(23)
k_b = 1.38064852 * 10**(-23)

# Volume Properties of water
M_H2O = 18.01528 * 10**(-3)
V_H2O = 30 * 10**(-30)

# Simulation properties
rho_star = 2.24 # Not exactly 3 due to simulation conditions
dt_star = 0.04
r_c_star = 1
eta_star = 1
m_star = 1
temp_star = 1
temp = 298 # Predetermined mapping of 1 : 298K for temperature
D_star = 1

# CG factor of 3 H2O molecules to 1 bead
N_m = 3

# Cutoff r_c_star was defined as 1
epsilon = temp * k_b/temp_star
r_c = r_c_star * (rho_star * N_m * V_H2O)**(1/3)
m = m_star * (N_m * M_H2O)/N_A
rho = rho_star * (m/r_c**3)
dt = dt_star * (m**0.5 * r_c)/((temp * k_b)**0.5)
eta = eta_star * (m/(r_c * dt))
D = D_star * (r_c**2)/dt
P = epsilon/r_c**3

#%% 
data = np.genfromtxt("diff_water_data.txt")

fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 4])

diffusion = diffusion_star * D
diffusion_um = diffusion * 10**9

for i in range(len(subfolders)):
    ax.scatter(25, diffusion_um[i])

plt.title("Diffusion Coefficient against Temperature", fontsize = 14)
ax.set_xlabel("T ($\circ C$)", fontsize = 12)
ax.set_ylabel("$D\ (\mu m^2s^{-1})$", fontsize = 12)

fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders + ["Diffusion Curve"]), loc = "center right")

ax.set_xlim([0, 100])

