#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import math

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% Functions
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

#%% Get foldernames
base = os.getcwd()
subfolders = [name for name in os.listdir(".") if 
              os.path.isdir(os.path.join(base, name))]

a_ij_values = [float(x.split("_")[-1]) for x in subfolders]

#%% Load datafiles - Average quantities
avpressure = {}
avtemp = {}
press_xz = {}
timestep_av_values = {}
timestep_press = {}

for i in range(len(subfolders)):
    data = np.loadtxt(subfolders[i] + "/" + "avpressure.txt")
    timestep_av_values[a_ij_values[i]] = data[:, 0]
    avpressure[a_ij_values[i]] = data[:, 1]
    
    data = np.loadtxt(subfolders[i] + "/" + "avtemp.txt")
    avtemp[a_ij_values[i]] = data[:, 1]
    
    data = np.loadtxt(subfolders[i] + "/" + "press_xz.txt")
    timestep_press[a_ij_values[i]] = data[:, 0]
    press_xz[a_ij_values[i]] = data[:, 1]

#%% Load datafiles - Velocity profiles
timestep_velo_profile = {}
vprofile_vx = {}

box_length = 10

for i in range(len(subfolders)):
    file_name = subfolders[i] + "/" + "velo_profile_vx.txt"
    data = np.genfromtxt(file_name, invalid_raise = False)
    timestep_velo_profile[a_ij_values[i]] = data[:, 0]
    v_vx = np.genfromtxt(file_name, invalid_raise = False, skip_header = 4)
    
    num_bins = data[0, 1].astype(int)
    num_entries = (len(v_vx)/num_bins).astype(int)
    
    v_vx = np.reshape(v_vx, (num_entries, num_bins, 4))
    v_vx[:, :, 1] = v_vx[:, :, 1] * box_length
    vprofile_vx[a_ij_values[i]] = v_vx    

#%% Load datafiles - z-direction velocities, same horizonal bins as x-direction
vprofile_vz = {}

box_length = 10

for i in range(len(subfolders)):
    file_name = subfolders[i] + "/" + "velo_profile_vz.txt"
    data = np.genfromtxt(file_name, invalid_raise = False)
    v_vz = np.genfromtxt(file_name, invalid_raise = False, skip_header = 4)
    
    num_bins = data[0, 1].astype(int)
    num_entries = (len(v_vz)/num_bins).astype(int)
    
    v_vz = np.reshape(v_vz, (num_entries, num_bins, 4))
    v_vz[:, :, 1] = v_vz[:, :, 1] * box_length
    vprofile_vz[a_ij_values[i]] = v_vz   

#%% Plotting pressure and temperature
fig, ax = plt.subplots(1, 2)

for i in range(len(subfolders)):
    entry = a_ij_values[i]
    ax[0].scatter(timestep_av_values[entry], avpressure[entry], marker = ".")
    ax[1].scatter(timestep_av_values[entry], avtemp[entry], marker = ".")


ax[0].set_ylim([0, 1.5])
ax[1].set_ylim([0, 1.5])

ax[0].set_title("Reduced Pressure", fontsize = 14)
ax[1].set_title("Reduced Temperature", fontsize = 14)

ax[0].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[0].set_ylabel("P*", fontsize = 12)

ax[1].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[1].set_ylabel("T*", fontsize = 12)

fig.legend(handles = (ax[0].get_children() ), 
           labels = (subfolders), loc = "center right")

fig.tight_layout()
plt.show()

#%% Plotting velocity profile -vx
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    x = vprofile_vx[entry][80, :, 1]
    y = vprofile_vx[entry][80, :, 3]
    ax.scatter(x, y)
    ax.plot(x, y)
    # ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("Vx*", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")

plt.show()

#%% Plotting velocity profile - vz
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    x = vprofile_vz[entry][88, :, 1]
    y = vprofile_vz[entry][88, :, 3]
    ax.scatter(x, y)

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("Vz*", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")

plt.show()

#%% Get viscosities, assuming couette flow in x-direction
shearrate_star = {}
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    shearrate_star[entry] = np.zeros(len(vprofile_vx[entry]))
    
    for i in range(len(shearrate_star[entry])):
        x = vprofile_vx[entry][i, :, 1]
        y = vprofile_vx[entry][i, :, 3]
        shearrate_star[entry][i] = np.polyfit(x, y, 1)[0]

#%% Plot shear rate w.r.t. timestep
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    ax.scatter(timestep_velo_profile[entry], shearrate_star[entry], 
               marker = ".")

ax.set_ylim([0, 0.5])
plt.title("Shear Rate w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$\gamma^*$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")

plt.show()

#%% Plot pressure in xz direction
fig, ax = plt.subplots()
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    ax.scatter(timestep_press[entry], press_xz[entry], marker = ".")
    
plt.title("Pressure (XZ component) w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$P_{xz}^*$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders), loc = "center right")
ax.set_ylim([-0.08, 0])

plt.show()

#%% Get viscosities from the moving average pressure and shearrate values
viscosity_star = np.zeros(len(subfolders))
for i in range(len(subfolders)):
    entry = a_ij_values[i]
    x_ave, y_ave = movingaverage(timestep_press[entry], press_xz[entry])
    x1_ave, y1_ave = movingaverage(timestep_velo_profile[entry], 
                                   shearrate_star[entry])
    viscosity_star[i] = -y_ave[-1]/y1_ave[-1]

#%% Performing conversion to real units
N_A = 6.02214 * 10**(23)
k_b = 1.38064852 * 10**(-23)

# Volume Properties of water
M_H2O = 18.01528 * 10**(-3)
V_H2O = 30 * 10**(-30)

# Simulation properties
rho_star = 3375/1000 # Not exactly 3 due to simulation conditions
dt_star = 0.04
r_c_star = 1
eta_star = 1
m_star = 1
temp_star = 1
temp = 298 # Predetermined mapping of 1 : 298K for temperature

# CG factor of 3 H2O molecules to 1 bead
N_m = 3

# Cutoff r_c_star was defined as 1
epsilon = temp * k_b/temp_star
r_c = r_c_star * (rho_star * N_m * V_H2O)**(1/3)
m = m_star * (N_m * M_H2O)/N_A
rho = rho_star * (m/r_c**3) # Results in grams!
dt = dt_star * (m**0.5 * r_c)/((temp * k_b)**0.5)
eta = eta_star * (m/(r_c * dt))

#%% Get actual viscosities and plot
data = np.genfromtxt("temp_viscosity_data.txt")

fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 2])

viscosity = eta * viscosity_star
t_simulation = np.zeros(len(subfolders))

for i in range(len(subfolders)):
    entry = a_ij_values[i]
    x_ave, y_ave = movingaverage(timestep_av_values[entry], avtemp[entry])
    t_simulation[i] = y_ave[-1] * (epsilon/k_b) - 273.15
    ax.scatter(t_simulation[i], viscosity[i]) 


plt.title("Dynamic Viscosity against Temperature", fontsize = 14)
ax.set_xlabel("T ($\circ C$)", fontsize = 12)
ax.set_ylabel("$\eta\ (Pa\ s)$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders + ["Viscosity Curve"]), loc = "center right")

plt.show()


