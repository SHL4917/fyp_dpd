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

#%% Plotting pressure and temperature
fig, ax = plt.subplots(1, 2)
labels = len(a_ij_values)

for key, entry in a_ij_values.items():
    ax[0].scatter(timestep_av_values[key], avpressure[key], marker = ".")
    ax[1].scatter(timestep_av_values[key], avtemp[key], marker = ".")


ax[0].set_ylim([0, 50])
ax[1].set_ylim([0, 1.5])

ax[0].set_title("Reduced Pressure", fontsize = 14)
ax[1].set_title("Reduced Temperature", fontsize = 14)

ax[0].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[0].set_ylabel("P*", fontsize = 12)

ax[1].set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax[1].set_ylabel("T*", fontsize = 12)

fig.legend(handles = (ax[0].get_children() ), 
           labels = (a_ij_values), loc = "center right")

fig.tight_layout()
plt.show()

#%% Plotting velocity profile -vx
fig, ax = plt.subplots()
for key, entry in a_ij_values.items():
    x = vprofile_vx[key][45, :, 1]
    y = vprofile_vx[key][45, :, 3]
    ax.scatter(x, y)
    ax.plot(x, y)
    # ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

plt.title("Velocity Profile in the z-direction", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("Vx*", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (a_ij_values), loc = "center right")

plt.show()

#%% Plotting number density profile
fig, ax = plt.subplots()
for key, entry in a_ij_values.items():
    x = rhoprofile[key][45, :, 1]
    y = rhoprofile[key][45, :, 3]
    ax.scatter(x, y)
    
plt.title("Number density profile", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel(r"$\rho^*$", fontsize = 12)
ax.set_ylim([0, 4])
fig.legend(handles = (ax.get_children() ), 
           labels = (a_ij_values), loc = "center right")

plt.show()

#%% Plotting temperature profile
fig, ax = plt.subplots()
for key, entry in a_ij_values.items():
    x = tprofile[key][45, :, 1]
    y = tprofile[key][45, :, 3]
    ax.scatter(x, y)
    
plt.title("Temperature profile", fontsize = 14)
ax.set_xlabel("z (r*)", fontsize = 12)
ax.set_ylabel("$T^*$", fontsize = 12)
ax.set_ylim([0, 2])
fig.legend(handles = (ax.get_children() ), 
           labels = (a_ij_values), loc = "lower right")

plt.show()


#%% Get viscosities, assuming couette flow in x-direction
shearrate_star = {}
for key, entry in a_ij_values.items():
    shearrate_star[key] = np.zeros(len(vprofile_vx[key]))
    
    for i in range(len(shearrate_star[key])):
        x = vprofile_vx[key][i, 3:7, 1]
        y = vprofile_vx[key][i, 3:7, 3]
        shearrate_star[key][i] = np.polyfit(x, y, 1)[0]

#%% Plot shear rate w.r.t. timestep
fig, ax = plt.subplots()
for key, entry in a_ij_values.items():
    ax.scatter(timestep_velo_profile[key], shearrate_star[key], 
               marker = ".")

ax.set_ylim([0, 0.5])
plt.title("Shear Rate w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$\gamma^*$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (a_ij_values), loc = "center right")

plt.show()

#%% Plot pressure in xz direction
fig, ax = plt.subplots()
for key, entry in a_ij_values.items():
    ax.scatter(timestep_press[key], press_xz[key], marker = ".")
    
plt.title("Pressure (XZ component) w.r.t. Timestep", fontsize = 14)
ax.set_xlabel("Steps ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("$P_{xz}^*$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (a_ij_values), loc = "center right")
ax.set_ylim([-0.2, 0])

plt.show()

#%% Get viscosities from the moving average pressure and shearrate values
viscosity_star = {}
for key, entry in a_ij_values.items():
    x_ave, y_ave = movingaverage(timestep_press[key], press_xz[key])
    x1_ave, y1_ave = movingaverage(timestep_velo_profile[key], 
                                   shearrate_star[key])
    viscosity_star[key] = -y_ave[-1]/y1_ave[-1]

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
    dt[key] = dt_star * (m**0.5 * r_c[key])/((temp * k_b)**0.5)
    eta[key] = eta_star * (m/(r_c[key] * dt[key]))
    P[key] = epsilon/r_c[key]**3

#%% Get actual viscosities and plot
data = np.genfromtxt("temp_viscosity_data.txt")

fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 2])


t_simulation  = {}

viscosity = {}
for key, entry in a_ij_values.items():
    viscosity[key] = eta[key] * viscosity_star[key]
    x_ave, y_ave = movingaverage(timestep_av_values[key], avtemp[key])
    t_simulation[key] = y_ave[-1] * (epsilon/k_b) - 273.15
    ax.scatter(t_simulation[key], viscosity[key]) 


plt.title("Dynamic Viscosity against Temperature", fontsize = 14)
ax.set_xlabel("T ($\circ C$)", fontsize = 12)
ax.set_ylabel("$\eta\ (Pa\ s)$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (subfolders + ["Viscosity Curve"]), loc = "center right")

plt.show()

#%% Plot certain viscosities
fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 2])


t_simulation  = {}


x_ave, y_ave = movingaverage(timestep_av_values[1], avtemp[1])
t_simulation[1] = y_ave[-1] * (epsilon/k_b) - 273.15
ax.scatter(t_simulation[1], viscosity[1])

x_ave, y_ave = movingaverage(timestep_av_values[1.01], avtemp[1.01])
t_simulation[1.01] = y_ave[-1] * (epsilon/k_b) - 273.15
ax.scatter(t_simulation[1.01], viscosity[1.01])

x_ave, y_ave = movingaverage(timestep_av_values[1.1], avtemp[1.1])
t_simulation[1.1] = y_ave[-1] * (epsilon/k_b) - 273.15
ax.scatter(t_simulation[1.1], viscosity[1.1])

x_ave, y_ave = movingaverage(timestep_av_values[1.2], avtemp[1.2])
t_simulation[1.2] = y_ave[-1] * (epsilon/k_b) - 273.15
ax.scatter(t_simulation[1.2], viscosity[1.2])

x_ave, y_ave = movingaverage(timestep_av_values[1.3], avtemp[1.3])
t_simulation[1.3] = y_ave[-1] * (epsilon/k_b) - 273.15
ax.scatter(t_simulation[1.3], viscosity[1.3])

plt.title("Dynamic Viscosity against Temperature", fontsize = 14)
ax.set_xlabel("T ($\circ C$)", fontsize = 12)
ax.set_ylabel("$\eta\ (Pa\ s)$", fontsize = 12)
fig.legend(handles = (ax.get_children() ), 
           labels = (["rc_1", "rc_1.01", "rc_1.1", "rc_1.2", "rc_1.3"] + 
                     ["Viscosity Curve"]), loc = "center right")

plt.show()

