#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import animation
import pandas as pd

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

#%% Load data
poly_data = np.load("polymer_coords.npy")

#%% Function definitions
def get_e2e(first_slice):
    e2e = np.zeros(np.shape(first_slice)[0])
    
    for i in range(np.shape(first_slice)[0]):
        vec = first_slice[i, 0, :] - first_slice[i, -1, :]
        e2e[i] = np.linalg.norm(vec)
    
    return e2e

def find_com(first_slice, polymer_struct, polymer_masses):
    totmass = 0
    for i in range(len(polymer_struct)):
        totmass += polymer_masses[polymer_struct[i] - 1]
    
    com = np.zeros((np.shape(first_slice)[0], 3))

    for i in range(np.shape(first_slice)[0]):
        denom = np.array([0, 0, 0])
        for j in range(np.shape(first_slice)[1]):
            denom = denom + polymer_masses[polymer_struct[j] - 1] * \
                first_slice[i, j, :]
        
        com[i, :] = denom/totmass
    
    return com

def find_r_g(first_slice, com_slice):
    r_g = np.zeros(np.shape(first_slice)[0])

    for i in range(np.shape(first_slice)[0]):
        vec = 0
        for j in range(np.shape(first_slice)[1]):
            vec = vec + np.linalg.norm(first_slice[i, j, :] - \
                                       com_slice[i, :]) ** 2
        
        r_g[i] = (vec/np.shape(first_slice)[1]) ** 0.5
    
    return r_g

#%% Set parameters - Ideally should have been saved along with the poly_coords
x_bounds = np.array([0, 40])
y_bounds = np.array([0, 40])
z_bounds = np.array([0, 40])

x_length = x_bounds[1] - x_bounds[0]
y_length = y_bounds[1] - y_bounds[0]
z_length = z_bounds[1] - z_bounds[0]

boxlengths = np.array([x_length, y_length, z_length])

#%% Get the end-to-end vector for all timesteps
e2e = np.zeros((np.shape(poly_data)[0], (np.shape(poly_data)[3])))
for i in range(np.shape(poly_data)[3]):
    e2e[:, i] = get_e2e(poly_data[:, :, :, i])

#%% Get average of end-to-end vector for each timestep
e2e_ave = np.zeros(np.shape(e2e)[1])
e2e_std = np.zeros(np.shape(e2e)[1])
dt = np.zeros(np.shape(e2e)[1])
for i in range(np.size(e2e_ave)):
    e2e_ave[i] = np.average(e2e[:, i])
    e2e_std[i] = np.std(e2e[:, i])
    dt[i] = i * 500

#%% Calculate COM of each polymer
polymer_struct = np.array([2, 2, 2, 2, 3, 3, 3, 3])
polymer_masses = np.array([1, 1.12, 0.866])

com = np.zeros((np.shape(poly_data)[0], 3, np.shape(poly_data)[3]))
for i in range(np.shape(poly_data)[3]):
    com[:, :, i] = find_com(poly_data[:, :, :, i], polymer_struct,
                            polymer_masses)

#%% Calculate radius of gyration
r_g = np.zeros((np.shape(poly_data)[0], np.shape(poly_data)[3]))
for i in range(np.shape(poly_data)[3]):
    r_g[:, i] = find_r_g(poly_data[:, :, :, i], com[:, :, i])
    
#%% Get average and std of radius of gyration for each timestep
r_g_ave = np.zeros(np.shape(r_g)[1])
r_g_std = np.zeros(np.shape(r_g)[1])
dt = np.zeros(np.shape(r_g)[1])
for i in range(np.size(r_g_ave)):
    r_g_ave[i] = np.average(r_g[:, i])
    r_g_std[i] = np.std(r_g[:, i])
    dt[i] = i * 500
    
#%% Plot average e2e with std
lo = 150
hi = 200
fig, ax = plt.subplots()
ax.scatter(dt[lo:hi], e2e_ave[lo:hi], marker = ".", label = "Mean")
ax.errorbar(dt[lo:hi], e2e_ave[lo:hi], yerr = e2e_std[lo:hi], ls = "none", 
            capsize = 3, label = "One Standard Deviation") 

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

plt.ylim([0, 5])

plt.title("Plot of Average End-to-End Polymer Distances")
ax.set_xlabel("Timestep (dt = 0.04)")
ax.set_ylabel("End-to-End Distance $R$")

plt.show()

#%% Plot the histograms?
fig, ax = plt.subplots()
ax.hist(e2e[:, 250], bins = 50)

plt.title("Histogram of End-to-End Polymer Distances")
ax.set_xlabel("Distance $R$")
ax.set_ylabel("Instances")
plt.show()

#%% Plot average gyration radius with std
lo = 150
hi = 200
fig, ax = plt.subplots()
ax.scatter(dt[lo:hi], r_g_ave[lo:hi], marker = ".", label = "Mean")
ax.errorbar(dt[lo:hi], r_g_ave[lo:hi], yerr = r_g_std[lo:hi], ls = "none", 
            capsize = 3, label = "One Standard Deviation") 

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

plt.ylim([0, 2])

plt.title("Plot of Average Radius of Gyrations for Polymers")
ax.set_xlabel("Timestep (dt = 0.04)")
ax.set_ylabel("Radius of Gyration $R_g$")

plt.show()

#%% Plot some histograms
fig, ax = plt.subplots()
ax.hist(r_g[:, 250], bins = 50)

plt.title("Histogram of Radius of Gyrations for Polymers")
ax.set_xlabel("Distance $R$")
ax.set_ylabel("Instances")
plt.show()

#%% Save relevant files
np.save("com.npy", com)