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

def get_com_all(poly_data, polymer_struct, 
                polymer_masses = np.array([1, 1.12, 0.866])):
    
    com = np.zeros((np.shape(poly_data)[0], 3, np.shape(poly_data)[3]))
    for i in range(np.shape(poly_data)[3]):
        com[:, :, i] = find_com(poly_data[:, :, :, i], polymer_struct,
                                polymer_masses)
        
    return com

def get_r_g_all(poly_data, com):
    r_g = np.zeros((np.shape(poly_data)[0], np.shape(poly_data)[3]))
    for i in range(np.shape(poly_data)[3]):
        r_g[:, i] = find_r_g(poly_data[:, :, :, i], com[:, :, i])
        
    return r_g

def get_r_g_stats(r_g):    
    r_g_ave = np.zeros(np.shape(r_g)[1])
    r_g_std = np.zeros(np.shape(r_g)[1])
    dt = np.zeros(np.shape(r_g)[1])
    for i in range(np.size(r_g_ave)):
        r_g_ave[i] = np.average(r_g[:, i])
        r_g_std[i] = np.std(r_g[:, i])
        dt[i] = i * 500
    
    return r_g_ave, r_g_std

#%% Iterating over the different folders
r_g_array = np.zeros(7)
mass_array = np.zeros(7)

for i in range(3, 10, 1):
    print(i)
    folder = "e" + str(i) + "/"
    polymer_struct = np.zeros(i, dtype = int) + 2
    
    poly_data = np.load(folder + "polymer_coords.npy")
    com = get_com_all(poly_data, polymer_struct)
    np.save(folder + "com.npy", com)
    
    r_g_ave, r_g_std = get_r_g_stats(get_r_g_all(poly_data, com))
    r_g_array[i - 3] = np.average(r_g_ave)
    mass_array[i - 3] = i * 396

#%% Plot 
fig, ax = plt.subplots()
ax.scatter(np.log(mass_array), np.log(r_g_array), marker = ".", 
           label = "Points")

c = np.polyfit(np.log(mass_array[1:]), np.log(r_g_array[1:]), deg = 1)
x_vals = np.array(ax.get_xlim())
y_vals = c[0] * x_vals + c[1]

ax.plot(x_vals, y_vals, ls = "--", label = f"Gradient = {c[0]:.3f}")

plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
plt.title("Radius of Gyration against Polymer Molar Mass")
ax.set_ylabel("$ln(R_g)$")
ax.set_xlabel("$ln(M)$")
plt.legend(loc = "lower right")

# fig.savefig('betterfit.png', format='png', dpi=1200, bbox_inches='tight')




