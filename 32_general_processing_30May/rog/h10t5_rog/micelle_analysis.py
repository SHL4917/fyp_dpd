import numpy as np
import sys
import os
from pathlib import Path
import re
import pickle

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append(base)
from fyp_functions_short import *

"""
To do:
    Get plot of num_clusters against time
    Get average cluster radius of gyration, averaged over all particles and 
        over the period of which equilibrium is reached
    Get COM movement over the equilibrium timesteps for individual polymers:
        Calculate estimate for individual diffusion coefficient
        Can generalize to micelle coefficient (although not the same), as the
            micelle can be assumed to have a low rotational component and thus
            movement of the COM of individual particles should mirror the bulk
            movement of the micelle
    
    Simulate micelles with different compositions and get a relationship
        between the composition and:
            equilibrium time, radius of gyration, diffusion coefficient

"""

#%%
dump = "dump.poly"


filename = base.split("/")[-1]
folder = filename.split("_")[0]
struct = re.split('(\d+)', folder)
struct = np.array([int(struct[1]), int(struct[3])])

polyfile = f"rog_{struct[1]:.0f}.npy"

polymer_struct = np.zeros(struct.sum(), dtype = int)
for i in range(len(polymer_struct)):
    polymer_struct[i] = 2 if i < struct[0] else 3

dt = 1000

#%%
poly_coord = get_poly_from_dump(dump, polyfile)
poly_com = get_com_all(poly_coord, polymer_struct)

np.save(f"poly_coord_{filename}.npy", poly_coord)
np.save(f"poly_com_{filename}.npy", poly_com)

#%%
poly_coord = np.load(f"poly_coord_{filename}.npy")
poly_com = np.load(f"poly_com_{filename}.npy")

#%%
num_clusters = np.zeros(np.shape(poly_com)[2])
timestep = np.zeros(np.shape(poly_com)[2])
rog_dict = {}

for i in range(np.shape(poly_com)[2]):
    print(f"Timestep: {i}")
    cluster_dict, poly_dict = get_clusters(poly_com[:, :, i], 
                                           poly_coord[:, :, :, i], 36, 2)
    
    gyration_dict = {}
    for key, val in poly_dict.items():
        rog = 0
        cluster_com = cluster_dict[key].sum(axis = 0)/len(cluster_dict[key])
        for j in range(np.shape(val)[0]):
            dist = ((val[j, :, :] - cluster_com)**2).sum(axis = 1).max()
            rog += dist
        
        rog = rog/np.shape(val)[0]
        rog = rog**0.5
        gyration_dict[key] = rog
                 
    num_clusters[i] = len(cluster_dict)
    timestep[i] = i * dt
    rog_dict[i] = gyration_dict

#%%
np.save(f"num_clusters_{filename}.npy", num_clusters)
np.save(f"timestep_{filename}.npy", timestep)

with open(f"rog_dict_{filename}.pickle", 'wb') as handle:
    pickle.dump(rog_dict, handle, protocol = pickle.HIGHEST_PROTOCOL)

# To read:
"""
with open("pickle name", "rb") as handle:
    variable = pickle.load(handle)
"""    

#%%




