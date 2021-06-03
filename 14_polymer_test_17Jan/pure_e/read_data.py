#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

#%% Function definitions
def unwrap_polymers(first_slice, boxlengths):
    chain_length = np.shape(first_slice)[1]
    for i in range(np.shape(first_slice)[0]):
        for j in range(chain_length):
            for k in range(3):
                diff = first_slice[i, j, k] - first_slice[i, 0, k]
                if np.abs(diff) > boxlengths[k] * 0.7 and np.sign(diff) > 0:
                    first_slice[i, j, k] = first_slice[i, j, k] - boxlengths[k]
                elif np.abs(diff) > boxlengths[k] * 0.7 and np.sign(diff) < 0:
                    first_slice[i, j, k] = first_slice[i, j, k] + boxlengths[k]
    
    return first_slice

def get_poly_from_dump(dump, polyfile):
    polymer_num = np.loadtxt(polyfile).astype("int")
    polymer_coord = np.zeros(np.shape(polymer_num) + (3,1))
    
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1
    
    x_length = x_bounds[1] - x_bounds[0]
    y_length = y_bounds[1] - y_bounds[0]
    z_length = z_bounds[1] - z_bounds[0]
    
    boxlengths = np.array([x_length, y_length, z_length])    
    
    # Load 
    i = 0 # Just for debugging line info
    timestep = 0
    skip_count = 0
    skip = False
    snapshot_coord = np.zeros(np.shape(polymer_num) + (3,1))
    with open(dump) as infile:
        for line in infile:
            i = i + 1
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot_coord):
                    polymer_coord = np.concatenate((polymer_coord, snapshot_coord)
                                                   ,axis = 3)
                    snapshot_coord = np.zeros(np.shape(polymer_num) + (3,1))
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            if data[1] == 1:
                continue
            elif data[1] == 2 or data[1] == 3:
                row, col = np.where(polymer_num == data[0].astype("int"))
                snapshot_coord[row, col, 0] = data[2] * x_length
                snapshot_coord[row, col, 1] = data[3] * y_length
                snapshot_coord[row, col, 2] = data[4] * z_length
    
    polymer_coord = polymer_coord[:, :, :, 1:]
    
    for i in range(np.shape(polymer_coord)[3]):
        polymer_coord[:, :, :, i] = unwrap_polymers(polymer_coord[:, :, :, i], 
                                                    boxlengths)
        print(i)
        
    return polymer_coord

#%% Get polymer coordinates for each folder and save in each folder
for i in range(3, 10, 1):
    polymer_coord = get_poly_from_dump("e" + str(i) + "/dump.atom", 
                                       "e" + str(i) + "/e" + str(i) + ".txt")
    np.save("e" + str(i) + "/polymer_coords.npy", polymer_coord)



        
            