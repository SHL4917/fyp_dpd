print("Importing modules...")
import numpy as np
import sys
import os
from pathlib import Path
import shutil

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append(base)
from fyp_functions_short import *

#%%
bound = 36

filenames = os.listdir()
struct_filename = [x for x in filenames if "conc_" in x]
dense = [int(x.split("_")[1]) for x in struct_filename]
dense = [x/100 for x in dense]

polymer_masses = np.array([1, 1, 1])

num_bins = 400

dump_all = "dump.allparticles"
poly_count = np.zeros(np.size(dense))

#%%
n = 0
for i in dense:
    directory = f"{base}/conc_{i * 100:.0f}"
    os.chdir(directory)
    
    print(f"Extracting Coordinates from {directory}...")
    
    try:
        poly_count[n] = np.loadtxt("polymer_count.txt")
        coord1 = np.load("coord1.npy")
        coord2 = np.load("coord2.npy")
        coord3 = np.load("coord3.npy")
    except:
        print(f"Data files not found in {directory}!")
        os.chdir(base)
        
    
    for j in range(np.shape(coord1)[2]):
        c1 = coord1[:, :, j]
        c2 = coord2[:, :, j]
        c3 = coord3[:, :, j]
    
        t1, t2, t3, hist_axis = get_bilayer_histogram_norotate(c1[:, 0:3], 
                                                               c2[:, 0:3], 
                                                               c3[:, 0:3], 
                                                               1500, bound)
        if j == 0:
            hist1 = t1.copy()
            hist2 = t2.copy()
            hist3 = t3.copy()
            continue   
        hist1 = hist1 + t1
        hist2 = hist2 + t2
        hist3 = hist3 + t3
        
        print(f"Timestep: {j}")
        
    hist1 = hist1 / (j + 1)
    hist2 = hist2 / (j + 1)
    hist3 = hist3 / (j + 1)
    
    for j in range(np.shape(coord1)[2]):
        c1 = coord1[:, :, j]
        c2 = coord2[:, :, j]
        c3 = coord3[:, :, j]
        
        t_profile, p_axis = get_p_profile_norotate(np.concatenate([c1, c2, c3]), 
                                             num_bins, bound)
        
        if j == 0:
            p_profile = t_profile.copy()
            continue
        
        p_profile = p_profile + t_profile
        
        print(f"Timestep: {j}")
    
    p_profile = p_profile/(j + 1)   
    
    os.chdir(f"{base}/processed_data")
    
    np.save(f"hist1_conc_{i * 100:.0f}", hist1)
    np.save(f"hist2_conc_{i * 100:.0f}", hist2)
    np.save(f"hist3_conc_{i * 100:.0f}", hist3)
    np.save(f"hist_axis_conc_{i * 100:.0f}", hist_axis)
    
    np.save(f"p_profile_conc_{i * 100:.0f}", p_profile)
    np.save(f"p_profile_axis_conc_{i * 100:.0f}", p_axis)
    
    os.chdir(base)
    
    n += 1

os.chdir(f"{base}/processed_data")
np.savetxt("poly_count.txt", poly_count)
    
    