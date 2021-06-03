import numpy as np
import sys
import os
from pathlib import Path
import re
import pickle5 as pickle
import matplotlib.pyplot as plt
from scipy import stats

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%%
def outlier_elim(np_array):
    mad = stats.median_abs_deviation(np_array)
    median = np.median(np_array)
    filtered = np_array[(0.6745 * abs(np_array - median))/mad < 3.5]
    return filtered

#%%
folder_path = Path.cwd().parents[0]
foldernames = [x for x in  os.listdir(folder_path) if "conc" in x]

#%%
fig, ax = plt.subplots()
i = 0
valid = [1, 2, 3, 6, 8]

for filename in foldernames:
    if i not in valid:
        i += 1
        continue
    
    clusters = np.load(folder_path / filename / f"num_clusters_{filename}.npy")
    timestep = np.load(folder_path / filename / f"timestep_{filename}.npy")
    
    ax.scatter(timestep/1000, clusters, s = 8, label = f"{filename}")
    
    i += 1

plt.figtext(.5,.95, f"Number of Clusters against Timesteps", 
                fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1, 0.9))
ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
ax.set_ylabel(r"Count", fontsize = 12)

plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, timestep.max()/1000])

fig.savefig("clusters_against_timestep_selected.png", format='png', dpi=1200, 
                  bbox_inches='tight')  

#%%
equilibrium_timestep = np.array([100, 200, 300, 520, 500])
valid_foldernames = [foldernames[i] for i in valid]

i = 0
for filename in valid_foldernames:
    #if i != 0:
    #    i += 1
    #    continue
    
    struct = filename.split("_")[0]
    conc = float(filename.split("_")[-1])/100
    
    with open(folder_path/filename/f"rog_dict_{filename}.pickle", 
              "rb") as handle:
        rog = pickle.load(handle)
    
    cluster_rog = []
    for key, item in rog.items():
        if key < equilibrium_timestep[i]:
            continue
        
        for subkey, subitem in item.items():
            cluster_rog += [subitem]
    
    cluster_rog = np.array(cluster_rog)
    cluster_rog_filtered = outlier_elim(cluster_rog)
    
    fig, ax = plt.subplots(1, 2)
    ax[0].hist(cluster_rog, bins = 100)
    ax[1].hist(cluster_rog_filtered, bins = 100)
    
    plt.figtext(.5,1.17, f"Histogram of Radius of Gyration of Micelles", 
            fontsize=14, ha='center')
    
    plt.figtext(.5, 1.11, f"Micelle Structure: {struct}, $c={conc}$", 
            fontsize=12, ha='center')
    
    plt.figtext(.5, 1.05, f"Mean: {cluster_rog.mean():.3f}, "
                f"{cluster_rog_filtered.mean():.3f}", 
            fontsize=10, ha='center')
    
    plt.figtext(.5, 1, f"Std: {np.std(cluster_rog):.3f}, "
                f"{np.std(cluster_rog_filtered):.3f}", 
            fontsize=10, ha='center')
    
    plt.figtext(.3,.9, f"Raw Values", 
            fontsize=10, ha='center')
    
    plt.figtext(.73,.9, f"After Outlier Detection", 
            fontsize=10, ha='center')
    
    plt.figtext(.5,.0, "Non-Dimensional Radius of Gyration $R^*$", 
            fontsize=12, ha='center')
    ax[0].set_ylabel(r"Count", fontsize = 12)
    
    
    fig.savefig(f'rog_histogram_{filename}.png', format='png', dpi=1200, 
              bbox_inches='tight')
    
    i += 1

#%%

