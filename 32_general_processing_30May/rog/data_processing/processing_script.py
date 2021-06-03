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

sys.path.append("../../../")
from fyp_functions import *

#%%
def outlier_elim(np_array):
    mad = stats.median_abs_deviation(np_array)
    median = np.median(np_array)
    filtered = np_array[(0.6745 * abs(np_array - median))/mad < 3.5]
    return filtered

#%%
folder_path = Path.cwd().parents[0]
foldernames = [x for x in  os.listdir(folder_path) if "rog" in x]

struct = [x.split("_")[0] for x in foldernames]
struct = [int(re.split('(\d+)', x)[3]) for x in struct]

for i in range(len(struct)):
    if struct[0] > struct[-1]:
        struct = struct[1:] + [struct[0]]
        foldernames = foldernames[1:] + [foldernames[0]]

struct_latex = ["$H_{" + str(2 * x) + "} T_{" + str(x) + "}$" for x in struct]

#%%
fig, ax = plt.subplots()
i = 0
valid = [0, 1, 2, 3, 4]

for filename in foldernames:
    if i not in valid:
        i += 1
        continue
    
    clusters = np.load(folder_path / filename / f"num_clusters_{filename}.npy")
    timestep = np.load(folder_path / filename / f"timestep_{filename}.npy")
    
    ax.scatter(timestep/1000, clusters, s = 5, label = struct_latex[i])
    
    i += 1

plt.figtext(.5,.95, f"Number of Clusters against Timesteps", 
                fontsize=14, ha='center')

plt.legend(bbox_to_anchor = (1, 0.402995))
ax.set_xlabel("Timesteps (Thousands)", fontsize = 12)
ax.set_ylabel(r"Count", fontsize = 12)

plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([0, 400])
ax.set_ylim([0, 100])

# fig.savefig("clusters_against_timestep_justone.png", format='png', dpi=1200, bbox_inches='tight')  

#%%
equilibrium_timestep = np.array([150, 150, 150, 150, 150])
valid_foldernames = [foldernames[i] for i in valid]

#%%
i = 0
for filename in valid_foldernames:
    # if i != 0:
    #     i += 1
    #     continue
    
    conc = 0.07
    
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
    
    fig, ax = plt.subplots()
    ax.hist(cluster_rog, bins = 100)
    
    plt.figtext(.5,1.1, f"Histogram of Radius of Gyration of Micelles", 
            fontsize=14, ha='center')
    
    plt.figtext(.5, 1.03, f"Micelle Structure: " + struct_latex[i] + 
                f", $c={conc}$", 
            fontsize=12, ha='center')
    
    plt.figtext(.5, 0.98, f"Mean: {cluster_rog.mean():.3f}",
            fontsize=10, ha='center')
    
    plt.figtext(.5, 0.93, f"Std: {np.std(cluster_rog):.3f}", 
            fontsize=10, ha='center')
    
    plt.figtext(.5,.0, "Non-Dimensional Radius of Gyration $R^*_{micelle}$", 
            fontsize=12, ha='center')
    ax.set_ylabel(r"Count", fontsize = 12)
    
    
    fig.savefig(f'rog_histogram_{filename}.png', format='png', dpi=1200,bbox_inches='tight')
    
    i += 1

#%%
fig, ax = plt.subplots()
estimates = np.array([2.3, 3, 3.5, 2.5, 2.6])
ax.scatter(struct, estimates, s = 20)
ax.plot(struct, estimates, ls = "--", lw = "1")

plt.figtext(.5,.93, f"Estimated Mode of Radius of Gyration against Hydrophilic Bead Number", 
            fontsize=14, ha='center')

ax.set_xlabel("Number of Hydrophilic Beads", fontsize = 12)
ax.set_ylabel(r"$R^*_{micelle}$", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)

ax.set_xlim([1, 7])
ax.set_ylim([1, 4])    

#fig.savefig(f'rog_trend.png', format='png', dpi=1200, bbox_inches='tight')



