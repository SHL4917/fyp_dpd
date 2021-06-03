import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import shutil
import re
import scipy as sp

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../")
from fyp_functions import *

#%%
bound = 36
filenames = os.listdir()

struct = [x for x in filenames if "t4t4" in x]
for i in range(len(struct)):
    if len(struct[0]) > 6:
        struct = struct[1:] + [struct[0]]

struct_latex = ["$H_{" + re.split('(\d+)', x)[1] + \
                "} T_4 T_4$" for x in struct]
number_p_mol = [int(re.split('(\d+)', x)[1]) + 8 for x in struct]


#%%
for i in range(len(struct)):
    
    if i != 7:
        continue
    
    os.chdir(f"{base}\{struct[i]}\processed_data")
    filenames = os.listdir()
    
    conc = [x.split("_")[-1].split(".")[0] for x in filenames if "hist_axis" in x]
    conc = [int(x) for x in conc]
    conc = np.array(conc)

    
    # for j in range(len(conc)):
    #     p_profile = np.load(f"p_profile_conc_{conc[j]}.npy") 
    #     p_axis = np.load(f"p_profile_axis_conc_{conc[j]}.npy")
        
    #     poly_count = np.loadtxt(f"poly_count_{conc[j]}.txt")
        
    #     fig, ax = plt.subplots()
    #     ax.scatter(p_axis, p_profile, s = 15)
    #     ax.set_xlim([0, bound])
        
    #     plt.figtext(.5,.97, 
    #                 f"Pressure Profile of Bilayer - {struct_latex[i]}", 
    #                 fontsize=14, ha='center')
        
    #     plt.figtext(.5,.91, f"Concentration of "
    #                 f"{poly_count * number_p_mol[i]/(3 * bound**3):.3f}", 
    #                 fontsize=11, ha='center')
        
    #     fig.savefig(f"pressure_profile_{struct[i]}_"
    #                 f"{poly_count * number_p_mol[i]/(3 * bound**3):.3f}.png", 
    #                 format='png', dpi=1200, bbox_inches='tight')
    #     plt.clf()
    #     plt.close(fig)   
    
    tension = np.zeros(len(conc))
    preff_area = np.zeros(len(conc))
    
    fig, ax = plt.subplots()
    
    for j in range(len(conc)):   
    
        p_profile = np.load(f"p_profile_conc_{conc[j]}.npy") 
        p_axis = np.load(f"p_profile_axis_conc_{conc[j]}.npy")
        
        poly_count = np.loadtxt(f"poly_count_{conc[j]}.txt")
        
        lower = len(p_axis[p_axis < 0]) - 1
        upper = len(p_axis[p_axis < bound])
        
        tension[j] = sp.integrate.simps(p_profile[lower:upper], 
                                        p_axis[lower:upper])
        preff_area[j] = (bound**2) / (poly_count * 0.5)
        
    preff_area, tension = zip(*sorted(zip(preff_area, tension)))
    preff_area = np.array(preff_area)
    tension = np.array(tension)
        
    ax.scatter(preff_area, -tension, s = 20)
    ax.plot(preff_area, -tension, linestyle = "--")
        
    plt.figtext(.5,.97, "Plot of Surface Tension against Polymer Projected Area", 
            fontsize=14, ha='center')

    plt.figtext(.5,.91,
            f"Polymer Structure: {struct_latex[i]}", 
            fontsize=11, ha='center')
    
    ax.set_xlabel("Projected Area $A^*$", fontsize = 12)
    ax.set_ylabel(r"Surface Tension $\sigma^*_z$", fontsize = 12)
    plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
    
    fig.savefig(f'surface_tension_vs_area_{struct[i]}.png', format='png', dpi=1200, 
            bbox_inches='tight')
    
    np.savetxt(f"preff_area_vs_tension_{struct[i]}.txt", 
           np.concatenate(([preff_area], [tension]), axis = 0).T)
    
    os.chdir(base)