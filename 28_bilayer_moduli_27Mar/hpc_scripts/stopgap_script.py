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

sys.path.append("../../")
from fyp_functions import *

#%%
bound = 36
filenames = os.listdir()

struct = [x for x in filenames if "t4t4" in x]
for i in range(len(struct)):
    if len(struct[0]) > 6:
        struct = struct[1:] + [struct[0]]

if not os.path.exists("poly_num"):
    os.makedirs("poly_num")
                      
#%%
for i in range(len(struct)):
    os.chdir(f"{base}\{struct[i]}")
    filenames = os.listdir()
    
    conc = [x.split("_")[-1].split(".")[0] for x in filenames if "conc" in x]
    conc = [int(x) for x in conc]
    conc = np.array(conc)
    
    for j in conc:
        directory = f"{base}\{struct[i]}\conc_{j}"
        os.chdir(directory)
        
        poly_count = np.loadtxt("polymer_count.txt")
        
        os.chdir(f"{base}\poly_num")
        
        if not os.path.exists(f"{base}\poly_num\{struct[i]}"):
            os.makedirs(f"{base}\poly_num\{struct[i]}")
            
        os.chdir(f"{base}\poly_num\{struct[i]}")
        
        np.savetxt(f"poly_count_{j}.txt", [poly_count])
        
        
            
        