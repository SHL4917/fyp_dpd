#%% Imports and setup
import numpy as np
import os
from pathlib import Path
from input_func import *

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%%
poly_num = create_poly_datafile(40, 0.05, np.array([2, 2, 2, 2, 3, 3, 3, 3]), 
                     filename = "e4o4.txt")

#%%
np.savetxt("e4o4_num.txt", poly_num, fmt = "%.d")