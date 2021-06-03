#%% Imports and setup
import numpy as np
import os
from pathlib import Path
from input_func import create_poly_datafile

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%%
for index, val in enumerate([2, 3, 4]):
    for index2, conc in enumerate([0.01, 0.03, 0.05]):
        filename = str(val) + "chains_" + str(conc)
        if not os.path.exists(filename):
            os.makedirs(filename)
        
        os.chdir(filename)
        polystruct = np.ones(val) * 3
        polystruct[0] = 2
        
        create_poly_datafile(18, conc, polystruct)
        
        os.chdir(__dir__)

    