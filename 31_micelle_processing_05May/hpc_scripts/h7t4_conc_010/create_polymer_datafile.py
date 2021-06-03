import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
bound = 36
conc = 0.10
struct = np.array([2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3])
mass_array = np.array([1, 1, 1])

#%%
poly_num = create_poly_datafile(bound, conc, struct, mass_array = mass_array)
np.save("poly_num", poly_num)