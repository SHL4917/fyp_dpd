import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt


__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../")
from fyp_functions import *

#%%
poly_coords = get_poly_from_dump("dump.atom", "e4o4_num.txt")

#%%
np.save("poly_coords", poly_coords)

#%%
com = get_com_all(poly_coords, np.array([2, 2, 2, 2, 3, 3, 3, 3]))

#%%
np.save("com", com)