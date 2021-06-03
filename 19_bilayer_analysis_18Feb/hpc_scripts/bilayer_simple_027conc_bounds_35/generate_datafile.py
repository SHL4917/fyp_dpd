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
bounds = 35
conc = 0.27
head = np.array([2, 2, 2, 2, 2])
tail1 = np.array([3, 3, 3, 3])
tail2 = np.array([3, 3, 3, 3])

#%%
polymer_coords = create_bilayer_datafile(bounds, conc, head, tail1, tail2)

#%%
np.save("polymer_coords", polymer_coords)
