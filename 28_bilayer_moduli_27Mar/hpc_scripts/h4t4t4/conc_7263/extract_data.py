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

bounds = 36

dump_all = "dump.allparticles"

coord1, coord2, coord3 = get_allparticles_from_dump2(dump_all)

np.save("coord1.npy", coord1)
np.save("coord2.npy", coord2)
np.save("coord3.npy", coord3)