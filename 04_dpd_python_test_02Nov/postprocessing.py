#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)

#%% Load datafiles - rdf histogram
data = np.loadtxt("rdf.txt", skiprows = 4)
bin_centers = data[:, 1]
rdf = data[:, 2]

#%% Set bin edges 
bins = bin_centers - bin_centers[0]

#%% Plot histogram
fig, ax = plt.subplots()
ax.scatter(bin_centers, rdf, marker = ".")

plt.title("Histogram of the Radial Distribution Function", fontsize = 14)
ax.set_xlabel("$r^*$", fontsize = 12)
ax.set_ylabel("$g(r^*)$", fontsize = 12)
plt.show()

#%% Plot histogram - line graph
fig, ax = plt.subplots()
ax.plot(bin_centers, rdf)

plt.title("Histogram of the Radial Distribution Function", fontsize = 14)
ax.set_xlabel("$r^*$", fontsize = 12)
ax.set_ylabel("$g(r^*)$", fontsize = 12)
plt.show()
