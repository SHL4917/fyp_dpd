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
    
#%% Extract data
snap_7 = np.load("temp_snap_7.npy")

#%% Get COM of all polymers
head = np.array([2, 2, 2, 2, 2])
tail1 = np.array([3, 3, 3, 3])
tail2 = np.array([3, 3, 3, 3])
polymer_struct = np.concatenate((head, tail1, tail2))
polymer_masses = np.array([1, 1, 1])
bounds = 30

#%% Get COM of all polymers - 2
com_poly = find_com(snap_7, polymer_struct, polymer_masses)
    
#%%
planes = get_planes(com_poly, 200, bounds, 5)

#%% Visualise
plane1 = generate_plane_points(bounds, planes[0, :], 30)
to_plot = np.concatenate((com_poly, plane1), axis = 0)
poly_colors = np.zeros(np.shape(com_poly))
plane_colors = np.zeros(np.shape(plane1))
plane_colors[:] = plane_colors[:] + np.array([1, 0, 0])

colors = np.concatenate((poly_colors, plane_colors), axis = 0)

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(to_plot)
pcd.colors = o3d.utility.Vector3dVector(colors)

o3d.visualization.draw_geometries([pcd, get_bounding_box(bounds)])

#%%
np.savetxt("plane_coeff.txt", planes[0, :])
    
    
    