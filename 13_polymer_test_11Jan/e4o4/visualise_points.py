#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import animation
import pandas as pd
import open3d as o3d

from sklearn import cluster
from sklearn.neighbors import NearestNeighbors

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

#%% Functions
def unwrap_pbc(coords, bounds):    
    to_add = coords.copy()   
    for i in range(3):
        for j in range(3):
            for k in range(3):
                copied = to_add.copy()
                copied[:, 0] = copied[:, 0] + i * bounds
                copied[:, 1] = copied[:, 1] + j * bounds
                copied[:, 2] = copied[:, 2] + k * bounds
                coords = np.concatenate((coords, copied), axis = 0)
    
    coords = coords[np.shape(to_add)[0]:, :]
    coords[:, 0] = coords[:, 0] - bounds
    coords[:, 1] = coords[:, 1] - bounds
    coords[:, 2] = coords[:, 2] - bounds
    
    return coords

#%% Load data
poly_data = np.load("polymer_coords.npy")
com = np.load("com.npy")

bounds = 40
#%% Wrap boundary conditions
coords = com[:, :, 200]
coords = unwrap_pbc(coords, 40)

#%% Visualize
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(coords)

pcd.paint_uniform_color([0, 0, 0])
o3d.visualization.draw_geometries([pcd])

#%% 

