#%% Imports and setup
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import open3d as o3d

from sklearn import cluster

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

#%% Functions
def unwrap_pbc_com(coords, bounds):    
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

def unwrap_pbc(coords, bounds):
    to_add = poly_coords.copy()
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                copied = to_add.copy()
                copied[:, :, 0] = copied[:, :, 0] + i * bounds
                copied[:, :, 1] = copied[:, :, 1] + j * bounds
                copied[:, :, 2] = copied[:, :, 2] + k * bounds
                coords = np.concatenate((coords, copied), axis = 0)
    
    coords = coords[np.shape(to_add)[0]:, :, :]
    coords[:, :, 0] = coords[:, :, 0] - bounds
    coords[:, :, 1] = coords[:, :, 1] - bounds
    coords[:, :, 2] = coords[:, :, 2] - bounds
    
    return coords

def get_clusters(com_coords, poly_coords, bounds):
    coords = unwrap_pbc_com(com_coords, bounds)
    poly = unwrap_pbc(poly_coords, bounds)
    
    cl = cluster.DBSCAN(eps = 2.5)
    cl.fit(coords)
    labels = cl.labels_
    
    cluster_com = np.zeros([labels.max() + 1, 3])

    for i in range(labels.max()):
        xyz = coords[labels == i, :]
        cluster_com[i, :] = np.sum(xyz, axis = 0)/np.shape(xyz)[0]

    x = cluster_com[:, 0]
    y = cluster_com[:, 1]
    z = cluster_com[:, 2]
    
    label_id = np.arange(0, labels.max() + 1, 1)
    
    valid_com = cluster_com[np.all([x > 0, x < bounds, y > 0, y < bounds,
                       z > 0, z < bounds], axis = 0), :]
    valid_labels = label_id[np.all([x > 0, x < bounds, y > 0, y < bounds,
                       z > 0, z < bounds], axis = 0)]
    
    cluster_dict = {}
    poly_dict = {}
    
    for i in range(len(valid_labels)):
        cluster_dict[valid_labels[i]] = coords[labels == valid_labels[i], :]
        poly_dict[valid_labels[i]] = poly[labels == valid_labels[i], :, :]
        
    return cluster_dict, poly_dict

def visualise_clusters(cluster_dict, bounds):  
    color_array = np.linspace(0, 240, num = len(cluster_dict)).astype(int)
    cmap = cm.viridis(color_array)
    
    to_plot = np.zeros([0, 3])
    colors = np.zeros([0, 3])
    i = 0
    for key, item in cluster_dict.items():
        to_plot = np.concatenate((to_plot, item), axis = 0)
        add_color = item.copy()
        add_color = add_color - add_color + cmap[i, 0:3]
        colors = np.concatenate((colors, add_color), axis = 0)
        i = i + 1
    
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(to_plot)
    pcd.colors = o3d.utility.Vector3dVector(colors)
    
    o3d.visualization.draw_geometries([pcd, get_bounding_box(bounds)])    
    
def visualise(coords, color):
    # Color arg is an RGB array
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(coords)
    
    pcd.paint_uniform_color(color)
    o3d.visualization.draw_geometries([pcd, get_bounding_box(bounds)])
    
def get_bounding_box(bounds):
    points = np.array([
        [0, 0, 0],
        [0, 0, bounds],
        [0, bounds, 0],
        [0, bounds, bounds],
        [bounds, 0, 0],
        [bounds, 0, bounds],
        [bounds, bounds, 0],
        [bounds, bounds, bounds],    
        ])
    
    lines = np.array([
        [0, 1],
        [0, 2],
        [0, 4],
        [1, 3],
        [1, 5],
        [2, 3],
        [2, 6],
        [3, 7],
        [4, 5],
        [4, 6],
        [5, 7],
        [6, 7],
        ])
    
    line_set = o3d.geometry.LineSet(
        points=o3d.utility.Vector3dVector(points),
        lines=o3d.utility.Vector2iVector(lines),
    )
    
    return line_set

#%% Load data
poly_data = np.load("polymer_coords.npy")
com = np.load("com.npy")

#%% Set Slicing Parameters
bounds = 40
step = 200

#%% Get cluster dictionaries
coords = com[:, :, step]
poly_coords = poly_data[:, :, :, step]

cluster_dict, poly_dict = get_clusters(coords, poly_coords, bounds)

#%% Plot number of clusters per time step
timesteps = 225
nclusters = np.zeros(timesteps)
steps = np.arange(0, timesteps, 1)
for i in range(timesteps):
    coords = com[:, :, i]
    poly_coords = poly_data[:, :, :, i]
    
    cluster_dict, poly_dict = get_clusters(coords, poly_coords, bounds)
    print("Timestep " + str(i))
    
    nclusters[i] = len(cluster_dict)
    
#%% Plotting
fig, ax = plt.subplots()
ax.scatter(steps, nclusters, marker = ".")

plt.figtext(.5,.95,
            "Number of Clusters Detected Against Timestep", 
            fontsize=14, ha='center')
plt.figtext(.5,.90,
            "Snapshot every 1000 Timesteps", 
            fontsize=12, ha='center')

ax.set_xlabel("Timestep ($\Delta t = 0.04$)", fontsize = 12)
ax.set_ylabel("Clusters", fontsize = 12)
plt.grid(b = True, which = "major", linestyle = "--", lw = 1)
ax.set_xlim([0, 250])

# fig.savefig('clusters.png', format='png', dpi=1200, bbox_inches='tight')

#%%

