import numpy as np
import open3d as o3d
import math
from matplotlib import cm
from sklearn import cluster

"""
Functions used so far for FYP project

To do:
    List functions by category
    Provide documentation for functions input and output!
"""
def movingaverage (x, y, lump = 1):
    entries = len(x)
    num_aves = math.floor(entries/lump)
    
    y_ave = np.zeros(num_aves)
    x_ave = np.zeros(num_aves)
    
    for i in range(num_aves):
        ave = np.sum(y[i * lump: (i+1) * lump]) / lump
        y_ave[i] = (np.sum(y_ave) + ave)/(i + 1)
        x_ave[i] = x[i * lump]
    
    return (x_ave, y_ave)

def get_rotation_matrix(a, b):
    """
    Get the rotation matrix via Rodrigues's method
    
    Arguments:
        a - [3] sized initial vector
        b - [3] sized target vector
    
    Returns:
        R - [3, 3] rotation matrix such that Ra = b
    """
    v = np.cross(a, b)
    v = v/np.linalg.norm(v)
    theta = np.arccos(np.dot(a, b))
    
    A = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    
    if theta < 0.001:
        R = np.identity(3)
    else:
        R = np.identity(3) + np.sin(theta) * A + \
            (1 - np.cos(theta)) * np.matmul(A, A)
    return R

def rotate_stuff(plane, coord, rotate_all = True):
    """
    Assumes coord is a [N, 9] array corresponding to the N particles with the
    quantities along the column as:
        x, y, z, pxx, pyy, pzz, pxy, pxz, pyz
        
    Returns the coordinates rotated such that the plane specified is parallel
    to a cardinal plane (xy, xz, xy)
    """    
    norm = plane[0:3]
    target = np.zeros(3)
    index, = np.where(abs(norm) == abs(norm).max())
    target[index] = 1
    
    R = get_rotation_matrix(norm, target)
    
    
    xyz = np.einsum("ij,kj -> ki", R, coord[:, 0:3])
    
    if not rotate_all:
        return xyz
    
    press_tensor = np.zeros([len(coord), 3, 3])
        
    for i in range(len(coord)):
        press_tensor[i, 0, 0] = coord[i, 3]
        press_tensor[i, 1, 1] = coord[i, 4]
        press_tensor[i, 2, 2] = coord[i, 5]
        
        press_tensor[i, 0, 1] = coord[i, 6]
        press_tensor[i, 0, 2] = coord[i, 7]
        press_tensor[i, 1, 2] = coord[i, 8]
        press_tensor[i, 2, 0] = coord[i, 7]
        press_tensor[i, 2, 1] = coord[i, 8]
        press_tensor[i, 1, 0] = coord[i, 6] 

    p = np.einsum("ip, jq, kpq -> kij", R, R, press_tensor)
    
    coord_trans = np.zeros([len(coord), 9])
    
    for i in range(len(coord)):
        coord_trans[i, 0] = xyz[i, 0]
        coord_trans[i, 1] = xyz[i, 1]
        coord_trans[i, 2] = xyz[i, 2]
        coord_trans[i, 3] = p[i, 0, 0]
        coord_trans[i, 4] = p[i, 1, 1]
        coord_trans[i, 5] = p[i, 2, 2]
        coord_trans[i, 6] = p[i, 0, 1]
        coord_trans[i, 7] = p[i, 0, 2]
        coord_trans[i, 8] = p[i, 1, 2]
    
    return coord_trans

def unwrap_polymers(first_slice, boxlengths):
    chain_length = np.shape(first_slice)[1]
    for i in range(np.shape(first_slice)[0]):
        for j in range(chain_length):
            for k in range(3):
                diff = first_slice[i, j, k] - first_slice[i, 0, k]
                if np.abs(diff) > boxlengths[k] * 0.7 and np.sign(diff) > 0:
                    first_slice[i, j, k] = first_slice[i, j, k] - boxlengths[k]
                elif np.abs(diff) > boxlengths[k] * 0.7 and np.sign(diff) < 0:
                    first_slice[i, j, k] = first_slice[i, j, k] + boxlengths[k]
    
    return first_slice

def get_poly_from_dump(dump, polyfile):
    polymer_num = np.load(polyfile).astype("int")
    polymer_coord = np.zeros(np.shape(polymer_num) + (3,1))
    
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1
    
    x_length = x_bounds[1] - x_bounds[0]
    y_length = y_bounds[1] - y_bounds[0]
    z_length = z_bounds[1] - z_bounds[0]
    
    boxlengths = np.array([x_length, y_length, z_length])    
    
    # Load 
    i = 0 # Just for debugging line info
    timestep = 0
    skip_count = 0
    skip = False
    snapshot_coord = np.zeros(np.shape(polymer_num) + (3,1))
    with open(dump) as infile:
        for line in infile:
            i = i + 1
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot_coord):
                    polymer_coord = np.concatenate((polymer_coord, snapshot_coord)
                                                   ,axis = 3)
                    snapshot_coord = np.zeros(np.shape(polymer_num) + (3,1))
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            if data[1] == 1:
                continue
            elif data[1] == 2 or data[1] == 3:
                row, col = np.where(polymer_num == data[0].astype("int"))
                snapshot_coord[row, col, 0] = data[2] * x_length
                snapshot_coord[row, col, 1] = data[3] * y_length
                snapshot_coord[row, col, 2] = data[4] * z_length
    
    polymer_coord = polymer_coord[:, :, :, 1:]
    
    for i in range(np.shape(polymer_coord)[3]):
        polymer_coord[:, :, :, i] = unwrap_polymers(polymer_coord[:, :, :, i], 
                                                    boxlengths)
        print(i)
        
    return polymer_coord

def get_allparticles_from_dump(dump):
    """
    Assumes that there are 3 unique particles inside!
    """
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1
    
    x_length = x_bounds[1] - x_bounds[0]
    y_length = y_bounds[1] - y_bounds[0]
    z_length = z_bounds[1] - z_bounds[0]
    boxlengths = np.array([x_length, y_length, z_length])

    i = 0 # Just for debugging line info
    timestep = 0
    skip_count = 0
    skip = False
    
    snapshot_coord1 = np.zeros((1, 3, 1))
    snapshot_coord2 = np.zeros((1, 3, 1))
    snapshot_coord3 = np.zeros((1, 3, 1))
    coord1 = np.zeros((1, 3, 1))
    coord2 = np.zeros((1, 3, 1))
    coord3 = np.zeros((1, 3, 1))
     
    with open(dump) as infile:
        for line in infile:
            i = i + 1
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot_coord1):              
                    snapshot_coord1 = snapshot_coord1[1:, :, :]
                    snapshot_coord2 = snapshot_coord2[1:, :, :]
                    snapshot_coord3 = snapshot_coord3[1:, :, :]
                    if not np.count_nonzero(coord1):
                        coord1 = snapshot_coord1
                        coord2 = snapshot_coord2
                        coord3 = snapshot_coord3
                    coord1 = np.concatenate((coord1, snapshot_coord1)
                                                   ,axis = 2)
                    coord2 = np.concatenate((coord2, snapshot_coord2)
                                                   ,axis = 2)
                    coord3 = np.concatenate((coord3, snapshot_coord3)
                                                   ,axis = 2)    
                    snapshot_coord1 = np.zeros((1, 3, 1))
                    snapshot_coord2 = np.zeros((1, 3, 1))
                    snapshot_coord3 = np.zeros((1, 3, 1))
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            if data[1] == 1:
                xyz = np.array([[[data[2] * x_length], [data[3] * y_length], 
                                [data[4] * z_length]]])
                snapshot_coord1 = np.concatenate((snapshot_coord1, xyz), 
                                                 axis = 0)
            elif data[1] == 2:
                xyz = np.array([[[data[2] * x_length], [data[3] * y_length], 
                                [data[4] * z_length]]])
                snapshot_coord2 = np.concatenate((snapshot_coord2, xyz), 
                                                 axis = 0)
            elif data[1] == 3:
                xyz = np.array([[[data[2] * x_length], [data[3] * y_length], 
                                [data[4] * z_length]]])
                snapshot_coord3 = np.concatenate((snapshot_coord3, xyz), 
                                                 axis = 0)
    return coord1, coord2, coord3
    
def get_allparticles_from_dump2(dump):
    """
    Assumes that there are 3 unique particles inside!
    Also assumes the column values after the xyz coords are the 6 pressure tensors!
    """
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1
    
    x_length = x_bounds[1] - x_bounds[0]
    y_length = y_bounds[1] - y_bounds[0]
    z_length = z_bounds[1] - z_bounds[0]
    boxlengths = np.array([x_length, y_length, z_length])

    i = 0 # Just for debugging line info
    timestep = 0
    skip_count = 0
    skip = False
    
    snapshot_coord1 = []
    snapshot_coord2 = []
    snapshot_coord3 = []
    coord1 = np.zeros((1, 9, 1))
    coord2 = np.zeros((1, 9, 1))
    coord3 = np.zeros((1, 9, 1))
     
    with open(dump) as infile:
        for line in infile:
            i = i + 1
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot_coord1):              
                    snapshot_coord1 = np.array(snapshot_coord1)
                    snapshot_coord2 = np.array(snapshot_coord2)
                    snapshot_coord3 = np.array(snapshot_coord3)
                    if not np.count_nonzero(coord1):
                        coord1 = snapshot_coord1
                        coord2 = snapshot_coord2
                        coord3 = snapshot_coord3
                    coord1 = np.concatenate((coord1, snapshot_coord1)
                                                   ,axis = 2)
                    coord2 = np.concatenate((coord2, snapshot_coord2)
                                                   ,axis = 2)
                    coord3 = np.concatenate((coord3, snapshot_coord3)
                                                   ,axis = 2)    
                    snapshot_coord1 = []
                    snapshot_coord2 = []
                    snapshot_coord3 = []
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            xyz = [[[data[2] * x_length], [data[3] * y_length],
                             [data[4] * z_length], [data[5]], [data[6]], 
                             [data[7]], [data[8]], [data[9]], [data[10]]]]
            if data[1] == 1:
                snapshot_coord1.append(xyz[0])
            elif data[1] == 2:
                snapshot_coord2.append(xyz[0])
            elif data[1] == 3:
                snapshot_coord3.append(xyz[0])
                
    return coord1, coord2, coord3
    
def get_allparticles_from_dump3(dump):
    """
    Assumes that there are 5 unique particles inside!
    """
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1

        x_length = x_bounds[1] - x_bounds[0]
        y_length = y_bounds[1] - y_bounds[0]
        z_length = z_bounds[1] - z_bounds[0]
        boxlengths = np.array([x_length, y_length, z_length])

        i = 0 # Just for debugging line info
        timestep = 0
        skip_count = 0
        skip = False

        snapshot_coord1 = []
        snapshot_coord2 = []
        snapshot_coord3 = []
        snapshot_coord4 = []
        snapshot_coord5 = []
        coord1 = np.zeros((1, 3, 1))
        coord2 = np.zeros((1, 3, 1))
        coord3 = np.zeros((1, 3, 1))
        coord5 = np.zeros((1, 3, 1))
        coord5 = np.zeros((1, 3, 1))

    with open(dump) as infile:
        for line in infile:
            i = i + 1
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot_coord1):              
                    snapshot_coord1 = np.array(snapshot_coord1)
                    snapshot_coord2 = np.array(snapshot_coord2)
                    snapshot_coord3 = np.array(snapshot_coord3)
                    snapshot_coord4 = np.array(snapshot_coord4)
                    snapshot_coord5 = np.array(snapshot_coord5)
                    if not np.count_nonzero(coord1):
                        coord1 = snapshot_coord1
                        coord2 = snapshot_coord2
                        coord3 = snapshot_coord3
                        coord4 = snapshot_coord4
                        coord5 = snapshot_coord5
                    coord1 = np.concatenate((coord1, snapshot_coord1)
                                                   ,axis = 2)
                    coord2 = np.concatenate((coord2, snapshot_coord2)
                                                   ,axis = 2)
                    coord3 = np.concatenate((coord3, snapshot_coord3)
                                                   ,axis = 2)
                    coord4 = np.concatenate((coord4, snapshot_coord4)
                                                   ,axis = 2) 
                    coord5 = np.concatenate((coord5, snapshot_coord5)
                                                   ,axis = 2) 
                    snapshot_coord1 = []
                    snapshot_coord2 = []
                    snapshot_coord3 = []
                    snapshot_coord4 = []
                    snapshot_coord5 = []
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            xyz = [[[data[2] * x_length], [data[3] * y_length],
                             [data[4] * z_length]]]
            if data[1] == 1:
                snapshot_coord1.append(xyz[0])
            elif data[1] == 2:
                snapshot_coord2.append(xyz[0])
            elif data[1] == 3:
                snapshot_coord3.append(xyz[0])
            elif data[1] == 4:
                snapshot_coord4.append(xyz[0])
            elif data[1] == 5:
                snapshot_coord5.append(xyz[0])
                
    return coord1, coord2, coord3, coord4, coord5

def get_data_from_dump(dump, num_particles = 3, num_variables = 0):
    
    # Load the dump file line by line, read input parameters
    i = 0
    with open(dump) as infile:
        for line in infile:
            if i == 3:
                tot_particles = np.fromstring(line, dtype = "int", sep = " ")
            elif i == 5:
                x_bounds = np.fromstring(line, sep = " ")
            elif i == 6:
                y_bounds = np.fromstring(line, sep = " ")
            elif i == 7:
                z_bounds = np.fromstring(line, sep = " ")
            elif i == 8:
                break
            i = i + 1
    
        x_length = x_bounds[1] - x_bounds[0]
        y_length = y_bounds[1] - y_bounds[0]
        z_length = z_bounds[1] - z_bounds[0]
        boxlengths = np.array([x_length, y_length, z_length])

    timestep = 0
    skip_count = 0
    skip = False
    
    snapshot = {}
    coord = {}
    timestep_array = []
    
    for i in range(num_particles):
        snapshot[i + 1] = []
        coord[i + 1] = np.zeros((1, num_variables + 4, 1))

    with open(dump) as infile:
        for line in infile:
            if line == "ITEM: TIMESTEP\n":
                skip = True
                skip_count = 0
                print("Snapshot: " + str(timestep))
                if np.count_nonzero(snapshot[1]):
                    for i in range(num_particles):
                        snapshot[i + 1] = np.array(snapshot[i + 1])
                        
                    if not np.count_nonzero(coord[1]):
                        for i in range(num_particles):
                            coord[i + 1] = snapshot[i + 1]
                    
                    for i in range(num_particles):
                        coord[i + 1] = np.concatenate((coord[i + 1], 
                                                       snapshot[i + 1]), 
                                                      axis = 2)
                        snapshot[i + 1] = []
    
                timestep = timestep + 1
                continue
            elif skip:
                skip_count = skip_count + 1
                if skip_count == 1:
                    timestep_array += [int(line)]
                elif skip_count == 8:
                    skip = False
                continue
            
            data = np.fromstring(line, sep = " ")
            xyz = [[data[2] * x_length], [data[3] * y_length],
                             [data[4] * z_length], [data[0]]]
            for i in range(num_variables):
                xyz.append([data[-1 - i]])
            
            snapshot[data[1]].append(xyz)

    timestep_array = np.array(timestep_array)
    timesteps, index = np.unique(timestep_array, return_index = True)
    for i in range(num_particles):
        coord[i + 1] = coord[i + 1][:, :, index]
    
    return timesteps, coord

def get_e2e(first_slice):
    e2e = np.zeros(np.shape(first_slice)[0])
    
    for i in range(np.shape(first_slice)[0]):
        vec = first_slice[i, 0, :] - first_slice[i, -1, :]
        e2e[i] = np.linalg.norm(vec)
    
    return e2e

def find_com(first_slice, polymer_struct, polymer_masses):
    totmass = 0
    for i in range(len(polymer_struct)):
        totmass += polymer_masses[polymer_struct[i] - 1]
    
    com = np.zeros((np.shape(first_slice)[0], 3))

    for i in range(np.shape(first_slice)[0]):
        denom = np.array([0, 0, 0])
        for j in range(np.shape(first_slice)[1]):
            denom = denom + polymer_masses[polymer_struct[j] - 1] * \
                first_slice[i, j, :]
        
        com[i, :] = denom/totmass
    
    return com

def find_r_g(first_slice, com_slice):
    r_g = np.zeros(np.shape(first_slice)[0])

    for i in range(np.shape(first_slice)[0]):
        vec = 0
        for j in range(np.shape(first_slice)[1]):
            vec = vec + np.linalg.norm(first_slice[i, j, :] - \
                                       com_slice[i, :]) ** 2
        
        r_g[i] = (vec/np.shape(first_slice)[1]) ** 0.5
    
    return r_g

def get_com_all(poly_data, polymer_struct, 
                polymer_masses = np.array([1, 1, 1])):
    
    com = np.zeros((np.shape(poly_data)[0], 3, np.shape(poly_data)[3]))
    for i in range(np.shape(poly_data)[3]):
        com[:, :, i] = find_com(poly_data[:, :, :, i], polymer_struct,
                                polymer_masses)
        
    return com

def get_r_g_all(poly_data, com):
    r_g = np.zeros((np.shape(poly_data)[0], np.shape(poly_data)[3]))
    for i in range(np.shape(poly_data)[3]):
        r_g[:, i] = find_r_g(poly_data[:, :, :, i], com[:, :, i])
        
    return r_g

def get_r_g_stats(r_g):    
    r_g_ave = np.zeros(np.shape(r_g)[1])
    r_g_std = np.zeros(np.shape(r_g)[1])
    dt = np.zeros(np.shape(r_g)[1])
    for i in range(np.size(r_g_ave)):
        r_g_ave[i] = np.average(r_g[:, i])
        r_g_std[i] = np.std(r_g[:, i])
        dt[i] = i * 500
    
    return r_g_ave, r_g_std
    
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
    to_add = coords.copy()
    
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

def get_clusters(com_coords, poly_coords, bounds, cluster_dist = 2.5):
    coords = unwrap_pbc_com(com_coords, bounds)
    poly = unwrap_pbc(poly_coords, bounds)
    
    cl = cluster.DBSCAN(eps = cluster_dist)
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
    cmap = cm.plasma(color_array)
    
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
    
def visualise(coords, color, bounds):
    # Color arg is an RGB array
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(coords)
    
    pcd.paint_uniform_color(color)
    o3d.visualization.draw_geometries([pcd, get_bounding_box(bounds)])
    
def visualise_plane_with_points(plane, bounds, point_number, coords):
    plane_points = generate_plane_points(bounds, plane, point_number)
    
    plane_color = np.zeros([len(plane_points), 3])
    plane_color[:] = plane_color[:] + np.array([0, 0, 1])
    
    point_color = np.zeros([len(coords), 3])
    point_color[:] = point_color[:] + np.array([0, 0, 0])
    
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(np.concatenate((plane_points, 
                                                            coords)))
    pcd.colors = o3d.utility.Vector3dVector(np.concatenate((plane_color, 
                                                            point_color)))

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

def generate_plane_points(bounds, plane, num_points):
    A = plane[0]
    B = plane[1]
    C = plane[2]
    D = plane[3]
    
    x = np.linspace(0, bounds, num_points)
    y = np.linspace(0, bounds, num_points)
    
    gridx, gridy = np.meshgrid(x, y)
    grid = np.array([gridx, gridy]).reshape([2, num_points**2]).T
    
    z = (D - grid[:, 0] * A - grid[:, 1] * B)/C
    
    x = gridx.reshape([num_points**2]).T
    y = gridy.reshape([num_points**2]).T
    
    return np.array([x, y, z]).T
    
def create_poly_datafile(box_size, polymer_conc, polymer_structure, 
                         filename = "data.txt", mass_array = [1, 44/54, 56/54], 
                         rho_star = 3):
    def got_space(x, y, z, matrix, polymer_length):
        space = True
        for i in range(polymer_length):
            if matrix[x + i][y][z] != 0:
                space = False
        
        return space
    
    
    polymer_length = len(polymer_structure)
    polymer_structure = polymer_structure.astype(int)
    description = str(polymer_length) + "_" + str(polymer_conc)
    
    angles = 0
    dihedrals = 0
    impropers = 0
    
    atom_types = len(mass_array)
    bond_types = 1
    angle_types = 0
    dihedral_types = 0
    improper_types = 0
    
    # Calculate simulation parameters
    particle_count = rho_star * box_size**3
    
    polymer_beadcount = np.round((particle_count * polymer_conc)).astype(int)
    num_chains = np.round((polymer_beadcount/polymer_length)).astype(int)
    polymer_beadcount = num_chains * polymer_length
    
    num_solvent = particle_count - polymer_beadcount
    
    gridsize = np.ceil((particle_count ** (1/3))).astype(int)
    spacing = box_size / (gridsize - 1)
    
    matrix = np.zeros((gridsize, gridsize, gridsize))
    polymer_coords = np.zeros((num_chains, 3))
    # Set the coordinates for the polymer heads and other chains
    # Polymer chains are initially all oriented in the x-direction
    
    """
    Sets the polymer bead type array such that the head has value of 69, such that 
    when filtering the matrix values in writing to input file, the polymer head 
    will not be mis-attributed
    
    e.g. when polymer_structure = [2, 2, 3, 3]
    
    """
    polymer_encoding = polymer_structure.copy()
    polymer_encoding[0] = 69
    
    
    i = 0 
    while i < num_chains:
        x = np.random.randint(0, gridsize - polymer_length)
        y = np.random.randint(0, gridsize - 1)
        z = np.random.randint(0, gridsize - 1)
        
        if got_space(x, y, z, matrix, polymer_length):
            for j in range(polymer_length):
                matrix[x + j][y][z] = polymer_encoding[j]
                
            polymer_coords[i][0] = x
            polymer_coords[i][1] = y
            polymer_coords[i][2] = z
            
            i = i + 1
            
    # Populate remaining with solvent beads
    i = 0
    for x in range(gridsize):
        for y in range(gridsize):
            for z in range(gridsize):
                if (matrix[x][y][z] == 0) and (i < num_solvent):
                    matrix[x][y][z] = 1
                    i = i + 1
    
    
    with open(filename, mode = "w") as line:
        line.write("LAMMPS " + description + "\n \n")
        
        line.write(str(particle_count) + " atoms\n")
        line.write(str(num_chains * (polymer_length - 1)) + " bonds\n")
        line.write(str(angles) + " angles\n")
        line.write(str(dihedrals) + " dihedrals\n")
        line.write(str(impropers) + " impropers\n \n")
        
        line.write(str(atom_types) + " atom types\n")
        line.write(str(bond_types) + " bond types\n")
        line.write(str(angle_types) + " angle types\n")
        line.write(str(dihedral_types) + " dihedral types\n")
        line.write(str(improper_types) + " improper types\n \n")
        
        line.write("0 " + str(box_size) + " xlo xhi\n")
        line.write("0 " + str(box_size) + " ylo yhi\n")
        line.write("0 " + str(box_size) + " zlo zhi\n \n")
        
        line.write("Masses\n \n")
        for index, val in enumerate(mass_array):
            line.write(str(index + 1) + " " + str(val) + "\n")
            
        line.write("\n")
        
        line.write("Atoms\n \n")
        i = 1
        j = 0
        polymer_number = np.zeros((num_chains, polymer_length))
        for index, val in np.ndenumerate(matrix):
            x = index[0] * spacing
            y = index[1] * spacing
            z = index[2] * spacing
            if val == 1:
                line.write(str(i) + " 0 1 " + 
                           "{:.3f}".format(x) + " " + "{:.3f}".format(y) + " " 
                           + "{:.3f}".format(z) + "\n")
                i = i + 1
            elif val == 69:
                for k in range(polymer_length):
                     line.write(str(i) + " 0 " + str(polymer_structure[k]) + 
                                " " + "{:.3f}".format(x + spacing * k) + " " + 
                           "{:.3f}".format(y) + " " + "{:.3f}".format(z) + 
                           "\n")
                     polymer_number[j][k] = i
                     i = i + 1
                j = j + 1
                
            
        line.write("\n")
        
        line.write("Bonds\n \n")
        i = 1
        bond_per_chain = polymer_length - 1
        polymer_number = polymer_number.astype(int)
        for row in polymer_number:
            for j in range(bond_per_chain):
                line.write(str(i) + " 1 " + str(row[j]) + " " + str(row[j + 1])
                           + "\n")
                i = i + 1
                
    return polymer_number

def create_bilayer_datafile(box_size, 
                            polymer_conc, 
                            polymer_head, 
                            polymer_t1,
                            polymer_t2, 
                            filename = "data.txt", 
                            mass_array = [1, 1, 1], 
                            rho_star = 3, 
                            description = "Nothing"):
    
    def got_space_poly(x, y, z, matrix, poly_num):
        space = True
        for i in range(poly_num[0]):
            if matrix[x + i][y][z] != 0:
                space = False        
        for i in range(poly_num[1]):
            if matrix[x + poly_num[0] + i][y][z] != 0:
                space = False        
        for i in range(poly_num[2]):
            if matrix[x + poly_num[0] + i][y + 1][z] != 0:
                space = False        
        return space
        
    def get_angle(poly_num):
        count = 0
        for i in poly_num:
            if i > 2:
                count += i - 2        
        # Account for no tails, only head (AKA straight line polymer)
        if len(poly_num) != 1:
            count += 3            
            if poly_num[1] > 1:
                count += 1                
            if poly_num[2] > 1:
                count += 1        
        return count
    
    poly_struct = {}
    poly_struct["head"] = polymer_head.astype(int)
    poly_struct["t1"] = polymer_t1.astype(int)
    poly_struct["t2"] = polymer_t2.astype(int)    
    poly_num = np.array([len(polymer_head), len(polymer_t1), len(polymer_t2)])    
    polymer_length = 0
    
    for key, val in poly_struct.items():
        polymer_length += len(val)
    
    dihedrals = 0
    impropers = 0    
    atom_types = len(mass_array)
    bond_types = 1
    angle_types = 3
    dihedral_types = 0
    improper_types = 0
    
    # Calculate simulation parameters
    particle_count = rho_star * box_size**3   
    polymer_beadcount = np.round((particle_count * polymer_conc)).astype(int)
    num_chains = np.round((polymer_beadcount/polymer_length)).astype(int)
    polymer_beadcount = num_chains * polymer_length
    
    num_solvent = particle_count - polymer_beadcount
    angles = num_chains * get_angle(poly_num)
    bonds = num_chains * (polymer_length - 1)
    
    gridsize = np.ceil((particle_count ** (1/3))).astype(int)
    spacing = box_size / (gridsize - 1)
    
    # Create matrix to populate with particles
    matrix = np.zeros((gridsize, gridsize, gridsize))
    polymer_coords = np.zeros((num_chains, 3))
    
    # Populate the matrix with polymer beads first
    i = 0
    while i < num_chains:
        x = np.random.randint(0, gridsize - polymer_length)
        y = np.random.randint(0, gridsize - 1)
        z = np.random.randint(0, gridsize - 1)
        
        if got_space_poly(x, y, z, matrix, poly_num):
            for j in range(poly_num[0]):
                matrix[x + j][y][z] = poly_struct["head"][j]            
            for j in range(poly_num[1]):
                matrix[x + poly_num[0] + j][y][z] = poly_struct["t1"][j]            
            for j in range(poly_num[2]):
                matrix[x + poly_num[0] + j][y + 1][z] = poly_struct["t2"][j]            
            i = i + 1            
            matrix[x][y][z] = 69
    
    # Populate remaining with solvent beads
    i = 0
    for x in range(gridsize):
        for y in range(gridsize):
            for z in range(gridsize):
                if (matrix[x][y][z] == 0) and (i < num_solvent):
                    matrix[x][y][z] = 1
                    i = i + 1
    
    # Write to LAMMPS
    with open(filename, mode = "w") as line:
        line.write("LAMMPS " + description + "\n \n")
        
        line.write(str(particle_count) + " atoms\n")
        line.write(str(bonds) + " bonds\n")
        line.write(str(angles) + " angles\n")
        line.write(str(dihedrals) + " dihedrals\n")
        line.write(str(impropers) + " impropers\n \n")
        
        line.write(str(atom_types) + " atom types\n")
        line.write(str(bond_types) + " bond types\n")
        line.write(str(angle_types) + " angle types\n")
        line.write(str(dihedral_types) + " dihedral types\n")
        line.write(str(improper_types) + " improper types\n \n")
        
        line.write("0 " + str(box_size) + " xlo xhi\n")
        line.write("0 " + str(box_size) + " ylo yhi\n")
        line.write("0 " + str(box_size) + " zlo zhi\n \n")
        
        line.write("Masses\n \n")
        for index, val in enumerate(mass_array):
            line.write(f"{index + 1} {val}\n")
            
        line.write("\n")
        
        line.write("Atoms\n \n")
        i = 1
        j = 0
        polymer_number = np.zeros((num_chains, polymer_length), dtype = int)
        
        for index, val in np.ndenumerate(matrix):
            x = index[0] * spacing
            y = index[1] * spacing
            z = index[2] * spacing
            if val == 1:
                line.write(f"{i} 0 1 {x:.3f} {y:.3f} {z:.3f}\n")
                i = i + 1
            elif val == 69:
                for k in range(poly_num[0]):
                    bead_type = poly_struct["head"][k]
                    line.write(f"{i} 0 {bead_type} "
                               f"{(x + spacing * k):.3f} {y:.3f} {z:.3f}\n")
                    polymer_number[j][k] = i
                    i = i + 1                   
                for k in range(poly_num[1]):
                    bead_type = poly_struct["t1"][k]
                    line.write(f"{i} 0 {bead_type} "
                               f"{(x + spacing * (k + poly_num[0])):.3f} "
                               f"{y:.3f} {z:.3f}\n")
                    polymer_number[j][k + poly_num[0]] = i
                    i = i + 1                
                for k in range(poly_num[2]):
                    bead_type = poly_struct["t2"][k]
                    line.write(f"{i} 0 {bead_type} "
                               f"{(x + spacing * (k + poly_num[0])):.3f} "
                               f"{(y + spacing):.3f} {z:.3f}\n")
                    polymer_number[j][k + poly_num[0] + poly_num[1]] = i
                    i = i + 1
                j = j + 1
                
        line.write("\n")
        line.write("Bonds\n \n")
        
        i = 1
        for row in polymer_number:
            if poly_num[0] != 1:
                for j in range(poly_num[0] - 1):
                    line.write(f"{i} 1 {row[j]} {row[j + 1]}\n")
                    i = i + 1            
            if poly_num[1] != 1:
                for j in range(poly_num[1] - 1):
                    line.write(f"{i} 1 {row[j + poly_num[0]]} "
                               f"{row[j + 1 + poly_num[0]]}\n")
                    i = i + 1                
            if poly_num[2] != 1:
                for j in range(poly_num[2] - 1):
                    line.write(f"{i} 1 {row[j + poly_num[0] + poly_num[1]]} "
                               f"{row[j + 1 + poly_num[0] + poly_num[1]]}\n")
                    i = i + 1
            
            line.write(f"{i} 1 {row[poly_num[0] - 1]} {row[poly_num[0]]}\n")
            i = i + 1
            line.write(f"{i} 1 {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
        
        line.write("\n")
        line.write("Angles\n \n")
        
        i = 1
        for row in polymer_number:
            if poly_num[0] > 2:
                for j in range(poly_num[0] - 2):
                    line.write(f"{i} 1 {row[j]} {row[j + 1]} {row[j + 2]}\n")
                    i = i + 1            
            if poly_num[1] > 2:
                for j in range(poly_num[1] - 2):
                    line.write(f"{i} 1 {row[j + poly_num[0]]} "
                               f"{row[j + poly_num[0] + 1]} "
                               f"{row[j + poly_num[0] + 2]}\n")
                    i = i + 1            
            if poly_num[2] > 2:
                for j in range(poly_num[2] - 2):
                    line.write(f"{i} 1 {row[j + poly_num[0] + poly_num[1]]} "
                               f"{row[j + poly_num[0] + poly_num[1] + 1]} "
                               f"{row[j + poly_num[0] + poly_num[1] + 2]}\n")
                    i = i + 1
            
            line.write(f"{i} 2 {row[poly_num[0] - 2]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0]]}\n")
            i = i + 1
            line.write(f"{i} 2 {row[poly_num[0] - 2]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
            line.write(f"{i} 2 {row[poly_num[0]]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
            
            if poly_num[1] > 1:
                line.write(f"{i} 3 {row[poly_num[0] - 1]} {row[poly_num[0]]} "
                           f"{row[poly_num[0] + 1]}\n")
                i = i + 1            
            if poly_num[2] > 1:
                line.write(f"{i} 3 {row[poly_num[0] - 1]} "
                           f"{row[poly_num[0] + poly_num[1]]} "
                           f"{row[poly_num[0] + poly_num[1] + 1]}\n")
                i = i + 1
    return polymer_number

def create_bilayer_datafile_v2(box_size,
                               polymer_head,
                               polymer_t1,
                               polymer_t2,
                               dense,
                               filename = "data.txt",
                               mass_array = [1, 1, 1],
                               rho_star = 3,
                               description = "Nothing"):
    
    def got_space_poly(x, y, z, matrix, poly_num):
        space = True
        for i in range(poly_num[0]):
            if matrix[x + i][y][z] != 0:
                space = False        
        for i in range(poly_num[1]):
            if matrix[x + poly_num[0] + i][y][z] != 0:
                space = False        
        for i in range(poly_num[2]):
            if matrix[x + poly_num[0] + i][y + 1][z] != 0:
                space = False        
        return space
    
    def get_angle(poly_num):
        count = 0
        for i in poly_num:
            if i > 2:
                count += i - 2        
        # Account for no tails, only head (AKA straight line polymer)
        if len(poly_num) != 1:
            count += 3            
            if poly_num[1] > 1:
                count += 1                
            if poly_num[2] > 1:
                count += 1        
        return count

    poly_struct = {}
    poly_struct["head"] = polymer_head.astype(int)
    poly_struct["t1"] = polymer_t1.astype(int)
    poly_struct["t2"] = polymer_t2.astype(int)    
    poly_num = np.array([len(polymer_head), len(polymer_t1), len(polymer_t2)])    
    polymer_length = 0
    
    for key, val in poly_struct.items():
        polymer_length += len(val)

    dihedrals = 0
    impropers = 0    
    atom_types = len(mass_array)
    bond_types = 1
    angle_types = 3
    dihedral_types = 0
    improper_types = 0
    
    # Calculate simulation parameters
    print("Calculating Simulation Parameters...")
    particle_count = rho_star * box_size**3   
    
    gridsize = np.ceil((particle_count ** (1/3))).astype(int)
    spacing = box_size / (gridsize - 1)

    bilayer_top_chains = np.floor((gridsize * gridsize * dense)/2).astype(int)
    num_chains = np.floor(bilayer_top_chains * 2 * 1.1).astype(int)
    
    poly_index = np.zeros([num_chains, polymer_length, 3])
    
    head_thickness = np.ceil(len(polymer_head)/2).astype(int)
    tail_thickness = max(len(polymer_t1), len(polymer_t2))

    height = (gridsize/4).astype(int)
    even = True if len(polymer_head) % 2 == 0 else False

    # Create matrix to populate with particles
    print("Populating Space with Particles...")
    matrix = np.zeros((gridsize, gridsize, gridsize))
    
    cc = 0
    for i in range(gridsize):
        for j in range(0, gridsize, 2):
            if np.random.uniform() > dense:
                continue
        
            pc = 0
            for k in range(head_thickness):
                z = k + height
                z1 = height + 2 * (head_thickness + tail_thickness) - k - 1
                if k == 0 and even:
                    poly_index[cc, pc, :] = np.array([i, j, z])
                    matrix[i, j, z] = 2
                    poly_index[cc + 1, pc, :] = np.array([i, j, z1])
                    matrix[i, j, z1] = 2
                    
                    pc += 1
                    
                    poly_index[cc, pc, :] = np.array([i, j + 1, z])
                    matrix[i, j + 1, z] = 2
                    poly_index[cc + 1, pc, :] = np.array([i, j + 1, z1])
                    matrix[i, j + 1, z1] = 2                    
                    
                    pc += 1
                    
                elif k == 0:
                    poly_index[cc, pc, :] = np.array([i, j, z])
                    matrix[i, j, z] = 2
                    poly_index[cc + 1, pc, :] = np.array([i, j, z1])
                    matrix[i, j, z1] = 2
                    
                    pc += 1                    
                else:
                    poly_index[cc, pc, :] = np.array([i, j, z])
                    matrix[i, j, z] = 2
                    poly_index[cc + 1, pc, :] = np.array([i, j, z1])
                    matrix[i, j, z1] = 2
                    pc += 1
                    
                    poly_index[cc, pc, :] = np.array([i, j + 1, z])
                    matrix[i, j + 1, z] = 2
                    poly_index[cc + 1, pc, :] = np.array([i, j + 1, z1])
                    matrix[i, j + 1, z1] = 2 
                    pc += 1
                    
            for k in range(tail_thickness):
                z = k + height + head_thickness
                z1 = height + 2 * (tail_thickness)  + head_thickness - k - 1
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 3
                poly_index[cc + 1, pc, :] = np.array([i, j, z1])
                matrix[i, j, z1] = 3
                pc += 1
            
            
            for k in range(tail_thickness):
                z = k + height + head_thickness
                z1 = height + 2 * (tail_thickness)  + head_thickness - k - 1
                poly_index[cc, pc, :] = np.array([i, j + 1, z])
                matrix[i, j + 1, z] = 3
                poly_index[cc + 1, pc, :] = np.array([i, j + 1, z1])
                matrix[i, j + 1, z1] = 3
                pc += 1
                
            cc += 2
            pc = 0
    
    num_chains = cc
    angles = num_chains * get_angle(poly_num)
    bonds = num_chains * (polymer_length - 1)

    poly_index = poly_index[0:cc, :, :]

    num_solvent = particle_count - num_chains * polymer_length
    polymer_conc = (num_chains * polymer_length)/particle_count

    # Populate remaining with solvent beads
    i = 0
    for x in range(gridsize):
        for y in range(gridsize):
            for z in range(gridsize):
                if (matrix[x][y][z] == 0) and (i < num_solvent):
                    matrix[x][y][z] = 1
                    i = i + 1
                    
    print("Writing to File...")
    with open(filename, mode = "w") as line:
        line.write("LAMMPS " + description + "\n \n")
            
        line.write(str(particle_count) + " atoms\n")
        line.write(str(bonds) + " bonds\n")
        line.write(str(angles) + " angles\n")
        line.write(str(dihedrals) + " dihedrals\n")
        line.write(str(impropers) + " impropers\n \n")
        
        line.write(str(atom_types) + " atom types\n")
        line.write(str(bond_types) + " bond types\n")
        line.write(str(angle_types) + " angle types\n")
        line.write(str(dihedral_types) + " dihedral types\n")
        line.write(str(improper_types) + " improper types\n \n")
        
        line.write("0 " + str(box_size) + " xlo xhi\n")
        line.write("0 " + str(box_size) + " ylo yhi\n")
        line.write("0 " + str(box_size) + " zlo zhi\n \n")
        
        line.write("Masses\n \n")
        for index, val in enumerate(mass_array):
            line.write(f"{index + 1} {val}\n")
            
        line.write("\n")
        
        line.write("Atoms\n \n")
        
        i = 1
        polymer_number = np.zeros((num_chains, polymer_length), dtype = int)
        
        for index, val in np.ndenumerate(matrix):
            x = index[0] * spacing
            y = index[1] * spacing
            z = index[2] * spacing
            
            if val == 1:
                line.write(f"{i} 0 {val:.0f} {x:.3f} {y:.3f} {z:.3f}\n")
                i += 1
            elif val == 2 or val == 3:
                line.write(f"{i} 0 {val:.0f} {x:.3f} {y:.3f} {z:.3f}\n")
                trutharray = np.all([poly_index[:, :, 0] == index[0],
                                     poly_index[:, :, 1] == index[1],
                                     poly_index[:, :, 2] == index[2]], axis = 0)
                cc, pc = np.where(trutharray == True)
                
                polymer_number[cc, pc] = i
                i += 1
         
        line.write("\n")
        line.write("Bonds\n \n")
        
        i = 1
        for row in polymer_number:
            if poly_num[0] != 1:
                for j in range(poly_num[0] - 1):
                    line.write(f"{i} 1 {row[j]} {row[j + 1]}\n")
                    i = i + 1            
            if poly_num[1] != 1:
                for j in range(poly_num[1] - 1):
                    line.write(f"{i} 1 {row[j + poly_num[0]]} "
                               f"{row[j + 1 + poly_num[0]]}\n")
                    i = i + 1                
            if poly_num[2] != 1:
                for j in range(poly_num[2] - 1):
                    line.write(f"{i} 1 {row[j + poly_num[0] + poly_num[1]]} "
                               f"{row[j + 1 + poly_num[0] + poly_num[1]]}\n")
                    i = i + 1
            
            line.write(f"{i} 1 {row[poly_num[0] - 1]} {row[poly_num[0]]}\n")
            i = i + 1
            line.write(f"{i} 1 {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
        
        line.write("\n")
        line.write("Angles\n \n")
        
        i = 1
        for row in polymer_number:
            if poly_num[0] > 2:
                for j in range(poly_num[0] - 2):
                    line.write(f"{i} 1 {row[j]} {row[j + 1]} {row[j + 2]}\n")
                    i = i + 1            
            if poly_num[1] > 2:
                for j in range(poly_num[1] - 2):
                    line.write(f"{i} 1 {row[j + poly_num[0]]} "
                               f"{row[j + poly_num[0] + 1]} "
                               f"{row[j + poly_num[0] + 2]}\n")
                    i = i + 1            
            if poly_num[2] > 2:
                for j in range(poly_num[2] - 2):
                    line.write(f"{i} 1 {row[j + poly_num[0] + poly_num[1]]} "
                               f"{row[j + poly_num[0] + poly_num[1] + 1]} "
                               f"{row[j + poly_num[0] + poly_num[1] + 2]}\n")
                    i = i + 1
            
            line.write(f"{i} 2 {row[poly_num[0] - 2]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0]]}\n")
            i = i + 1
            line.write(f"{i} 2 {row[poly_num[0] - 2]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
            line.write(f"{i} 2 {row[poly_num[0]]} {row[poly_num[0] - 1]} "
                       f"{row[poly_num[0] + poly_num[1]]}\n")
            i = i + 1
            
            if poly_num[1] > 1:
                line.write(f"{i} 3 {row[poly_num[0] - 1]} {row[poly_num[0]]} "
                           f"{row[poly_num[0] + 1]}\n")
                i = i + 1            
            if poly_num[2] > 1:
                line.write(f"{i} 3 {row[poly_num[0] - 1]} "
                           f"{row[poly_num[0] + poly_num[1]]} "
                           f"{row[poly_num[0] + poly_num[1] + 1]}\n")
                i = i + 1
    
    return polymer_number, polymer_conc, num_chains

def get_planes(com_poly, resolution, bounds, top_number):
    def hough_transform(x, y, z, theta, phi):
        rho = x * np.cos(theta) * np.sin(phi)
        rho += y * np.sin(theta) * np.sin(phi)
        rho += z * np.cos(phi)
        
        return rho
    
    def hough_reverse(theta, phi, rho):
        A = np.cos(theta) * np.sin(phi)
        B = np.sin(theta) * np.sin(phi)
        C = np.cos(phi)
        D = rho
        
        return np.array([A, B, C, D])
    
    res_rho = resolution * 3
    theta = np.linspace(0, np.pi * 2, resolution)
    phi = np.linspace(-0.5 * np.pi, 0.5 * np.pi, resolution)
    
    theta_grid, phi_grid = np.meshgrid(theta, phi)
    big_grid = np.array([theta_grid, phi_grid]).reshape([2, resolution**2]).T
    
    rho = hough_transform(bounds, bounds, bounds, big_grid[:, 0], 
                          big_grid[:, 1])
    max_rho = rho.max() + 5
    min_rho = rho.min() + 5
    rho_range = max_rho - min_rho
    rho = np.linspace(min_rho, max_rho, res_rho)
    accum = np.zeros([resolution * resolution, res_rho])    
    
    for i in range(len(com_poly)):
        x = com_poly[i, 0]
        y = com_poly[i, 1]
        z = com_poly[i, 2]
        output = hough_transform(x, y, z, big_grid[:, 0], big_grid[:, 1])
        index = ((output - min_rho)/rho_range) * res_rho
        index = np.trunc(index).astype(int)
        for j in range(len(index)):
            accum[j, index[j]] = accum[j, index[j]] + 1
        
        if i % 10 == 0:
            print(f"Polymer count: {i}")
    
    accum = accum.reshape([resolution, resolution, res_rho], order = "F")
    accum_sorted = accum.flatten()
    accum_sorted = np.flip(np.sort(accum_sorted))
    
    planes = np.zeros([top_number, 4])
    for i in range(top_number):
        theta1, phi1, rho1 = np.where(accum == accum_sorted[i])
        planes[i] = hough_reverse(theta[theta1[0]], phi[phi1[0]], 
                                  rho[rho1[0]])
        
    return planes
    
def get_bilayer_histogram(snap1, snap2, snap3, num_bins, plane, bound, 
                           density = 3):
    """
    Gets the histogram of particle concentration in direction of the normal to
    the plane of the bilayer
    
    Arguments:
        snap1, snap2, snap3 - [N,3] coordinates of each particle type for a 
                                certain snapshot in time
        num_bins - Number of divisions
        plane - Array of the plane coefficients [A, B, C, D]
        bound - Bounds of the simulation, assuming a cube
    
    Returns:
        c1, c2, c3 - [num_bins] sized histogram array of the particle density
        axis - [num_bins] sized array denoting the bin edges
    """
    def get_distribution(coords, num_bins, bound, index):
        count = np.zeros(num_bins)
        height_range = 2 * bound
        for i in range(len(coords)):
            height = coords[i, index]
            if height > bound * 1.5 or height < -0.5 * bound:
                continue
            
            height_index = ((height + 0.5 * bound)/height_range) * num_bins
            height_index = np.trunc(height_index).astype(int) - 1
            count[height_index] += 1
        
        return count
    
    def get_distribution_density(*histograms, density = 3):
        max_val = 0
        to_return = ()
        for val in histograms:
            max_val = val.max() if max_val < val.max() else max_val
        
        for val in histograms:
            val = (val / max_val) * density
            to_return += (val,)
        
        return to_return
        
    index, = np.where(abs(plane[0:3]) == abs(plane[0:3]).max())
    
    snap1 = unwrap_pbc_com(snap1, bound)
    snap2 = unwrap_pbc_com(snap2, bound)
    snap3 = unwrap_pbc_com(snap3, bound)

    snap1_trans = rotate_stuff(plane, snap1, rotate_all = False)
    snap2_trans = rotate_stuff(plane, snap2, rotate_all = False)
    snap3_trans = rotate_stuff(plane, snap3, rotate_all = False)
    
    # Filter out irrelevant particles
    for i in range(3):
        if i == index:
            continue
        snap1_trans = snap1_trans[snap1_trans[:, i] < bound]
        snap1_trans = snap1_trans[snap1_trans[:, i] > 0]
        snap2_trans = snap2_trans[snap2_trans[:, i] < bound]
        snap2_trans = snap2_trans[snap2_trans[:, i] > 0]
        snap3_trans = snap3_trans[snap3_trans[:, i] < bound]
        snap3_trans = snap3_trans[snap3_trans[:, i] > 0]
        
    axis = np.linspace(-0.5 * bound, 1.5 * bound, num_bins)
    count1 = get_distribution(snap1_trans, num_bins, bound, index)
    count2 = get_distribution(snap2_trans, num_bins, bound, index)
    count3 = get_distribution(snap3_trans, num_bins, bound, index)
    
    bin_volume = bound * bound * (2 * bound/num_bins)
    c1 = count1/bin_volume
    c2 = count2/bin_volume
    c3 = count3/bin_volume
    
    return c1, c2, c3, axis

def get_bilayer_histogram_norotate(snap1, snap2, snap3, num_bins, bound, normal = 2):
    def get_distribution(coords, num_bins, bound, index):
        count = np.zeros(num_bins)
        height_range = 2 * bound
        for i in range(len(coords)):
            height = coords[i, index]
            if height > bound * 1.5 or height < -0.5 * bound:
                continue
            
            height_index = ((height + 0.5 * bound)/height_range) * num_bins
            height_index = np.trunc(height_index).astype(int) - 1
            count[height_index] += 1
        
        return count
        
    index = normal
    
    snap1 = unwrap_pbc_com(snap1, bound)
    snap2 = unwrap_pbc_com(snap2, bound)
    snap3 = unwrap_pbc_com(snap3, bound)
    
    # Filter out irrelevant particles
    for i in range(3):
        if i == index:
            continue
        snap1 = snap1[snap1[:, i] < bound]
        snap1 = snap1[snap1[:, i] > 0]
        snap2 = snap2[snap2[:, i] < bound]
        snap2 = snap2[snap2[:, i] > 0]
        snap3 = snap3[snap3[:, i] < bound]
        snap3 = snap3[snap3[:, i] > 0]
        
    axis = np.linspace(-0.5 * bound, 1.5 * bound, num_bins)
    count1 = get_distribution(snap1, num_bins, bound, index)
    count2 = get_distribution(snap2, num_bins, bound, index)
    count3 = get_distribution(snap3, num_bins, bound, index)
    
    bin_volume = bound * bound * (2 * bound/num_bins)
    c1 = count1/bin_volume
    c2 = count2/bin_volume
    c3 = count3/bin_volume
    
    return c1, c2, c3, axis

def get_p_profile(coords, num_bins, plane, bound):
    
    index, = np.where(abs(plane[0:3]) == abs(plane[0:3]).max())
    non_index = np.array([0, 1, 2])
    non_index = np.delete(non_index, index)
    
    coords = unwrap_pbc_com(coords, bound)
    coords = rotate_stuff(plane, coords)
    
    for i in range(3):
        if i == index:
            continue
        coords = coords[coords[:, i] < bound]
        coords = coords[coords[:, i] > 0]
        
    axis = np.linspace(-0.5 * bound, 1.5 * bound, num_bins)
    
    p_profile = np.zeros(num_bins)
    count = np.zeros(num_bins)
    height_range = 2 * bound
    
    for i in range(len(coords)):
        height = coords[i, index]
        if height > bound * 1.5 or height < -0.5 * bound:
            continue
        
        height_index = ((height + 0.5 * bound)/height_range) * num_bins
        height_index = np.trunc(height_index).astype(int) - 1
        
        p_profile[height_index] += coords[i, 3 + index] - \
            (coords[i, 3 + non_index[0]] + coords[i, 3 + non_index[1]]) * 0.5
        count[height_index] += 1
    
    p_profile = p_profile/count
    
    return p_profile, axis

def get_p_profile_norotate(coords, num_bins, bound, normal = 2):
    index = normal
    non_index = np.array([0, 1, 2])
    non_index = np.delete(non_index, index)
    
    coords = unwrap_pbc_com(coords, bound)
    
    for i in range(3):
        if i == index:
            continue
        coords = coords[coords[:, i] < bound]
        coords = coords[coords[:, i] > 0]
        
    axis = np.linspace(-0.5 * bound, 1.5 * bound, num_bins)
    
    p_profile = np.zeros(num_bins)
    count = np.zeros(num_bins)
    height_range = 2 * bound
    
    for i in range(len(coords)):
        height = coords[i, index]
        if height > bound * 1.5 or height < -0.5 * bound:
            continue
        
        height_index = ((height + 0.5 * bound)/height_range) * num_bins
        height_index = np.trunc(height_index).astype(int) - 1
        
        p_profile[height_index] += coords[i, 3 + index] - \
            (coords[i, 3 + non_index[0]] + coords[i, 3 + non_index[1]]) * 0.5
        count[height_index] += 1
    
    p_profile = p_profile/count
    
    return p_profile, axis
