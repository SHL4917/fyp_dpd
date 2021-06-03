import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../")
from fyp_functions import *
   
#%%
box_size = 30
polymer_head = np.array([2, 2, 2, 2, 2])
polymer_t1 = np.array([3, 3, 3, 3])
polymer_t2 = np.array([3, 3, 3, 3])
filename = "data.txt"
mass_array = [1, 1, 1]
rho_star = 3
description = "nothing"
dense = 0.7


#%%
    
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

gridsize = np.ceil((particle_count ** (1/3))).astype(int)
spacing = box_size / (gridsize - 1)

bilayer_top_chains = np.floor(((gridsize * gridsize)/2) * dense).astype(int)
num_chains = np.floor(bilayer_top_chains * 2 * 1.1).astype(int)
num_solvent = particle_count - num_chains * polymer_length
polymer_conc = (num_chains * polymer_length)/particle_count

poly_index = np.zeros([num_chains, polymer_length, 3])

head_thickness = np.ceil(len(polymer_head)/2).astype(int)
tail_thickness = max(len(polymer_t1), len(polymer_t2))

angles = num_chains * get_angle(poly_num)
bonds = num_chains * (polymer_length - 1)

height = (gridsize/4).astype(int)
even = True if len(polymer_head) % 2 == 0 else False

# Create matrix to populate with particles
matrix = np.zeros((gridsize, gridsize, gridsize))

#%%
cc = 0
for i in range(gridsize):
    for j in range(0, gridsize, 2):
        if np.random.uniform() > dense:
            continue
        
        pc = 0
        for k in range(head_thickness):
            z = k + height
            if k == 0 and even:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
                poly_index[cc, pc, :] = np.array([i, j + 1, z])
                matrix[i, j + 1, z] = 2
                pc += 1
            elif k == 0:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
            else:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
                poly_index[cc, pc, :] = np.array([i, j + 1, z])
                matrix[i, j + 1, z] = 2
                pc += 1
                
        for k in range(tail_thickness):
            z = k + height + head_thickness
            poly_index[cc, pc, :] = np.array([i, j, z])
            matrix[i, j, z] = 3
            pc += 1
        
        
        for k in range(tail_thickness):
            z = k + height + head_thickness
            poly_index[cc, pc, :] = np.array([i, j + 1, z])
            matrix[i, j + 1, z] = 3
            pc += 1
            
        cc += 1
        pc = 0
        
#%%
for i in range(gridsize):
    for j in range(0, gridsize, 2):
        if np.random.uniform() > dense:
            continue
        
        pc = 0
        for k in range(head_thickness):
            z = height + 2 * (head_thickness + tail_thickness) - k - 1
            if k == 0 and even:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
                poly_index[cc, pc, :] = np.array([i, j + 1, z])
                matrix[i, j + 1, z] = 2
                pc += 1
            elif k == 0:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
            else:
                poly_index[cc, pc, :] = np.array([i, j, z])
                matrix[i, j, z] = 2
                pc += 1
                
                poly_index[cc, pc, :] = np.array([i, j + 1, z])
                matrix[i, j + 1, z] = 2
                pc += 1
        
        for k in range(tail_thickness):
            z = height + 2 * (tail_thickness)  + head_thickness - k - 1
            poly_index[cc, pc, :] = np.array([i, j, z])
            matrix[i, j, z] = 3
            pc += 1
        
        for k in range(tail_thickness):
            z = height + 2 * (tail_thickness)  + head_thickness - k - 1
            poly_index[cc, pc, :] = np.array([i, j + 1, z])
            matrix[i, j + 1, z] = 3
            pc += 1
        
        cc += 1
        pc = 0

#%%
while cc < np.shape(poly_index)[0]:
    x = np.random.randint(0, gridsize - polymer_length)
    y = np.random.randint(0, gridsize - 1)
    z = np.random.randint(0, gridsize - 1)
    if got_space_poly(x, y, z, matrix, poly_num):
        pc = 0
        for i in range(poly_num[0]):
            poly_index[cc, pc, :] = np.array([x + i, y, z])
            matrix[x + i, y, z] = 2
            pc += 1
            
        for i in range(poly_num[1]):
            poly_index[cc, pc, :] = np.array([x + i + poly_num[0], y, z])
            matrix[x + i + poly_num[0], y, z] = 3
            pc += 1
            
        for i in range(poly_num[2]):
            poly_index[cc, pc, :] = np.array([x + i + poly_num[0], y + 1, z])
            matrix[x + i + poly_num[0], y + 1, z] = 3
            pc += 1
        
        cc += 1
#%%

# Populate remaining with solvent beads
i = 0
for x in range(gridsize):
    for y in range(gridsize):
        for z in range(gridsize):
            if (matrix[x][y][z] == 0) and (i < num_solvent):
                matrix[x][y][z] = 1
                i = i + 1
                
#%%
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
