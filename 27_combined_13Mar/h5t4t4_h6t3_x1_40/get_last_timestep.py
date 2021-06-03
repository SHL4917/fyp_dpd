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
dump = "bilayer_conc_014/dump.allparticles"
input_file = "bilayer_conc_014/data.txt"

dump_micelle = "micelle_conc_010/dump.all"
input_file_micelle = "micelle_conc_010/data.txt"
poly_num_micelle = np.load("micelle_conc_010/poly_num.npy")
bound = 36

#%%
print("Getting number of particles from bilayer dump...")
i = 0
with open(dump) as infile:
    for line in infile:
        if i == 3:
            tot_particles = np.fromstring(line, dtype = "int", sep = " ")[0]
        elif i == 8:
            break
        i = i + 1

# First column is id, second column is type, last three are the xyz coords
coords = np.zeros((tot_particles, 5))

#%% Get latest timestep
print("Getting bilayer coordinates...")
timestep_line = 0
i = 0
with open(dump) as infile:
    for line in infile:
        i += 1
        if line == "ITEM: TIMESTEP\n":
            timestep_line = i + 1
            continue
        elif i == timestep_line:
            timestep = int(line)

i = 0
j = 0
with open(dump) as infile:
    for line in infile:
        i += 1
        if i > timestep_line + 7:
            data = np.fromstring(line, sep = " ")
            coords[j, :] = data[0:5]
            j += 1

coords[:, 2:5] *= bound

coords1 = coords.copy()
coords1[:, 0] = coords1[:, 0] + tot_particles
coords1[:, 4] = coords1[:, 4] + bound + bound
coords = np.concatenate((coords, coords1), axis = 0)

#%%
print("Getting bilayer angle and bond information...")
i = 0
with open(input_file) as infile:
    for line in infile:
        i += 1
        if i == 4:
            bond = int(line.split(" ")[0])
        elif i == 5:
            angle = int(line.split(" ")[0])
        elif line == "Bonds\n":
            bond_line = i
        elif line == "Angles\n":
            angle_line = i
            
bond_coords = np.zeros((bond, 4))
angle_coords = np.zeros((angle, 5))
i = 0
j = 0
k = 0

with open(input_file) as infile:
    for line in infile:
        i += 1
        if i > bond_line + 1 and i < angle_line - 1:
            data = np.fromstring(line, sep = " ")
            bond_coords[j, :] = data
            j += 1
        if i > angle_line + 1:
            data = np.fromstring(line, sep = " ")
            angle_coords[k, :] = data
            k += 1

bond_coords1 = bond_coords.copy()
bond_coords1[:, 0] += bond
bond_coords[:, 2:4] += tot_particles

bond_coords = np.concatenate((bond_coords, bond_coords1), axis = 0)

angle_coords1 = angle_coords.copy()
angle_coords1[:, 0] += angle
angle_coords1[:, 2:5] += tot_particles

angle_coords = np.concatenate((angle_coords, angle_coords1), axis = 0)

#%%
print("Getting number of particles from micelle dump...")
i = 0
with open(dump_micelle) as infile:
    for line in infile:
        if i == 3:
            tot_particles_micelle = np.fromstring(line, dtype = "int", 
                                                  sep = " ")[0]
        elif i == 8:
            break
        i = i + 1

# First column is id, second column is type, last three are the xyz coords
coords_micelle = np.zeros((tot_particles_micelle, 5))

#%% Get latest timestep
print("Getting micelle coordinates...")
timestep_line = 0
i = 0
with open(dump_micelle) as infile:
    for line in infile:
        i += 1
        if line == "ITEM: TIMESTEP\n":
            timestep_line = i + 1
            continue
        elif i == timestep_line:
            timestep = int(line)

i = 0
j = 0
with open(dump_micelle) as infile:
    for line in infile:
        i += 1
        if i > timestep_line + 7:
            data = np.fromstring(line, sep = " ")
            coords_micelle[j, :] = data[0:5]
            j += 1

coords_micelle[:, 2:5] *= bound

#%% Rectifying the PBC conditions for micelles
for i in range(len(poly_num_micelle)):
    for j in range(np.shape(poly_num_micelle)[1]):
        if j == 0:
            bead = poly_num_micelle[i, j]
            id_num, = np.where(coords_micelle[:, 0] == bead)[0]
            xyz_head = coords_micelle[id_num, 2:5]
            continue
        
        bead = poly_num_micelle[i, j]
        id_num, = np.where(coords_micelle[:, 0] == bead)[0]
        for k in range(3):
            diff = xyz_head[k] - coords_micelle[id_num, k + 2]
            if abs(diff) > bound * 0.7 and np.sign(diff) > 0:
                coords_micelle[id_num, k + 2] += bound
            elif abs(diff) > bound * 0.7 and np.sign(diff) < 0:
                coords_micelle[id_num, k + 2] -= bound
            

#%%
coords_micelle[:, 4] += bound
coords_micelle[:, 0] += np.shape(coords)[0]
coords_micelle[coords_micelle[:, 1] > 1, 1] += 2

#%%
print("Getting micelle angle and bond information...")
i = 0
with open(input_file_micelle) as infile:
    for line in infile:
        i += 1
        if i == 4:
            bond = int(line.split(" ")[0])
        elif i == 5:
            angle = int(line.split(" ")[0])
        elif line == "Bonds\n":
            bond_line = i
            
bond_coords_micelle = np.zeros((bond, 4))
i = 0
j = 0

with open(input_file_micelle) as infile:
    for line in infile:
        i += 1
        if i > bond_line + 1:
            data = np.fromstring(line, sep = " ")
            bond_coords_micelle[j, :] = data
            j += 1

bond_coords_micelle[:, 0] += np.shape(bond_coords)[0]
bond_coords_micelle[:, 2:4] += np.shape(coords)[0]

#%%
all_coords = np.concatenate((coords, coords_micelle), axis = 0)
all_bond_coords = np.concatenate((bond_coords, bond_coords_micelle), axis = 0)
all_angle_coords = angle_coords

#%% Write to a datafile for now
print("Writing to file...")
new_input_file = "combined_input.txt"

mass_array = np.array([1, 1, 1, 1, 1])

particle_count = np.shape(all_coords)[0]
bonds = np.shape(all_bond_coords)[0]
angles = np.shape(all_angle_coords)[0]
dihedrals = 0
impropers = 0

atom_types = 5
bond_types = all_bond_coords[:, 1].max()
angle_types = all_angle_coords[:, 1].max()
dihedral_types = 0
improper_types = 0

with open(new_input_file, mode = "w") as line:
    line.write("LAMMPS " + "\n \n")
            
    line.write(str(particle_count) + " atoms\n")
    line.write(str(bonds) + " bonds\n")
    line.write(str(angles) + " angles\n")
    line.write(str(dihedrals) + " dihedrals\n")
    line.write(str(impropers) + " impropers\n \n")
    
    line.write(f"{atom_types:.0f} atom types\n")
    line.write(f"{bond_types:.0f} bond types\n")
    line.write(f"{angle_types:.0f} angle types\n")
    line.write(f"{dihedral_types:.0f} dihedral types\n")
    line.write(f"{improper_types:.0f} improper types\n \n")
    
    line.write("0 " + str(bound) + " xlo xhi\n")
    line.write("0 " + str(bound) + " ylo yhi\n")
    line.write("0 " + str(bound * 3) + " zlo zhi\n \n")
    
    line.write("Masses\n \n")
    for index, val in enumerate(mass_array):
        line.write(f"{index + 1} {val}\n")
        
    line.write("\n")
    
    line.write("Atoms\n \n")
    
    for i in range(len(all_coords)):
        data = all_coords[i, :]
        line.write(f"{data[0]:.0f} 0 {data[1]:.0f} "
                   f"{data[2]:.3f} "
                   f"{data[3]:.3f} "
                   f"{data[4]:.3f}\n")
        
    line.write("\n")
    line.write("Bonds\n \n")
    
    for i in range(len(all_bond_coords)):
        data = all_bond_coords[i, :]
        line.write(f"{data[0]:.0f} {data[1]:.0f} "
                   f"{data[2]:.0f} "
                   f"{data[3]:.0f}\n")
    
    line.write("\n")
    line.write("Angles\n \n")
    
    for i in range(len(all_angle_coords)):
        data = all_angle_coords[i, :]
        line.write(f"{data[0]:.0f} {data[1]:.0f} "
                   f"{data[2]:.0f} "
                   f"{data[3]:.0f} "
                   f"{data[4]:.0f}\n")
        
        