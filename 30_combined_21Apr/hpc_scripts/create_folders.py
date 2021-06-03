import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import re
from shutil import copy2

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

#%%
x1_array = np.linspace(14, 22, 10)
x1_array = list(x1_array)
x1_array = [f"{x:.2f}" for x in x1_array]

folder_name = [f"h5t4t4_h6t3_x1_{float(x) * 100:.0f}" for x in x1_array]

#%%
combined_first = []
combined_next = []

for i, aij in enumerate(x1_array):
    combined_first += [f"""
units			lj
comm_modify		cutoff 8 vel yes
atom_style		angle
pair_style		dpd 1 1 22445
neighbor		5 bin
neigh_modify		every 1 delay 0 check yes one 3000

timestep		0.04

read_data		combined_input.txt

bond_style		harmonic
bond_coeff		1 100 0.7

angle_style		harmonic
angle_coeff		1 6 180
angle_coeff		2 3 120
angle_coeff		3 4.5 120

pair_style		dpd 1 1 696969

velocity		all create 1 6969
pair_coeff		1 1 25 4.5 1
pair_coeff		1 2 25 4.5 1
pair_coeff		1 3 100 4.5 1
pair_coeff		2 2 25 4.5 1
pair_coeff		2 3 100 4.5 1
pair_coeff		3 3 25 4.5 1

pair_coeff		4 4 25 4.5 1
pair_coeff		5 5 25 4.5 1
pair_coeff		1 4 25 4.5 1
pair_coeff		1 5 100 4.5 1
pair_coeff		2 5 100 4.5 1
pair_coeff		3 4 100 4.5 1
pair_coeff		3 5 25	4.5 1
pair_coeff		4 5 100 4.5 1

pair_coeff		2 4 {aij} 4.5 1

fix			1 all nvt temp 1 1 4

run			20000

reset_timestep		0
thermo 			20000
thermo_style		custom step temp press density

compute		ke all ke/atom
compute		pe all pe/atom

dump			allparticles all custom 20000 dump.allparticles id type xs ys zs c_ke c_pe

run			200000
write_restart	restart.*

print	"$(step)" file last_timesteps.txt
    """]
    
    combined_next += [f"""
read_restart	restart.*

comm_modify		cutoff 8 vel yes
neighbor		5 bin
neigh_modify		every 1 delay 0 check yes one 3000

timestep		0.04

bond_style		harmonic
bond_coeff		1 100 0.7

angle_style		harmonic
angle_coeff		1 6 180
angle_coeff		2 3 120
angle_coeff		3 4.5 120

pair_coeff		1 1 25 4.5 1
pair_coeff		1 2 25 4.5 1
pair_coeff		1 3 100 4.5 1
pair_coeff		2 2 25 4.5 1
pair_coeff		2 3 100 4.5 1
pair_coeff		3 3 25 4.5 1

pair_coeff		4 4 25 4.5 1
pair_coeff		5 5 25 4.5 1
pair_coeff		1 4 25 4.5 1
pair_coeff		1 5 100 4.5 1
pair_coeff		2 5 100 4.5 1
pair_coeff		3 4 100 4.5 1
pair_coeff		3 5 25	4.5 1
pair_coeff		4 5 100 4.5 1

pair_coeff		2 4 {aij} 4.5 1

fix			1 all nvt temp 1 1 4

thermo 			20000
thermo_style		custom step temp press density

compute		ke all ke/atom
compute		pe all pe/atom

dump			allparticles all custom 20000 dump.allparticles id type xs ys zs c_ke c_pe
dump_modify		allparticles append yes

run			200000
write_restart	restart.*

print	"$(step)" append last_timesteps.txt
                      """]


#%%
for i, folder in enumerate(folder_name):
    directory = f"{base}\{folder}"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    copy2("combined_first.pbs", directory)
    copy2("combined_next.pbs", directory)
    copy2("combined_input.txt", directory)
    
    print(combined_first[i], file = open(f"{directory}\in.combined_first", 
                                         "w"))
    print(combined_next[i], file = open(f"{directory}\in.combined_subsequent", 
                                        "w"))
        