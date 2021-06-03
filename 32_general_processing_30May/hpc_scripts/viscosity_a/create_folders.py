import numpy as np
import sys
import os
from pathlib import Path
import re
from shutil import copy2

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
base = os.getcwd()

sys.path.append("../../../")
from fyp_functions import *

#%%
aij_values = [15, 20, 25, 30, 35, 40, 50]
folder_name = [f"aij_{x:.0f}" for x in aij_values]

#%%
input_script = []

for i, aij in enumerate(aij_values):
    input_script += [f"""
units		lj
timestep	0.04
comm_modify	vel yes
atom_style	atomic
pair_style	dpd 1 1 22445
neighbor	2 bin
neigh_modify	every 1 delay 0 check yes one 3000

region		box block 0 10 0 10 0 10 units box
create_box	1 box
create_atoms	1 random 3000 696969 box

mass		1 1

pair_coeff	1 1 {aij} 4.5 1

fix		nve all nve

run		100000

unfix		nve

change_box	all triclinic
fix		1 all nvt temp 1 1 4
fix		2 all deform 1 xz erate 0.15 remap v

run		50000
reset_timestep	0

thermo_style	custom step temp press density vol
thermo		1000

compute		temp all temp
variable	avetemp equal c_temp
fix		avtemp all ave/time 50 200 10000 v_avetemp file avtemp.txt

compute		pressure all pressure temp
variable	avepressure equal c_pressure
fix		avpressure all ave/time 50 200 10000 v_avepressure file avpressure.txt

compute		chunk_vx all chunk/atom bin/1d z center 0.1 units reduced
fix		velo_profile_vx all ave/chunk 50 200 10000 chunk_vx vx file velo_profile_vx.txt

fix		rhostar_profile all ave/chunk 50 200 10000 chunk_vx density/number file rhostar_profile.txt
fix		tstar_profile all ave/chunk 50 200 10000 chunk_vx temp file tstar_profile.txt

variable	pxz equal pxz
fix		press_xz all ave/time 50 200 10000 v_pxz file press_xz.txt

dump mydump all atom 1000 dump.atom

run		800000
                     """]
                     
#%%
for i, folder in enumerate(folder_name):
    directory = f"{base}\{folder}"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    copy2("viscosity.pbs", directory)
    
    print(input_script[i], file = open(f"{directory}\in.viscosity", "w"))
    