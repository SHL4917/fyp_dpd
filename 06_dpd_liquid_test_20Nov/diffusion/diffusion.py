#%% Imports and setup
from lammps import PyLammps
from pathlib import Path
import os

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
L = PyLammps()
base_dir = os.getcwd()

#%% Set some variables
box_length = 10
num_density = 2.7 # num_density being the number of particles in one unit

dt = 0.04
pressure = 24
temp = 1

a_ij= [50]
sigma = 3
r_c = 1
mass = 1.0

#%% Simulation Settings
L.timestep(dt)
L.comm_modify("vel", "yes")
L.units("lj")
L.atom_style("dpd")
L.pair_style("dpd/fdt", 1, 1, 22445)
L.neighbor(5.0, "bin")
L.neigh_modify("every", 2, "one", 3375)

#%% Create simulation box and populate with particles
L.lattice("sc", num_density)
L.region("box", "block", 0, box_length, 0, box_length, 0, box_length,
         "units", "box")
L.create_box(1, "box")
L.create_atoms(1, "box")

#%% Run simulation for different values of a_ij
for i in range(len(a_ij)):
    save_folder = "/a_ij_" + str(a_ij[i]) + "/"

    if not os.path.exists(os.getcwd() + save_folder):
        os.makedirs(os.getcwd() + save_folder)

    os.chdir(os.getcwd() + save_folder)
    
    # code, save
    L.mass(1, mass)
    L.velocity("all", "create", 1, 696969)
    L.pair_coeff(1, 1, a_ij[i], sigma, r_c)
    
    # Fixes isothermal-isobaric with a fixed z-height
    L.fix("1", "all", "npt", "temp", 1, 1, dt * 100, "z ", 
          str(pressure), str(pressure), str(dt * 1000))
    
    L.thermo_style("custom", "step", "temp", "press", "etotal", "density", 
                   "vol")
    L.thermo(1000)
    
    L.reset_timestep(0)    
    L.run(50000)
    
    L.unfix("1")
    L.fix("2", "all", "nvt", "temp", temp, temp, dt * 100)
    
    L.reset_timestep(0)
    
    L.variable("density equal density")
    L.fix("dense", "all", "ave/time", 200, 1, 200, "v_density", 
          "file", "density.txt")
    
    L.compute("msd", "all", "msd", "com", "yes")
    L.variable("msd_scalar equal c_msd[4]")
    L.fix("msd_vec", "all", "ave/time", 200, 1, 200, 
          "v_msd_scalar", "file", "msd_vec.txt")
    
    L.run(200000)
    os.chdir(base_dir)
    
    L.unfix("2")
    L.unfix("dense")
    L.uncompute("msd")
    L.unfix("msd_vec")
    
    
#%%