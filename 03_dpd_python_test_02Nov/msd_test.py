#%% Imports and setup
from lammps import PyLammps
from pathlib import Path
import os

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
L = PyLammps()

#%% Set some variables
box_length = 10
num_density = 3 # num_density being the number of particles in one unit

dt = 0.04
pressure = 1
temp = 1

#%% Simulation Settings
L.timestep(dt)
L.comm_modify("vel", "yes")
L.units("lj")
L.atom_style("dpd")
L.pair_style("dpd/fdt", 1, 1, 22445)
L.neighbor(5.0, "bin")
L.neigh_modify("every", 2)

#%% Create simulation box and populate with particles
L.lattice("sc", num_density)
L.region("box", "block", 0, box_length, 0, box_length, 0, box_length,
         "units", "box")
L.create_box(1, "box")
L.create_atoms(1, "box")

#%% Set mass, pair coefficients, velocities
L.mass(1, 1.0)
# Velocity in direction of x-axis, ramping in z- direction from low to high
L.velocity("all", "create", 1, 696969)
# a_ij = 100, sigma = 3, r_c = 1??
L.pair_coeff(1, 1, 100, 3, 1)

#%% Set fixes
# Fixes isothermal-isobaric with a fixed z-height
L.fix("1", "all", "npt", "temp", 1, 1, dt * 100, "z ", str(pressure), 
      str(pressure), str(dt * 1000))

# Fixes lambda = 0.65 for Velocity-Verlet integration
# Needs atom style mdpd???
# L.fix("2", "all", "mvv/dpd", 0.65)

#%% Set thermo log output
L.thermo_style("custom", "step", "temp", "press", "etotal", "density", "vol")
L.thermo(1000)

#%% Equilibrium run
L.run(50000)

#%% Change to nvt
L.unfix("1")
L.fix("2", "all", "nvt", "temp", temp, temp, dt * 100)

#%% Data gathering
L.reset_timestep(0)

L.variable("density equal density")
L.fix("dense", "all", "ave/time", 20, 1, 20, "v_density", 
      "file", "density.txt")

L.compute("msd", "all", "msd", "com", "yes")
L.variable("msd_scalar equal c_msd[4]")
L.fix("msd_vec", "all", "ave/time", 20, 1, 20, "v_msd_scalar", 
      "file", "msd_vec.txt")

#%%
L.run(500000)
