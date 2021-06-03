#%% Imports and setup
from lammps import PyLammps
import os
from pathlib import Path

__dir__ = Path(globals().get("__file__", "./_")).absolute().parent
os.chdir(__dir__)
L = PyLammps()

#%% Set some variables
box_length = 10
num_density = 2.7 # num_density being the number of particles in one unit

v_max = 1.5
s_rate = v_max/box_length

dt = 0.04
pressure = 24

a_ij= 50
sigma = 3
r_c = 1
mass = 1.0

save_folder = "/a_ij_" + str(a_ij) + "/"

if not os.path.exists(os.getcwd() + save_folder):
    os.makedirs(os.getcwd() + save_folder)

os.chdir(os.getcwd() + save_folder)

#%% Simulation Settings
L.timestep(dt)
L.comm_modify("vel", "yes")
L.units("lj")
L.atom_style("dpd")
L.pair_style("dpd/fdt", 1, 1, 22445)
L.neighbor(5.0, "bin")
L.neigh_modify("every", 2)

#%% Create simulation box and populate with particles
L.region("box", "block", 0, box_length, 0, box_length, 0, box_length,
         "units", "box")
L.create_box(1, "box")
L.create_atoms(1, "random", 3000, 696969, "box")

#%% Set mass, pair coefficients, velocities
L.mass(1, mass)
# Velocity in direction of x-axis, ramping in z- direction from low to high
L.velocity("all", "ramp", "vx", 0.0, v_max, "z", 0, box_length)
# a_ij = 100, sigma = 3, r_c = 1??
L.pair_coeff(1, 1, a_ij, sigma, r_c)

#%% Set nve fix
L.fix("nve", "all", "nve")

#%% Equilibrium run
L.run(50000)

#%% Set fixes
L.unfix("nve")
# Must put change_box command before the fixes!
L.change_box("all", "triclinic")
# Fixes isothermal-isobaric with a fixed z-height
L.fix("1", "all", "npt", "temp", 1, 1, dt * 100, "z ", pressure, 
      pressure, dt * 1000)
# Sets a deforming box equivalent to Lee-Edwards BC
L.fix("2", "all", "deform", 1, "xz", "erate", s_rate, "remap", "v")

#%% Set thermo log output
L.thermo_style("custom", "step", "temp", "press", "etotal", "density", "vol",
               "nbuild", "ndanger")
L.thermo(1000)

#%% Equilibrium run
L.run(50000)

#%% Set computes
# Get average temperature
L.compute("temp", "all", "temp")
L.variable("avetemp equal c_temp")
L.fix("avtemp", "all", "ave/time", 10, 200, 2000, "v_avetemp", 
      "file", "avtemp.txt")

# Get average pressure
L.compute("pressure", "all", "pressure", "temp")
L.variable("avepressure equal c_pressure")
L.fix("avpressure", "all", "ave/time", 10, 200, 2000, "v_avepressure", 
      "file", "avpressure.txt")

# Get velocity profiles in x-direction
L.compute("chunk_vx", "all", "chunk/atom", "bin/1d", "z", "center", "0.1", 
          "units", "reduced")
L.fix("velo_profile_vx", "all", "ave/chunk", 10, 200, 2000, "chunk_vx", "vx", 
      "file", "velo_profile_vx.txt")

# Get velocity profiles in z-direction, should average to zero
L.compute("chunk_vz", "all", "chunk/atom", "bin/1d", "x", "center", "0.1", 
          "units", "reduced")
L.fix("velo_profile_vz", "all", "ave/chunk", 10, 200, 2000, "chunk_vz", "vz", 
      "file", "velo_profile_vz.txt")

# Get bulk pressure values in six unique directions
L.variable("pxz equal pxz")
L.fix("press_xz", "all", "ave/time", 50, 200, 10000, "v_pxz", 
      "file", "press_xz.txt")

# Get density values
L.variable("dense equal density")
L.fix("density", "all", "ave/time", 50, 200, 10000, "v_dense", 
      "file", "density.txt")

# Dump particle positions
L.dump(1, "all", "atom", 100, "dump.viscosity" + str(a_ij))

#%% Actual run
L.reset_timestep(0)
L.run(500000)

#%%
print("Finished!")
