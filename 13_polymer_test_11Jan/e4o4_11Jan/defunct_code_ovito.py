from ovito.data import *
from ovito.pipeline import *

#%%
coords = com[:, :, 200]

# Create the data collection containing a Particles object:
particles = Particles()
data = DataCollection()
data.objects.append(particles)

# Create the particle position property:
pos_prop = particles.create_property('Position', data = coords)

# Create the particle type property and insert two atom types:
type_prop = particles.create_property('Particle Type')
type_prop.types.append(ParticleType(id = 1, name = 'Cu', color = (0.0,1.0,0.0)))
type_prop.types.append(ParticleType(id = 2, name = 'Ni', color = (0.0,0.5,1.0)))
type_prop[0] = 1  # First atom is Cu
type_prop[1] = 2  # Second atom is Ni
type_prop[2] = 2  # Third atom is Ni

# Create the simulation box:
cell = SimulationCell(pbc = (True, True, True))
# First three columns represent cell vectors, last column position of origin
cell[...] = [[bounds,0,0,0],
             [0,bounds,0,0],
             [0,0,bounds,0]]
cell.vis.line_width = 0.1
data.objects.append(cell)

#%%
# Create a pipeline, set source and insert it into the scene:
pipeline = Pipeline(source = StaticSource(data = data))
pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=2.8, 
                                                  sort_by_size=True,
                                                  compute_com=True,
                                                  unwrap_particles = True))

#%%
# Get computed data
cluster_data = pipeline.compute()