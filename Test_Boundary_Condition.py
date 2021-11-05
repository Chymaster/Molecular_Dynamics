import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=10, rho=0.7, sigma=1, T=0.1)

# Calculate force for each particle instances
"""for i in range(len(simulation.Particles)):
    simulation.update_force(
        simulation.Particles[i], epsilon=1, alpha=0, r_cut=2.5)"""

simulation.update_force(
    simulation.Particles[2], epsilon=1, alpha=0, r_cut=2.5)
print(simulation.Particles[2].position, [
      i.position for i in simulation.Particles[2].neighbour_list])
simulation.map(particles=[simulation.Particles[2]],
               neighbours=simulation.Particles[2].neighbour_list)

# Map the particles in the grid
#simulation.map()
