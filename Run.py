import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=20, rho=0.1, sigma=1, T=0.1)

# Calculate force for a particle instance
#simulation.update_force(
#    simulation.Particles[15], epsilon=1, alpha=0, r_cut=2.5)


# Map the particles in the grid
#simulation.map()
