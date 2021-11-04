import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=10, rho=0.5, sigma=1, T=0.1)

# Calculate force for each particle instances
simulation.update_force(simulation.Particles[2], epsilon=1, alpha=0, r_cut=2.5)

# Map the particles in the grid
"""
plt.figure(1)
plt.scatter(
    *np.array([i.position for i in simulation.Particles]).T, color="blue")
"""
for i in range(10):
    simulation.move(0.1, m=1)
    print(simulation.Particles[14].position)
