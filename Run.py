import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=10, rho=0.1, sigma=1, T=0.1)

# Calculate force for each particle instances
simulation.update_force(simulation.Particles[2], epsilon=1, alpha=1, r_cut=2.5)

# Map the particles in the grid

plt.figure(1)
plt.scatter(
    *np.array([i.position for i in simulation.Particles]).T, color="blue")

for i in range(1000):
    simulation.move(0.001, m=1)
    print("velocity is", max(np.absolute([max(np.absolute(i.velocity))
                                          for i in simulation.Particles])), "distance between two particle is", min(i.distance for i in simulation.Particles))
