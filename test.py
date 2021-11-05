import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=5, rho=0.5, sigma=1, T=0.1)

# Calculate force for each particle instances
simulation.update_force(simulation.Particles[2], epsilon=1, alpha=1, r_cut=2.5)

fig, ax = plt.subplots()
sc = ax.scatter(
    *np.array([i.position for i in simulation.Particles]).T, s=10)
plt.xlim(0, simulation.size)
plt.ylim(0, simulation.size)


def animate(i):
    simulation.move(0.001, m=1)
    sc.set_offsets(np.array([i.position for i in simulation.Particles]))
    return sc


ani = animation.FuncAnimation(fig, animate, interval=200)
plt.show()
