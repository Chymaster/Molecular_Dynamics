import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=10, rho=0.3, sigma=0.5, T=5, dt=0.001)

# Calculate force for each particle instances
simulation.update_force(
    simulation.Particles[2], epsilon=1, alpha=1, r_cut=2.5)

fig, ax = plt.subplots()
sc = ax.scatter(
    *np.array([i.position for i in simulation.Particles]).T, s=50)
plt.xlim(0, simulation.size)
plt.ylim(0, simulation.size)


def animate(i):
    simulation.move(m=1, log_file="log_test.txt", f_log=0.01)
    sc.set_offsets(np.array([i.position for i in simulation.Particles]))
    return sc


ani = animation.FuncAnimation(fig, animate, interval=20)
plt.show()
