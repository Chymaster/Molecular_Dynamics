import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=20, rho=0.1, sigma=0.5, T=5, dt=0.001)


# Making N moves
N = 10000

temperature = np.zeros(N)
kinetic_energy = np.zeros(N)
potential_energy = np.zeros(N)
total_energy = np.zeros(N)

for i in range(N):
    temperature[i] = simulation.T
    kinetic_energy[i] = simulation.kinetic_E
    potential_energy[i] = simulation.potential_E
    total_energy[i] = simulation.E
    simulation.move(m=1, alpha=0, log_file="", f_log=0.0001)
    if simulation.kinetic_E > 3e5:
        print("T", simulation.T, "kinetic_energy", simulation.kinetic_E)
        print("potential_energy", simulation.potential_E, "total_energy", simulation.E)
        print(max([np.linalg.norm(i.velocity) for i in simulation.Particles]))

# Ploting Temperature, kinetic energy, potential energy and Total energy
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.title.set_text('Temperature')
ax1.plot(temperature)

ax2 = fig.add_subplot(222)
ax2.title.set_text('Kinetic energy')
ax2.plot(kinetic_energy)

ax3 = fig.add_subplot(223)
ax3.title.set_text('Potential energy')
ax3.plot(potential_energy)

ax4 = fig.add_subplot(224)
ax4.title.set_text('Total energy')
ax4.plot(total_energy)

plt.tight_layout()

#plt.show()
plt.savefig('plot_remote.png')
