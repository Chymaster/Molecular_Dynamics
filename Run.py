import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md

# Name of the log files
name = "dt001"
# Alpha values
alpha = 1
# dt values
delta_t = 0.01


# Initialise MD system:
simulation = md(lx=20, rho=0.05, sigma=1, T=5, dt=delta_t, r_cut=2.5)


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
    simulation.move(m=1, alpha=alpha, log_file="", f_log=0.0001)
    if i == 100:
        velocities = np.array([np.linalg.norm(i.velocity) for i in simulation.Particles])
        ms_displacement = np.mean((velocities/simulation.dt)**2)
        D = ms_displacement/(4*simulation.dt)
        output = open(name+".txt", "w")
        output.write("Diffusion coefficient at final stage is "+str(D))
        output.close()



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

# Saving figures
# Name of the parameter
plt.savefig(name+'.png')

fig = plt.figure()
plt.hist(velocities, density=True, bins=30)
plt.ylabel('Probability')
plt.xlabel('Velocity')
plt.savefig(name+"_velocity_hist.png")
