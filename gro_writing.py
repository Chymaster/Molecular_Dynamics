import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=20, rho=0.05, sigma=1, T=500, dt=0.001, r_cut=2.5)


# Making N moves
N = 1000

# Logging frequency (once every f_traj times)
f_traj = 1

file_name = "traj.gro"
log = open(file_name, "w")


def log_event(simulation):
    size = simulation.size
    n_particles = len(simulation.Particles)
    particle_properties = []
    time = simulation.t
    for i in simulation.Particles:
        particle_properties.append([str(i.position[0])[:6], str(i.position[1])[:6], "0.000", str(i.velocity[0])[:6], str(i.velocity[1])[:6], "0.000"])
    log.write("MD of "+str(size)+" particles, ")
    log.write("t="+str(time)[:5]+"\n")
    log.write(str(n_particles)[:].rjust(5)+"\n")
    for i in range(len(particle_properties)):
        log.write(str(i).rjust(5))
        log.write("LJ".rjust(5))
        log.write("C".rjust(5))
        log.write(str(i).rjust(5))
        for j in particle_properties[i]:
            log.write(j.rjust(8))
        log.write("\n")
    log.write(str(size).rjust(10))
    log.write(str(size).rjust(10))
    log.write("0.002".rjust(10)+"\n")

for i in range(N):
    simulation.move(m=1, alpha=1, log_file="", f_log=0.0001)
    if i%f_traj == 0:
        log_event(simulation)
log.close()
