import numpy as np
import matplotlib.pyplot as plt
from MD_System_Class import Simulated_System as md


# Initialise MD system:
simulation = md(lx=20, rho=0.1, sigma=1, T=0.1)
#simulation.map()
simulation.update_neighbour_list(simulation.Particles[4])
print(len(simulation.Particles[4].neighbour_list))
