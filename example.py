from MD_System_Class import Simulated_System as md

# Initialise MD system:
simulation = md(lx=20, rho=0.05, sigma=1, T=5, dt=0.001, r_cut=2.5)

# Making N moves
N = 1000

for i in range(N):
        simulation.move(m=1, alpha=1, log_file="", f_log=0.0001)

simulation.map()
