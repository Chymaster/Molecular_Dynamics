import numpy as np
from MD_Particle_Class import Particle


class Simulated_System:
    def __init__(self, lx=20, rho=0.1, sigma=1, T=0.1):

        # Calculate number of particles
        number_of_particles = int((lx**2)*rho)

        # Generate position lists:
        positions = []

        def overlaping(suggested_position):
            for existing_position in positions:
                if np.linalg.norm(suggested_position
                                  - existing_position) < sigma*2**(1/6):
                    return True
            return False

        def random_position():
            return np.random.random(size=2)*lx

        # Generating positions
        for i in range(number_of_particles):
            suggested_position = random_position()
            while overlaping(suggested_position):
                suggested_position = random_position()
            positions.append(suggested_position)

        # Generate velocities lists
        velocities = np.random.random((number_of_particles, 2))
        vx_average = np.average(velocities.T[0])
        vy_average = np.average(velocities.T[1])
        v_average_squared = np.average(np.reshape(velocities, -1))**2
        fs = np.sqrt(3*T/v_average_squared)
        for i in range(len(velocities)):
            velocities[i] = [(velocities[i][0]-vx_average)*fs,
                             (velocities[i][1]-vy_average)*fs]

        # Create particle instances
        self.Particles = []
        for i in range(number_of_particles):
            self.Particles.append(Particle(positions[i], velocities[i]))

    # Mapping the Particles in grid
    def map(self):
        import matplotlib.pyplot as plt
        plt.scatter(*np.array([i.position for i in self.Particles]).T)
        plt.show()


###############################################################################
# Region of testing

    def update_neighbour_list(self, particle_instance, r_cut=2.5):
        particle_instance.neighbour_list = []
        for particles in self.Particles:
            if np.linalg.norm(particle_instance.position
                              - particles.position) < r_cut:
                particle_instance.neighbour_list.append(particles)
