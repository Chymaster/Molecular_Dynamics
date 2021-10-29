import numpy as np
from MD_Particle_Class import Particle


class Simulated_System:
    def __init__(self, lx=20, rho=0.1, sigma=1, T=0.1):

        self.sigma = sigma
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
            self.Particles.append(Particle(positions[i], velocities[i], i))

    # Mapping the Particles in grid
    def map(self):
        import matplotlib.pyplot as plt
        plt.scatter(*np.array([i.position for i in self.Particles]).T)
        plt.show()

    # Energy of each particle calculating
    def update_force(self, particle_instance, epsilon=1, alpha=0, r_cut=2.5):

        # Neighbour list updating
        def update_neighbour_list():
            particle_instance.neighbour_list = []
            for particles in self.Particles:
                if np.linalg.norm(particle_instance.position
                                  - particles.position) < r_cut and particle_instance.serial != particles.serial:
                    particle_instance.neighbour_list.append(particles)

        # Calculate force according to LJ
        def calc_force(neighbours):
            sigma = self.sigma
            r = np.linalg.norm(neighbours.position-particle_instance.position)
            force = 4*epsilon * (((sigma**12)/(12*(r**13)))
                                 - alpha*(sigma**6)/(6*(r**7)))
            r_hat = (neighbours.position-particle_instance.position)/r
            force = force*r_hat
            return force

        update_neighbour_list()
        print(particle_instance.position)
        print([i.position for i in particle_instance.neighbour_list])

        # Calculating forces and assigning them to particle instance
        for neighbours in particle_instance.neighbour_list:
            particle_instance.force += calc_force(neighbours)


###############################################################################
# Region of testing
