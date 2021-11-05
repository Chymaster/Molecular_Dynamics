import numpy as np
from MD_Particle_Class import Particle


class Simulated_System:
    def __init__(self, lx=20, rho=0.1, sigma=1, T=0.1):

        self.sigma = sigma

        self.size = lx

        self.t = 0

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
        v_average_squared = np.average(np.reshape(velocities, -1)**2)
        fs = np.sqrt(3*T/v_average_squared)
        for i in range(len(velocities)):
            velocities[i] = [(velocities[i][0]-vx_average)*fs,
                             (velocities[i][1]-vy_average)*fs]

        # Create particle instances
        self.Particles = []
        for i in range(number_of_particles):
            self.Particles.append(Particle(positions[i], velocities[i], i))

    # Mapping the Particles in grid
    def map(self, particles=None, neighbours=None):
        if particles is None:
            particles = self.Particles
        import matplotlib.pyplot as plt
        plt.scatter(*np.array([i.position for i in particles]).T, color="red")
        if neighbours is not None:
            plt.scatter(
                *np.array([i.position for i in neighbours]).T, color="blue")
            plt.scatter(
                *np.array([i.position for i in self.Particles]).T, color="green", alpha=0.3)

        plt.xlim([0, self.size])
        plt.ylim([0, self.size])
        plt.show()

    # Distance between two particles, with boundary condition applied
    def dist(self, particle1, particle2):
        r = np.linalg.norm(particle1.position - particle2.position)
        r_boundary_x = np.linalg.norm(
            abs(particle1.position - particle2.position) - np.array([self.size, 0]))
        r_boundary_y = np.linalg.norm(
            abs(particle1.position - particle2.position) - np.array([0, self.size]))
        r_boundary_xy = np.linalg.norm(
            abs(particle1.position - particle2.position) - np.array([self.size, self.size]))
        return min(r, r_boundary_x, r_boundary_y, r_boundary_xy)

    # Energy of each particle calculating

    def update_force(self, target_particle, epsilon=1, alpha=0, r_cut=2.5):

        # Neighbour list updating
        def update_neighbour_list():
            target_particle.neighbour_list = []
            for other_particles in self.Particles:
                if self.dist(target_particle, other_particles) < r_cut and target_particle.serial != other_particles.serial:
                    target_particle.neighbour_list.append(other_particles)
                    """# Making sure particles don't overlap
                    while self.dist(target_particle, other_particles) <= 0.4:
                        other_particles.position += (
                            other_particles.position - target_particle.position)"""
        # Calculate force according to LJ

        def calc_force(neighbours):
            sigma = self.sigma
            r = self.dist(neighbours, target_particle)
            force = - 24*epsilon * (2*((sigma**12)/(r**13))
                                    - alpha*(sigma**6)/(r**7))
            r_hat = (neighbours.position-target_particle.position) / \
                np.sqrt(np.sum((neighbours.position-target_particle.position)**2))
            force = force*r_hat
            """Test"""
            target_particle.distance = min(
                abs(r), abs(target_particle.distance))
            """/Test"""

            return force

        update_neighbour_list()

        # Calculating forces and assigning them to particle instance
        for neighbours in target_particle.neighbour_list:
            target_particle.force += calc_force(neighbours)


###############################################################################
# Region of testing #
###############################################################################


    def move(self, dt, m=1):
        for particle in self.Particles:
            """test"""
            # Periodic boundary condition
            while particle.position[0] < 0:
                particle.position[0] += self.size
            while particle.position[1] < 0:
                particle.position[1] += self.size
            while particle.position[0] > self.size:
                particle.position[0] -= self.size
            while particle.position[1] > self.size:
                particle.position[1] -= self.size
            """/test"""

            self.update_force(particle)
            # Position change
            if np.linalg.norm(particle.position_last) > self.size**2:
                print("last position before moving = ", particle.position_last)
            particle.position = 2*particle.position - \
                particle.position_last + particle.force*(dt**2)/m
            # Velocity change
            particle.velocity = (particle.position
                                 - particle.position_last) / 2*dt
            if np.linalg.norm(particle.velocity) > 2:
                print("current position = ", particle.position)
                print("force=", particle.force)
                print("last position = ", particle.position_last)
                print("velocity=", particle.velocity)
            # Periodic boundary condition
            while particle.position[0] < 0:
                particle.position[0] += self.size
            while particle.position[1] < 0:
                particle.position[1] += self.size
            while particle.position[0] > self.size:
                particle.position[0] -= self.size
            while particle.position[1] > self.size:
                particle.position[1] -= self.size

            # Update last position
            particle.position_last = particle.position
