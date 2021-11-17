import numpy as np
from MD_Particle_Class import Particle


class Simulated_System:
    def __init__(self, lx=20, rho=0.1, sigma=1, T=0.1, dt=0.01):

        self.sigma = sigma

        self.size = lx

        self.t = 0

        self.dt = dt

        self.potential_E = 0.

        self.kinetic_E = 0.

        self.E = 0.

        self.T = T

        # Calculate number of particles
        number_of_particles = int((lx**2)*rho)
        self.N = number_of_particles

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
            self.Particles.append(Particle(positions[i], velocities[i], i, dt))

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
    def dist(self, t_particle, o_particle):
        """test"""
        """if max(self.p_diff(t_particle, o_particle)) > 5.:
            print("t", t_particle.position, "o", o_particle.position)
            print("r", self.p_diff(t_particle, o_particle))
            print("distance", np.linalg.norm(
                self.p_diff(t_particle, o_particle)))"""
        """test"""
        return np.linalg.norm(self.p_diff(t_particle, o_particle))

    def p_diff(self, t_particle, o_particle):
        r = np.array([0., 0.])

        r[0] = o_particle.position[0] - t_particle.position[0]
        if abs(self.size - abs(o_particle.position[0] - t_particle.position[0])) < abs(r[0]):
            r[0] = (self.size - abs(o_particle.position[0] - t_particle.position[0]))
        r[1] = o_particle.position[1] - t_particle.position[1]
        if abs(self.size - abs(o_particle.position[1] - t_particle.position[1])) < abs(r[1]):
            r[1] = (self.size - abs(o_particle.position[1] - t_particle.position[1]))

        return r
    # Energy of each particle calculating

    def update_force(self, target_particle, epsilon=1, alpha=0, r_cut=2.5):

        # Neighbour list updating
        def update_neighbour_list():
            target_particle.neighbour_list = []
            for other_particles in self.Particles:
                if self.dist(target_particle, other_particles) < r_cut and target_particle.serial != other_particles.serial:
                    target_particle.neighbour_list.append(other_particles)
                    """test""" """this is CHEATING!!!"""
                    # Making sure particles don't overlap
                    while self.dist(target_particle, other_particles) <= 0.2:
                        other_particles.position += (
                            other_particles.position - target_particle.position)
                    """test"""
        # Calculate force according to LJ

        def calc_force(neighbours):
            change_in_position = self.p_diff(target_particle, neighbours)
            sigma = self.sigma
            r = self.dist(neighbours, target_particle)
            twelve_term = 2*(sigma**12) / (r**13)
            six_term = alpha * (sigma**6)/(r**7)

            force = (-1) * 24 * epsilon * (twelve_term - six_term)
            r_hat = change_in_position / \
                np.sqrt(np.sum((change_in_position)**2))
            force_vector = force*r_hat

            """test"""
            # Checking force calculation
            """if np.linalg.norm(force) > 100000:
                print("epsilon=", epsilon, "twelve_term",
                      twelve_term, "six_term", six_term)
                print("target position", target_particle.position)
                print("neighbours position", neighbours.position)
                print("r", r)
                print("force_vector", force_vector)
                print("change_in_position", change_in_position)
                quit()"""

            return force_vector

        def calc_energy(neighbours):
            change_in_position = self.p_diff(target_particle, neighbours)
            sigma = self.sigma
            r = self.dist(neighbours, target_particle)
            twelve_term = (sigma / r)**12
            six_term = alpha * (sigma / r)**7

            energy = 4 * epsilon * (twelve_term - six_term)
            return energy

        update_neighbour_list()

        # Calculating forces and assigning them to particle instance
        target_particle.force = np.array([0., 0.])
        target_particle.energy = 0
        for neighbours in target_particle.neighbour_list:
            """test"""
            """if np.linalg.norm(calc_force(neighbours)) > 10:
                print("force", calc_force(neighbours))
                print("target particle position", target_particle.position)
                print("neighhbours", neighbours.position)
                print("distance", np.linalg.norm(
                    neighbours.position - target_particle.position))
                quit()"""
            """test"""

            target_particle.force += calc_force(neighbours)
            target_particle.energy += calc_energy(neighbours)

    def move(self, m=1, alpha=0, log_file="", f_log=0.0001):
        dt = self.dt
        for particle in self.Particles:

            self.update_force(particle, alpha=alpha)

            """Test"""
            """# If the last position is too large
            if np.linalg.norm(particle.force) > 10:
                print("current position = ", particle.position)
                print("force=", particle.force)
                print("last position = ", particle.position_last)
                print("velocity=", particle.velocity)
                print("neighbours", [
                      i.position for i in particle.neighbour_list])
                print("distance between neighbours", [np.linalg.norm(
                    i.position - particle.position_last) for i in particle.neighbour_list])
                self.map(particles=[particle],
                         neighbours=particle.neighbour_list)
                quit()"""
            """/Test"""

            # Set current position
            current_position = np.array([0., 0.])
            current_position = particle.position

            # Position change
            new_position = 2*particle.position - \
                particle.position_last + particle.force*(dt**2)/m

            particle.position = new_position
            # Velocity change
            particle.velocity = (particle.position
                                 - particle.position_last) / (2*dt)
            """test"""
            """if np.linalg.norm(particle.velocity) > 10:
                print("current position = ", particle.position)
                print("force=", particle.force)
                print("last position = ", particle.position_last)
                print("velocity=", particle.velocity)
                print("neighbours", [
                      i.position for i in particle.neighbour_list])"""
            """end test"""
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
            particle.position_last = current_position

            # Update time
            self.t += dt

            # Update potential energy
            self.potential_E = sum(i.energy for i in self.Particles)

            # Update kinetic Energy
            self.kinetic_E = 0.5 * m * \
                sum((np.linalg.norm(i.velocity))**2 for i in self.Particles)

            # Update total Energy
            self.E = self.potential_E + self.kinetic_E

            # Update Temperature
            kb = 1.38e-23
            self.T = self.kinetic_E / (self.N)

            # Logging
            if log_file is not "" and int(self.t / self.dt) % int(1/f_log) == 0:
                log = open(log_file, "a")
                log.write("".join(str(i)+" " for i in [self.t, self.T, self.kinetic_E,
                          self.potential_E, self.E, "\n"]))
                log.close()

###############################################################################
# Region of testing #
###############################################################################
    def gro_log(self, gro_file="MD.gro", f_traj=0.01):
        gro = open(gro_file)
        if int(self.t / self.dt) % int(1/f_traj) == 0:
            gro.write()
