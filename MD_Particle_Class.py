import numpy as np


class Particle:
    def __init__(self, position, velocity, serial):
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.serial = serial

        # Neighbour list will be updated in MD_System_Class
        self.neighbour_list = []

        # Force from other particles will be recorded
        self.force = np.array([0., 0.])

        # Last position
        self.position_last = self.position + self.velocity*0.1

################################################################################
# Region of test
################################################################################
        self.distance = 5.
