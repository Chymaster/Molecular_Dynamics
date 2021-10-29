import numpy as np


class Particle:
    def __init__(self, position, velocity):
        self.position = np.array(position)
        self.velocity = np.array(velocity)

        # Neighbour list will be updated in MD_System_Class
        self.neighbour_list = []
