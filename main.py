#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 20, 2025"

from src.gdmpidc import *
import numpy as np

inside_points = np.array([point for point in vec_r if np.linalg.norm(point) <= cloud.radius])

def func(entry_point, exit_point, t):
    """
    Given the entry and exit points of the sphere at time 𝑡, this function calculates the slope and 𝑦-intercept of the
    line that passes through these points.

    Args:
        entry_point (list): (x, y, z) coordinates of the entry point
        exit_point (list): (x, y, z) coordinates of the exit point
        t (float): Time since impact [s]

    Returns:
        tuple: (slope, y_intercept)
    """
    



if __name__ == "__main__":
    pass
