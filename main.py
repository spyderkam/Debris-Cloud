#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 20, 2025"

from src.gdmpidc import *
import numpy as np


cloud = Cloud(characteristic_length=0.05, num_fragments=1000, breakup_type="collision")
all_initial_positions = cloud.sample_positions()
inside_points = [point for point in all_initial_positions if np.linalg.norm(point) <= cloud.radius]

def cloud_intersect_line(entry_point, exit_point, t):
    """
    Given the entry and exit points of the cloud at time 𝑡, compute the slope 
    and 𝑦-intercept of the line that passes through these points.

    Args:
        entry_point (list): (x, y, z) coordinates of the entry point
        exit_point (list): (x, y, z) coordinates of the exit point
        t (float): Time since impact [s]

    Returns:
        tuple: (slope, y_intercept)
    """
    



if __name__ == "__main__":
    pass
