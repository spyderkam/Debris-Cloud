#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 20, 2025"

from src.gdmpidc import *
import numpy as np


def get_entry_exit(radius, center=(0, 0, 0), diameter=False):
    """
    Generate two random points (entry and exit) on the surface of a sphere.
    
    Parameters:
    - radius (float): The radius of the sphere.
    - center (tuple): Center of the sphere in (x, y, z) coordinates. Default is (0, 0, 0).
    - diameter (bool): If True, the entry and exit points will be the endpoints of a diameter
                      (i.e., they will be diametrically opposite).
                      If False, the entry and exit points will not be the endpoints of a diameter
                      (i.e., they will not be diametrically opposite). Default is False.
    
    Returns:
    - tuple: ((x1, y1, z1), (x2, y2, z2)) - entry and exit points on the sphere's surface.
    """
    
    # Generate a random point on the unit sphere using the Gaussian method
    # This creates a uniform distribution of points on the sphere
    vec1 = np.random.normal(0, 1, 3)
    vec1 = vec1 / np.linalg.norm(vec1)
    
    if diameter:
        # For diameter, the opposite point is simply negated
        vec2 = -vec1
    else:
        # Generate another random point on the unit sphere
        # Make sure it's not diametrically opposite
        vec2 = np.random.normal(0, 1, 3)
        vec2 = vec2 / np.linalg.norm(vec2)
        
        # Check if the points are nearly diametrically opposite (dot product close to -1)
        dot_product = np.dot(vec1, vec2)
        
        # If the dot product is close to -1, the points are nearly diametrically opposite
        if dot_product < -0.9999:
            # Regenerate the second point until it's not diametrically opposite
            while dot_product < -0.9999:
                vec2 = np.random.normal(0, 1, 3)
                vec2 = vec2 / np.linalg.norm(vec2)
                dot_product = np.dot(vec1, vec2)
    
    # Scale to the desired radius and shift to the desired center
    entry_point = tuple(center[i] + radius * vec1[i] for i in range(3))
    exit_point = tuple(center[i] + radius * vec2[i] for i in range(3))
    
    return entry_point, exit_point


def cloud_intersect_line(entry_point, exit_point, t):
    """
    Given the entry and exit points of the cloud at time 𝑡, compute the slope 
    and 𝑦-intercept of the line that passes through these points.

    Args:
        entry_point (list): (𝑥,𝑦,𝑧) coordinates of the entry point
        exit_point (list): (𝑥,𝑦,𝑧) coordinates of the exit point
        t (float): Time since impact [s]

    Returns:
        tuple: (slope, y_intercept)
    """

    return None    
    



if __name__ == "__main__":
    cloud0 = Cloud(characteristic_length=0.05, num_fragments=10)
    #all_initial_positions = cloud.sample_positions()
    #inside_points = [point for point in all_initial_positions if np.linalg.norm(point) <= cloud.radius]

    cloud1 = Cloud(0.15, num_fragments=100)

    print(np.linalg.norm(get_entry_exit(cloud0.radius)[0]))
    print(np.linalg.norm(get_entry_exit(cloud1.radius)[0]))
    print(np.linalg.norm(get_entry_exit(cloud0.radius)[1]))
    print(np.linalg.norm(get_entry_exit(cloud1.radius)[1]))
