#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh, Claude 3.7 Sonnet, Grok 3"
__date__ = "May 20, 2025"

from src.gdmpidc import *
import numpy as np

def get_entry_exit(radius, center=(0, 0, 0), diameter=False):
    """
    Generate two random points (entry and exit) on the surface of a sphere.
    
    Parameters:
    - radius (float): The radius of the sphere.
    - center (tuple): Center of the sphere in (𝑥,𝑦,𝑧) coordinates. Default is (0,0,0).
    - diameter (bool): If True, the points will be diametrically opposite.
                      If False, the points will not be diametrically opposite. Default is False.
    
    Returns:
    - tuple: ((x1, y1, z1), (x2, y2, z2)) - entry and exit points on the sphere's surface.
    """
    
    # Convert center to numpy array for vectorized operations
    center = np.array(center)
    
    # Generate first random point on unit sphere using Gaussian method
    vec1 = np.random.normal(0, 1, 3)
    vec1 = vec1 / np.linalg.norm(vec1)
    
    if diameter:
        # For diameter, the opposite point is simply negated
        vec2 = -vec1
    else:
        # Keep generating random points until we get one that's not diametrically opposite
        while True:
            vec2 = np.random.normal(0, 1, 3)
            vec2 = vec2 / np.linalg.norm(vec2)
            
            # Check if points are not diametrically opposite
            # We use -0.98 as threshold rather than -1 to avoid numerical precision issues
            # This allows for more variation while still preventing diameter-like configurations
            if np.dot(vec1, vec2) > -0.98:
                break
    
    # Scale by radius and shift to center (vectorized)
    entry_point = tuple(center + radius * vec1)
    exit_point = tuple(center + radius * vec2)
    
    return entry_point, exit_point


def line_parametric_3d(p1, p2):
    """
    Compute the parametric equation of the 3D line through points p1 and p2.
    
    Args:
        p1 (tuple): First point as (x1, y1, z1)
        p2 (tuple): Second point as (x2, y2, z2)
    
    Returns:
        dict: Contains 'point' (p1), 'direction' (vector from p1 to p2), and
              a function 'evaluate(t)' that returns the point on the line at parameter t
    """
    
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    
    # Direction vector: (x2 - x1, y2 - y1, z2 - z1)
    direction = (x2 - x1, y2 - y1, z2 - z1)
    
    # Function to evaluate the line at parameter t
    def evaluate(t):
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        z = z1 + t * (z2 - z1)
        return (x, y, z)
    
    return {
        'point': p1,
        'direction': direction,
        'evaluate': evaluate
    }


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
    p1, p2 = get_entry_exit(cloud0.radius)

    # Example usage of line_parametric_3d
    print("p1 = ", p1)
    print(f"p2 = {p2}")
    line = line_parametric_3d(p1, p2)
    
    # Access components
    print(f"Starting point: {line['point']}")
    print(f"Direction vector: {line['direction']}")
    # Parametric equation: r(t) = (1, 2, 3) + t (3, 3, 3)
    
    # Evaluate at specific t values
    print(f"Point at t = 0: {line['evaluate'](0)}")  # Should return P1: (1, 2, 3)
    print(f"Point at t = 1: {line['evaluate'](1)}")  # Should return P2: (4, 5, 6)
    print(f"Point at t = 0.5: {line['evaluate'](0.5)}")  # Midpoint: (2.5, 3.5, 4.5)
