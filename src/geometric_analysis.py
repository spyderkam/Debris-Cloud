#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh, Claude 3.7 Sonnet, Grok 3, Claude 3.5 Sonnet V2"
__date__ = "May 20, 2025 - May 21, 2025"

try:
    # When imported as a module from parent directory
    from src.gdmpidc import *
except ImportError:
    # When run directly from src directory
    from gdmpidc import *

def get_entry_exit(radius, center=(0, 0, 0), diameter=False):
    """
    Generate two random points (entry and exit) on the surface of a sphere.
    
    Parameters:
    - radius (float): The radius of the sphere.
    - center (tuple): Center of the sphere in (𝑥,𝑦,𝑧) coordinates. Default is (0,0,0).
    - diameter (bool): If True, the points will be diametrically opposite.
                       If False, the points will not be diametrically opposite. Default is False.
    
    Returns:
    - tuple: ((𝑥1,𝑦1,𝑧1), (𝑥2,𝑦2,𝑧2)) - entry and exit points on the sphere's surface.
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
    Compute the parametric equation of the 3D line through points 𝑝1 and 𝑝2.
    
    Args:
        𝑝1 (tuple): First point as (𝑥1,𝑦1,𝑧1)
        𝑝2 (tuple): Second point as (𝑥2,𝑦2,𝑧2)
    
    Returns:
        callable: Function evaluate (𝜆) that returns the point (𝑥,𝑦,𝑧) on the line at parameter 𝜆
    """
    
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    
    def evaluate(λ):
        x = x1 + λ*(x2 - x1)
        y = y1 + λ*(y2 - y1)
        z = z1 + λ*(z2 - z1)
        return (x, y, z)
    
    return evaluate

def count_points_near_line(line, points, distance) -> int:
    """
    Count number of points that are within given distance from a parametric line.
    
    Args:
        line (callable): Parametric line function 𝑟(𝜆) = 𝑝1 + 𝜆(𝑝2-𝑝1)
        points (np.ndarray): Array of points to check, shape (𝑁,3)
        distance (float): Maximum distance from line
        
    Returns:
        int: Number of points within the distance from the line
    """
    
    # Initialize count
    count = 0
    
    # For each point
    for point in points:
        # Find λ that minimizes distance to line by projecting point onto line
        # Let k = p2-p1, then λ = (point-p1)·k / |k|^2
        p1 = line(0)  # Get p1 by evaluating at λ=0
        p2 = line(1)  # Get p2 by evaluating at λ=1
        k = tuple(b-a for a,b in zip(p1,p2))
        k_mag_sq = sum(x*x for x in k)
        
        # Vector from p1 to point
        q = tuple(b-a for a,b in zip(p1,point))
        
        # Projection coefficient
        λ = sum(a*b for a,b in zip(q,k)) / k_mag_sq
        
        # Get closest point on line
        closest = line(λ)
        
        # Calculate distance from point to line
        dist = sum((a-b)**2 for a,b in zip(point,closest))**0.5
        
        if dist <= distance:
            count += 1
            
    return count


if __name__ == "__main__":
    cloud = SubCloud(characteristic_length=0.05, parent_rad=4.1439e+11, num_fragments=1000)
    
    all_points = cloud.sample_positions()
    inside_points = [point for point in all_points if np.linalg.norm(point) <= cloud.radius]
    inside_points = np.array(inside_points)

    p1, p2 = get_entry_exit(cloud.radius, diameter=False)
    line = line_parametric_3d(p1, p2)

    x = count_points_near_line(line, inside_points, 1)
    print(x)
