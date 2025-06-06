#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh, Claude 3.7 Sonnet, Grok 3, Claude 3.5 Sonnet V2"
__date__ = "May 20, 2025 - June 6, 2025"

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

def importance_sample_entry_exit(R_c, center=None, mu=0.6, sigma_IS=None, avoid_diameter=False, max_attempts=10000):
    """
    Generate entry/exit points biased toward trajectories passing near mu*R_c
    using rejection sampling.
    
    Parameters:
        - R_c: cloud radius
        - center: sphere center coordinates (default: [0,0,0])
        - mu: peak density location as fraction of R_c (default 0.6)
        - sigma_IS: importance sampling width (default: 0.2*R_c)
        - avoid_diameter: if True, ensure trajectory doesn't pass through sphere
        - max_attempts: maximum rejection sampling attempts
    
    Returns:
        - entry_point, exit_point: 3D coordinates on sphere
    """

    if center is None:
        center = np.array([0.0, 0.0, 0.0])
    else:
        center = np.array(center)
        
    if sigma_IS is None:
        sigma_IS = 0.2 * R_c
    
    peak_radius = mu * R_c
    max_importance = 1.0
    
    for _ in range(max_attempts):
        # Step 1: Generate uniform entry/exit points on unit sphere
        u1 = np.random.uniform(-1, 1)
        theta1 = np.arccos(u1)
        phi1 = np.random.uniform(0, 2*np.pi)
        entry_unit = np.array([
            np.sin(theta1) * np.cos(phi1),
            np.sin(theta1) * np.sin(phi1),
            np.cos(theta1)
        ])
        
        u2 = np.random.uniform(-1, 1)
        theta2 = np.arccos(u2)
        phi2 = np.random.uniform(0, 2*np.pi)
        exit_unit = np.array([
            np.sin(theta2) * np.cos(phi2),
            np.sin(theta2) * np.sin(phi2),
            np.cos(theta2)
        ])
        
        # Scale and translate to actual sphere
        entry = center + R_c * entry_unit
        exit_point = center + R_c * exit_unit
        
        # Check diameter constraint if needed
        if avoid_diameter:
            # Check if trajectory passes through the sphere
            # This happens when the closest point to center is inside the sphere
            direction = exit_point - entry
            direction_norm = direction / np.linalg.norm(direction)
            
            # Vector from entry to center
            entry_to_center = center - entry
            
            # Find closest point on trajectory to center
            t_closest = np.dot(entry_to_center, direction_norm)
            
            # Only check if closest point is between entry and exit
            if 0 < t_closest < np.linalg.norm(direction):
                closest_point = entry + t_closest * direction_norm
                dist_to_center = np.linalg.norm(closest_point - center)
                
                if dist_to_center < R_c:
                    continue  # Reject this sample, try again
        
        # Calculate minimum distance from trajectory to peak sphere
        # (Note: peak sphere is centered at 'center', not origin)
        entry_relative = entry - center
        exit_relative = exit_point - center
        direction = exit_relative - entry_relative
        direction_norm = direction / np.linalg.norm(direction)
        
        # Find closest point on trajectory to center (in relative coordinates)
        t_closest = -np.dot(entry_relative, direction_norm)
        closest_point_relative = entry_relative + t_closest * direction_norm
        
        # Distance from closest point to center
        dist_to_center = np.linalg.norm(closest_point_relative)
        
        # Minimum distance to peak sphere
        l_min = abs(dist_to_center - peak_radius)
        
        # Importance function value
        importance = np.exp(-l_min**2 / (2 * sigma_IS**2))
        
        # Accept/reject based on importance
        if np.random.uniform(0, 1) < importance / max_importance:
            return entry, exit_point
    
    # If we fail, return uniform sample (with diameter check)
    print(f"Warning: Importance sampling failed after {max_attempts} attempts")
    # Could fall back to uniform sampling here
    if avoid_diameter:
        # Use your existing get_entry_exit function as fallback
        return get_entry_exit(R_c, center, diameter=avoid_diameter)
    else:
        return entry, exit_point


if __name__ == "__main__":
    cloud = Cloud(100e3, 10,)
    p1, p2 = importance_sample_entry_exit(cloud.radius)
    line = line_parametric_3d(p1, p2)

    x = count_points_near_line(line, cloud.all_points, 1)
    print(round(x/len(cloud.all_points)*100, 2))
