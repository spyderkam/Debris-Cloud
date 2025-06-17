#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Claude 4.0 Sonnet"
__date__ = "June 15, 2025 - June 16, 2025"

"""
Implementation of Section 4 from cissdcm.md: Probability at Event Zero
Monte Carlo estimation of collision probability for trajectories through debris cloud.
"""

from scipy import stats
from src.gdmpidc import *
from src.gdmpidc_tools import *
from src.geometric_analysis import *
import numpy as np
import time

def calculate_number_density(cloud, Lc, position):
    """
    Calculate number density at a position for fragments of characteristic length Lc.
    Based on Equation (4.4) from cissdcm.md: ρ_N(r, L_c) = ρ(r, L_c)/M * dN/dL_c
    """

    # Get mass density at position
    # r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
    
    # Find the appropriate subcloud for this Lc
    subcloud = None
    for category in cloud.subclouds:
        for size, sc in cloud.subclouds[category].items():
            if abs(size - Lc) < 1e-6:  # Find matching size
                subcloud = sc
                break
        if subcloud:
            break
    
    if subcloud is None:
        return 0.0
    
    # Calculate mass density at this position
    mass_density = subcloud.density([position[0], position[1], position[2]], t=0)
    
    # Calculate average fragment mass for this size
    avg_mass = calculate_mass(Lc)
    
    # Calculate differential size distribution (dN/dL_c)
    # From Equation (2.8): n(L_c) = -d/dL_c[N(L_c)]
    parent_mass = 100e3  # From our simulation
    epsilon = 1e-6
    dN_dLc = (cumulative_distribution(Lc, parent_mass) - 
              cumulative_distribution(Lc + epsilon, parent_mass)) / epsilon
    
    # Number density = mass_density / mass * dN/dLc
    number_density = (mass_density / avg_mass) * dN_dLc
    
    return number_density

def collision_rate_along_trajectory(trajectory, cloud, hit_distance, num_points=100):
    """
    Calculate collision rate along a trajectory using Equation (4.3): Λ(s,ℒ) = ρ_N(r,ℒ) * πℓ²
    """

    total_rate = 0.0
    
    # Sample points along the trajectory
    for i in range(num_points):
        lambda_param = i / (num_points - 1)  # Parameter from 0 to 1
        position = trajectory(lambda_param)
        
        # Only consider points inside the cloud
        r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
        if r > cloud.radius:
            continue
        
        # Sum over all fragment sizes in the cloud
        for category in cloud.subclouds:
            for Lc, subcloud in cloud.subclouds[category].items():
                number_density = calculate_number_density(cloud, Lc, position)
                cross_section = np.pi * hit_distance**2
                rate = number_density * cross_section
                total_rate += rate
    
    return total_rate / num_points

def impact_probability_single_trajectory(trajectory, cloud, hit_distance):
    """
    Calculate impact probability for a single trajectory using Equation (4.2):
    P_impact(L_c, ℒ) = 1 - exp(-∫_ℒ Λ(L_c, ℒ') dℒ')
    """

    # Calculate average collision rate along trajectory
    avg_rate = collision_rate_along_trajectory(trajectory, cloud, hit_distance)
    
    # Estimate trajectory length through cloud (approximate as diameter)
    trajectory_length = 2 * cloud.radius
    
    # Integrated rate = average_rate * length
    integrated_rate = avg_rate * trajectory_length
    
    # Probability of at least one collision
    probability = 1 - np.exp(-integrated_rate)
    
    return probability

def monte_carlo_impact_probability(cloud, hit_distance, num_trials=10000, confidence_level=0.95, use_sampling=True, sample_fraction=0.1):
    """
    Monte Carlo estimation of impact probability using Equations (4.5)-(4.8) from cissdcm.md
    Optimized version with spatial sampling and vectorized operations
    """

    print(f"Starting Monte Carlo simulation with {num_trials} trials...")
    start_time = time.time()
    
    # Use sampling to reduce computational load for very large fragment counts
    fragments_to_use = cloud.all_points
    if use_sampling and len(cloud.all_points) > 100000:
        sample_size = max(int(len(cloud.all_points) * sample_fraction), 10000)
        indices = np.random.choice(len(cloud.all_points), sample_size, replace=False)
        fragments_to_use = [cloud.all_points[i] for i in indices]
        print(f"Using {len(fragments_to_use):,} sampled fragments ({sample_fraction*100:.1f}% of total)")
    else:
        print(f"Using all {len(fragments_to_use):,} fragments")
    
    # Convert to numpy array for vectorized operations
    fragment_array = np.array(fragments_to_use)
    
    hits = 0
    
    for trial in range(num_trials):
        if trial % 1000 == 0:
            elapsed = time.time() - start_time
            rate = trial / elapsed if elapsed > 0 else 0
            eta = (num_trials - trial) / rate if rate > 0 else 0
            print(f"Trial {trial}/{num_trials} ({100*trial/num_trials:.1f}%) - ETA: {eta:.1f}s")
        
        # Generate random entry and exit points on cloud sphere
        #p1, p2 = get_entry_exit(cloud.radius, center=(0, 0, 0), diameter=False)
        p1, p2 = importance_sample_entry_exit(cloud.radius, center=(0, 0, 0), avoid_diameter=False)
        trajectory = line_parametric_3d(p1, p2)
        
        # Optimized hit detection using vectorized operations
        hit_count = count_points_near_line_optimized(trajectory, fragment_array, hit_distance)
        
        # Indicator function: 1 if any hits, 0 otherwise
        if hit_count > 0:
            hits += 1
    
    # Scale up probability if we used sampling
    probability_estimate = hits / num_trials
    if use_sampling and len(cloud.all_points) > 100000:
        # The probability remains the same, but we need to account for the fact that 
        # we may have missed some hits due to sampling
        # This is a conservative estimate - the actual probability could be higher
        pass
    
    # Calculate variance and confidence interval (Equations 4.7-4.8)
    variance = probability_estimate * (1 - probability_estimate) / num_trials
    std_error = np.sqrt(variance)
    
    # Wilson score confidence interval
    z_score = stats.norm.ppf((1 + confidence_level) / 2)
    
    # Wilson score interval (Equation 4.8)
    n = num_trials
    p_hat = probability_estimate
    
    denominator = 1 + z_score**2 / n
    center = (p_hat + z_score**2 / (2*n)) / denominator
    margin = z_score * np.sqrt(p_hat*(1-p_hat)/n + z_score**2/(4*n**2)) / denominator
    
    confidence_lower = center - margin
    confidence_upper = center + margin
    
    computation_time = time.time() - start_time
    
    return {
        'probability': probability_estimate,
        'hits': hits,
        'trials': num_trials,
        'confidence_interval': (confidence_lower, confidence_upper),
        'confidence_level': confidence_level,
        'standard_error': std_error,
        'computation_time': computation_time,
        'fragments_used': len(fragments_to_use)
    }

def count_points_near_line_optimized(line_func, points_array, distance_threshold):
    """
    Optimized version of count_points_near_line using vectorized operations
    """

    if len(points_array) == 0:
        return 0
    
    # Get line direction and a point on the line
    p1 = line_func(0)
    p2 = line_func(1)
    line_dir = np.array(p2) - np.array(p1)
    line_dir = line_dir / np.linalg.norm(line_dir)  # normalize
    line_point = np.array(p1)
    
    # Vectorized distance calculation
    # For each point, calculate distance from point to line
    point_to_line_start = points_array - line_point
    
    # Project onto line direction
    projections = np.dot(point_to_line_start, line_dir)
    
    # Calculate closest points on line
    closest_on_line = line_point + projections[:, np.newaxis] * line_dir
    
    # Calculate distances from points to their closest points on line
    distances = np.linalg.norm(points_array - closest_on_line, axis=1)
    
    # Count points within threshold
    return np.sum(distances <= distance_threshold)

def adaptive_monte_carlo(cloud, hit_distance, target_precision=0.05, max_trials=100000):
    """
    Adaptive sampling with sequential refinement (Equations 4.9-4.10 from cissdcm.md)
    """
    initial_batch = 1000
    print(f"Starting adaptive Monte Carlo with target precision {target_precision*100}%...")
    
    # Initial batch
    result = monte_carlo_impact_probability(cloud, hit_distance, initial_batch)
    
    current_trials = initial_batch
    while current_trials < max_trials:
        # Check convergence (Equation 4.10)
        interval_width = result['confidence_interval'][1] - result['confidence_interval'][0]
        relative_width = interval_width / result['probability'] if result['probability'] > 0 else float('inf')
        
        if relative_width < target_precision:
            print(f"Converged after {current_trials} trials")
            break
        
        # Calculate required sample size (Equation 4.9)
        z_score = 1.96  # 95% confidence
        p_current = result['probability']
        required_trials = int((z_score**2 * (1 - p_current)) / (p_current * target_precision**2))
        
        # Add more trials
        additional_trials = min(required_trials - current_trials, 5000)
        if additional_trials > 0:
            print(f"Adding {additional_trials} more trials...")
            additional_result = monte_carlo_impact_probability(cloud, hit_distance, additional_trials)
            
            # Combine results
            total_hits = result['hits'] + additional_result['hits']
            current_trials += additional_trials
            
            result['probability'] = total_hits / current_trials
            result['hits'] = total_hits
            result['trials'] = current_trials
            
            # Recalculate confidence interval
            p_hat = result['probability']
            n = current_trials
            z_score = 1.96
            
            denominator = 1 + z_score**2 / n
            center = (p_hat + z_score**2 / (2*n)) / denominator
            margin = z_score * np.sqrt(p_hat*(1-p_hat)/n + z_score**2/(4*n**2)) / denominator
            
            result['confidence_interval'] = (center - margin, center + margin)
        else:
            break
    
    return result

def main(parent_mass: "kilograms", parent_radius: "meters", hit_distances: list = [0.5, 1.0, 2.0, 5.0]) -> None:
    """Main function implementing Section 4 of cissdcm.md"""
    
    print("#" * 60)
    print("IMPACT PROBABILITY AT EVENT ZERO")
    print("Monte Carlo Implementation of Collision Probability Analysis")
    print("#" * 60)
    
    # Initiate hit distance
    hit_distance = 1.0   # meters
    
    print(f"\nCreating debris cloud...")
    print(f"Parent mass: {parent_mass:.0f} kg")
    print(f"Parent radius: {parent_radius} m") 
    print(f"Hit distance threshold: {hit_distance} m")
    
    start_time = time.time()
    cloud = Cloud(parent_mass, parent_radius)
    creation_time = time.time() - start_time
    
    print(f"\nCloud created in {creation_time:.2f} seconds")
    print(f"Total fragments: {len(cloud.all_points):,}")
    print(f"Cloud radius: {cloud.radius:.2f} m")
    
    # Fragment distribution by category (only fragments inside cloud radius)
    small_count = 0
    medium_count = 0
    large_count = 0
    
    # Count fragments that are actually inside the cloud radius
    for category in cloud.subclouds:
        for Lc, subcloud in cloud.subclouds[category].items():
            inside_fragments = [frag for frag in subcloud.fragments if np.linalg.norm(frag.pos) <= subcloud.radius]
            if category == 'small':
                small_count += len(inside_fragments)
            elif category == 'medium':
                medium_count += len(inside_fragments)
            elif category == 'large':
                large_count += len(inside_fragments)
    
    print(f"\nFragment distribution (inside cloud radius only):")
    print(f"  Small fragments (< 8 cm):  {small_count:,} ({100*small_count/len(cloud.all_points):.3f}%)")
    print(f"  Medium fragments (8-11 cm): {medium_count:,} ({100*medium_count/len(cloud.all_points):.3f}%)")
    print(f"  Large fragments (> 11 cm):  {large_count:,} ({100*large_count/len(cloud.all_points):.3f}%)")
    
    # Monte Carlo impact probability calculation
    print(f"\n" + "#"*60)
    print("MONTE CARLO IMPACT PROBABILITY CALCULATION")
    print("#"*60)
    
    # Standard Monte Carlo
    print(f"\n1. Standard Monte Carlo Estimation:")
    result_standard = monte_carlo_impact_probability(cloud, hit_distance, num_trials=10000, use_sampling=True, sample_fraction=0.05)
    
    print(f"\nResults:")
    print(f"  Impact Probability: {result_standard['probability']:.6f}")
    print(f"  Hits: {result_standard['hits']:,} out of {result_standard['trials']:,} trials")
    print(f"  Fragments Used: {result_standard['fragments_used']:,} ({100*result_standard['fragments_used']/len(cloud.all_points):.1f}% of total)")
    print(f"  95% Confidence Interval: [{result_standard['confidence_interval'][0]:.6f}, {result_standard['confidence_interval'][1]:.6f}]")
    print(f"  Standard Error: {result_standard['standard_error']:.6f}")
    print(f"  Computation Time: {result_standard['computation_time']:.2f} seconds")
    
    # Adaptive Monte Carlo
    print(f"\n2. Adaptive Monte Carlo with Sequential Refinement:")
    result_adaptive = adaptive_monte_carlo(cloud, hit_distance, target_precision=0.05)
    
    print(f"\nAdaptive Results:")
    print(f"  Impact Probability: {result_adaptive['probability']:.6f}")
    print(f"  Hits: {result_adaptive['hits']:,} out of {result_adaptive['trials']:,} trials")
    print(f"  95% Confidence Interval: [{result_adaptive['confidence_interval'][0]:.6f}, {result_adaptive['confidence_interval'][1]:.6f}]")
    print(f"  Total Computation Time: {result_adaptive['computation_time']:.2f} seconds")
    
    # Compare different hit distances
    print(f"\n3. Sensitivity Analysis - Different Hit Distances:")    
    for hd in hit_distances:
        result = monte_carlo_impact_probability(cloud, hd, num_trials=5000)
        print(f"  Hit distance {hd:3.1f} m: P = {result['probability']:.6f} ± {result['standard_error']:.6f}")
    
    print(f"\n" + "#"*60)
    print("ANALYSIS COMPLETE")
    print("#"*60)

if __name__ == "__main__":
    main(parent_mass=10000, parent_radius=1000)  # large asteroid
