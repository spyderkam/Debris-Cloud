#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "June 2, 2025"

"""
    Simple simulation of a debris cloud. 
    Random event zero simulation.
"""

from src.gdmpidc import *
from src.gdmpidc_tools import *
from src.geometric_analysis import *
import time

start_time = time.time()

parent_mass = 100e3  # [kilograms]
parent_rad = 10      # [meters]
cloud = Cloud(parent_mass, parent_rad)

print(f"\nComputational PID cloud creation time: {time.time() - start_time:.2f} seconds\n")

def pretty_output(p):
    return f"({p[0]:.2f}, {p[1]:.2f}, {p[2]:.2f})"

print(f"  • A mass of {9.38e4} kg has been broken up into {len(cloud.all_points)} (relavent) fragments between 0.001 and 1.0 meters")
print(f"  • The radius of the debris cloud at collision is {cloud.radius:.2f} meters")

p1, p2 = get_entry_exit(cloud.radius, center=(0,0,0), diameter=False)
line = line_parametric_3d(p1, p2)
hit_distance = 1.0

print(f"  • The number of fragments within {hit_distance} meters of a line passing through {pretty_output(p1)} and {pretty_output(p2)} is ⏱︎", end="\r")
point_count_start_time = time.time()
nHits = count_points_near_line(line, cloud.all_points, hit_distance)
point_count_end_time = time.time()

print(f"  • The number of fragments within {hit_distance} meters of a line passing through {pretty_output(p1)} and {pretty_output(p2)} is \033[4m{nHits}\033[0m")
print(f"      ◦ This is {round(nHits/len(cloud.all_points)*100, 2)}% of the total fragments within the cloud")
print(f"      ◦ Computation time for point count: {round(point_count_end_time - point_count_start_time, 2)} seconds\n")
