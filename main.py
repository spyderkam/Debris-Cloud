#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 27, 2025"


from src.geometric_analysis import *
from src.gdmpidc import *
from src.gdmpidc_tools import *
import numpy as np
import time

start_time = time.time()

ϵ = 0.00001              # [m]
#parent_mass = 9.39e20    # [kg]
parent_mass = 100
N = cumulative_distribution

def subrange_count(L1, L2, breakup_type: str = "collision") -> int:
    """Number of fragments between lengths L1 and L2 where L2 > L1."""
    if breakup_type == "collision":
        number_of_fragments = N(L1, parent_mass) - N(L2, parent_mass)
    elif breakup_type == "explosion":
        number_of_fragments = N(L1, parent_mass, "explosion") - N(L2, parent_mass, "explosion")
    return int(number_of_fragments)

def point_count(Lc: float, size_step: float = ϵ, breakup_type: str = "collision") -> int:
    if breakup_type == "collision":
        number_of_fragments = N(Lc, parent_mass) - N(Lc + size_step, parent_mass)
    elif breakup_type == "explosion":
        number_of_fragments = N(Lc, parent_mass, "explosion") - N(Lc + size_step, parent_mass, "explosion")
    return int(number_of_fragments)

small_cloud, ΔLc_small = dict(), 0.0001
med_cloud, ΔLc_med     = dict(), 0.0005
large_cloud, ΔLc_large = dict(), 0.001

Lc = 0.001
while Lc < 0.08:
    small_cloud[Lc] = SubCloud(Lc, point_count(Lc))
    Lc += ΔLc_small
while 0.08 <= Lc <= 0.11:
    med_cloud[Lc] = SubCloud(Lc, point_count(Lc))
    Lc += ΔLc_med
while 0.11 < Lc <= 1.0:
    large_cloud[Lc] = SubCloud(Lc, point_count(Lc))
    Lc += ΔLc_large

cloud = {"small": small_cloud, "medium": med_cloud, "large": large_cloud}

relavent_points = []
for category in cloud:
    for Lc, sub_cloud in cloud[category].items():
        sub_points = sub_cloud.sample_positions()
        inside_points = [point for point in sub_points if np.linalg.norm(point) <= sub_cloud.radius]
        relavent_points.extend(inside_points)

print(len(relavent_points))
print(f"Processing time: {time.time() - start_time:.2f} seconds")
