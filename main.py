#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 27, 2025"

from src.geometric_analysis import *
from src.gdmpidc import *
from src.gdmpidc_tools import *

ϵ = 0.00001              # [m]
parent_mass = 9.39e20    # [kg]

def subrange_count(L1, L2, breakup_type: str = "collision") -> int:
    """Number of fragments between lengths L1 and L2 where L2 > L1."""
    if breakup_type == "collision":
        return cumulative_distribution(L1, parent_mass) - cumulative_distribution(L2, parent_mass)
    elif breakup_type == "explosion":
        return cumulative_distribution(L1, parent_mass, "explosion") - cumulative_distribution(L2, parent_mass, "explosion")

def point_count(Lc: float, size_step: float = ϵ, breakup_type: str = "collision") -> int:
    if breakup_type == "collision":
        return cumulative_distribution(Lc, parent_mass) - cumulative_distribution(Lc + size_step, parent_mass)
    elif breakup_type == "explosion":
        return cumulative_distribution(Lc, parent_mass, "explosion") - cumulative_distribution(Lc + size_step, parent_mass, "explosion")

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
