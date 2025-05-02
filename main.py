#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 1, 2025"

"""
    Let the origin be the center of the parent body at impact (𝑡 = 0).
"""

from am_distribution_calculator import calculate_mass, cross_sectional_area, get_AM_value
import numpy as np

M_EARTH = 5.972e+24       # [kg]
GRAV_CONST = 6.67430e-11  # [m^3·kg^-1·s^-2]
parent_vel = 0            # [m·s^-1]
# ejec_vel = ???

class Fragment:
    """PID fragments"""

    def __init__(self, characteristic_length: float, pos: list, breakup_type="collision") -> None:
        self.characteristic_length = characteristic_length        # characteristic length
        self.pos = pos                                            # (x, y, z)
        self.radial = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)  # r = √(x^2 + y^2 + z^2)
        
        self.mass = calculate_mass(characteristic_length)
        self.area = cross_sectional_area(characteristic_length)
        
        if breakup_type == "collision":
            ejection_vel = 0.9*np.log10(get_AM_value(np.log10(characteristic_length))) + 2.9
        elif breakup_type == "explosion":
            ejection_vel = 0.2*np.log10(get_AM_value(np.log10(characteristic_length))) + 1.85
            
        self.vel = parent_vel + ejection_vel


if __name__ == "__main__":
    sample_fragment = Fragment(0.01, [0, 0, 0])
    #(vars(sample_fragment))

    μ = 2/3     # varies with fragment area-to-mass ratio
    