#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 1, 2025"

from am_distribution_calculator import get_AM_value, calculate_mass
import numpy as np

M_EARTH = 5.972e+24       # [kg]
GRAV_CONST = 6.67430e-11  # [m^3·kg^-1·s^-2]
parent_vel = 0            # [m·s^-1]
# ejec_vel = ???

class Fragment:
    """PID fragments"""

    def __init__(self, size: float, pos: list, vel: float, breakup_type="collision") -> None:
        self.size = size                                          # characteristic length
        self.pos = pos                                            # (x, y, z)
        self.radial = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)  # r = √(x^2 + y^2 + z^2)
        self.velocity = vel                                       # v_fragment = v_parent + v_ejection
        # Should velocity be a vector
        