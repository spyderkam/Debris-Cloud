#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 1, 2025 - May 5, 2025"

"""
    Let the origin be the center of the parent body at impact (𝑡 = 0).
"""

from gdmpidc_tools import calculate_mass, cross_sectional_area, empirical_parameters, expansion_velocity, get_AM_value, packing_density
import numpy as np

# M_EARTH = 5.972e+24       # [kg]
# GRAV_CONST = 6.67430e-11  # [m^3·kg^-1·s^-2]
parent_vel = 0              # [m·s^-1]
parent_rad = 4.1439e+11     # [m]


class Fragment:
    """PID fragments."""

    def __init__(self, characteristic_length: float, pos: list, creation_type="collision") -> None:
        self.size = characteristic_length                         # characteristic length
        self.pos = pos                                            # (x, y, z)
        self.radial = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)  # r = √(x^2 + y^2 + z^2)

        self.mass = calculate_mass(characteristic_length)
        self.area = cross_sectional_area(characteristic_length)

        if creation_type == "collision":
            ejection_vel = 0.9*np.log10(get_AM_value(np.log10(characteristic_length))) + 2.9
        elif creation_type == "explosion":
            ejection_vel = 0.2*np.log10(get_AM_value(np.log10(characteristic_length))) + 1.85
        self.vel = parent_vel + ejection_vel

        # Get empirical parameters
        self.μ, self.ρ0, self.σ0, self.α, self.γ = empirical_parameters(self.size)

    def dispersion(self, t: float) -> float:
        """
            Spacial Distribution of Fragment; Equation (4.11) of gdmpidc.md. 

            Args:
                t (float): Time since impact [s]

            Returns:
                float: Spacial distribution of fragment with specific charachteristic length at time t.
        """
        return self.σ0*self.size**(-self.α) + self.γ*t*self.size**(-self.α)


class Cloud:
    def __init__(self, characteristic_length, num_fragments, breakup_type="collision") -> None:
        self.radius = parent_rad*packing_density(characteristic_length)**(-1/3)
        #self.radius = initRadius
        self.nFrag = num_fragments
        self.fragSize = characteristic_length
        self.breakup_type = breakup_type

    def cloud_radius(self, t: float) -> float:
        """
            Cloud radius at time t; Equation (3.1) of gdmpidc.md.
            
            Args:
                t (float): Time since impact [s]
            
            Returns:
                float: Cloud radius at time t [m]
        """
        return self.radius + t*expansion_velocity(parent_mass=1000, L_min=0.001, L_max=15.0)

    def sample_fragment(self, pos) -> Fragment:
        return Fragment(self.fragSize, pos, creation_type=self.breakup_type)

    #def density(self) -> float:
    #    ρ0*np.exp(-0.5 * ((r - μ*self.radius)/(σ*self.radius))**2)


if __name__ == "__main__":
    #sample_fragment = Fragment(0.01, [0, 0, 0])
    #print(type(sample_fragment))

    cloud = Cloud(10,12)
    
