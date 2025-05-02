#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 1, 2025"

"""
    Let the origin be the center of the parent body at impact (𝑡 = 0).
"""

from am_distribution_calculator import calculate_mass, cross_sectional_area, get_AM_value
import numpy as np

#M_EARTH = 5.972e+24       # [kg]
#GRAV_CONST = 6.67430e-11  # [m^3·kg^-1·s^-2]
parent_vel = 0            # [m·s^-1]

class Fragment:
    """PID fragments."""

    def __init__(self, characteristic_length: float, pos: list, breakup_type="collision") -> None:
        self.size = characteristic_length                         # characteristic length
        self.pos = pos                                            # (x, y, z)
        self.radial = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)  # r = √(x^2 + y^2 + z^2)
        
        self.mass = calculate_mass(characteristic_length)
        self.area = cross_sectional_area(characteristic_length)
        
        if breakup_type == "collision":
            ejection_vel = 0.9*np.log10(get_AM_value(np.log10(characteristic_length))) + 2.9
        elif breakup_type == "explosion":
            ejection_vel = 0.2*np.log10(get_AM_value(np.log10(characteristic_length))) + 1.85
            
        self.vel = parent_vel + ejection_vel

def empirical_parameters(self):
    """
        Define empirical parameters for the fragment based on its characteristic length.

        Returns:
            tuple: (μ, ρ0, σ0, α, γ)
    """
    
    # Small fragments (≲ 1 cm)
    if self.size <= 0.01:
        ρ0 = 2.5           # Higher normalization constant for smaller fragments
        μ = 0.70           # Peak density slightly further out due to higher mobility
        γ = 0.005          # Faster evolution due to SRP and drag effects
        σ0 = 0.4           # Wider initial distribution due to higher ejection velocities
        α = 1.2            # Stronger size dependency
    # Medium fragments (1-10 cm)
    elif 0.01 < self.size <= 0.1:
        ρ0 = 2.0
        μ = 0.65
        γ = 0.003
        σ0 = 0.3
        α = 1.0
    # Large fragments (≳ 10 cm)
    else:
        ρ0 = 1.5           # Lower normalization constant for larger fragments
        μ = 0.60           # Peak closer to origin (less affected by dispersion)
        γ = 0.001          # Slower evolution due to smaller perturbation effects
        σ0 = 0.2           # Narrower distribution (less affected by ejection)
        α = 0.8            # Weaker size dependency
        
    return μ, ρ0, σ0, α, γ


if __name__ == "__main__":
    sample_fragment = Fragment(0.01, [0, 0, 0])
    #print(vars(sample_fragment))
    
    def dispersion(Lc, t):
        return σ0*Lc**(-α) + γ*t*Lc**(-α)
    