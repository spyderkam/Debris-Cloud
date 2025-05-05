#!/usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh, Claude 3.7 Sonnet V2"
__date__ = "May 1, 2025"

"""
Implementation of Gaussian Distribution Model of Post-Impact Debris Cloud (GDMPIDC)
Origin is the center of the parent body at impact (t = 0)
"""

import numpy as np
from am_distribution_calculator import calculate_mass, cross_sectional_area, get_AM_value

# Constants
M_EARTH = 5.972e24       # Earth mass [kg]
G = 6.67430e-11          # Gravitational constant [m^3·kg^-1·s^-2]
P_SRP = 4.56             # Solar radiation pressure at Earth orbit [N/m^2]
C_R = 1.44               # Radiation pressure coefficient for debris

class Fragment:
    """PID (Post-Impact Debris) fragment model"""
    
    def __init__(self, characteristic_length: float, pos: np.ndarray, 
                 breakup_type: str = "collision", parent_vel: float = 0) -> None:
        """
        Initialize fragment with NASA Standard Breakup Model parameters
        
        Args:
            characteristic_length: Average of three orthogonal dimensions [m]
            pos: Initial position vector [x, y, z] in meters
            breakup_type: Either "collision" or "explosion"
            parent_vel: Parent object velocity at impact [m/s]
        """
        self.L_c = characteristic_length
        self.pos = np.array(pos)
        self.r = np.linalg.norm(self.pos)  # Radial distance
        
        # Physical properties
        self.mass = calculate_mass(self.L_c)
        self.area = cross_sectional_area(self.L_c)
        self.am_ratio = self.area / self.mass
        
        # Calculate ejection velocity (Eq 2.4 & 2.5)
        log_am = np.log10(get_AM_value(np.log10(self.L_c)))
        if breakup_type == "collision":
            ejection_vel = 0.9 * log_am + 2.9
        else:  # explosion
            ejection_vel = 0.2 * log_am + 1.85
            
        self.vel = parent_vel + ejection_vel
        
    def calculate_density(self, r: float, t: float, 
                         R_c: float, sigma: float, mu: float = 2/3) -> float:
        """
        Calculate Gaussian density distribution (Eq 4.12)
        
        Args:
            r: Radial distance from center [m]
            t: Time since impact [s]
            R_c: Cloud radius at time t [m]
            sigma: Standard deviation parameter
            mu: Peak density radial position as fraction of R_c
        
        Returns:
            Density at position r [kg/m^3]
        """
        # Time-dependent normalization (Eq 3.2)
        rho_0_t = self.mass * (R_c / r)**3
        
        # Gaussian distribution (Eq 4.12)
        exponent = -0.5 * ((r - mu*R_c)/(sigma*R_c))**2
        return rho_0_t * np.exp(exponent)
    
    def calculate_srp_acceleration(self, sun_direction: np.ndarray) -> np.ndarray:
        """
        Calculate solar radiation pressure acceleration (Eq 4.7)
        
        Args:
            sun_direction: Unit vector pointing from fragment to Sun
        
        Returns:
            Acceleration vector [m/s^2]
        """
        return -P_SRP * self.am_ratio * C_R * sun_direction
    
    def calculate_gravity_acceleration(self) -> np.ndarray:
        """Calculate gravitational acceleration (Eq 4.3)"""
        unit_r = self.pos / self.r
        return -G * M_EARTH / (self.r**2) * unit_r

if __name__ == "__main__":
    # Example usage
    fragment = Fragment(0.01, [0, 0, 0])
    
    # Calculate density at various radial positions
    R_c = 100  # Cloud radius [m]
    sigma = 0.3  # Standard deviation
    t = 60  # Time since impact [s]
    
    r_values = np.linspace(0, R_c, 100)
    densities = [fragment.calculate_density(r, t, R_c, sigma) for r in r_values]

