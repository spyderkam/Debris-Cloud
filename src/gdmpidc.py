#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 1, 2025 - May 5, 2025"

"""
    Ceres asteroid; let the origin be the center of the parent body at impact (𝑡 = 0).
"""

from gdmpidc_tools import calculate_mass, cross_sectional_area, empirical_parameters, expansion_velocity, get_AM_value, packing_density
import numpy as np

# M_EARTH = 5.972e+24       # [kg]
# GRAV_CONST = 6.67430e-11  # [m^3·kg^-1·s^-2]
parent_vel = 0              # [m·s^-1]
parent_rad = 4.1439e+11     # [m]


class Fragment:
    """PID fragments."""

    def __init__(self, characteristic_length: float, pos: list, creation_type: str = "collision") -> None:
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

    '''  Not sure if I need this method, since it is already defined for Cloud..?
    def dispersion(self, t: float) -> float:
        """
            Spacial Distribution of Fragment; Equation (4.11) of gdmpidc.md. 

            Args:
                𝑡 (float): Time since impact [s]

            Returns:
                float: Spacial distribution of fragment with specific charachteristic length at time 𝑡.
        """
        return self.σ0*self.size**(-self.α) + self.γ*t*self.size**(-self.α)
    '''


class Cloud:
    def __init__(self, characteristic_length, num_fragments, breakup_type: str = "collision") -> None:
        self.nFrag = num_fragments
        self.fragSize = characteristic_length
        self.breakup_type = breakup_type

        self.radius = parent_rad * packing_density(characteristic_length)**(-1/3)  # Eqn. (1.1) of gdmpidc.md
        self.fragments = self._initialize_fragments()                              # initialize fragments with sampled positions

    def _initialize_fragments(self) -> list:
        """Initialize nFrag fragments with positions sampled from the cloud's density distribution."""
        positions = self.sample_positions(self.nFrag)
        return [Fragment(self.fragSize, pos, creation_type=self.breakup_type) for pos in positions]

    def sample_positions(self, n: int) -> list:
        """
        Sample n positions from the cloud's Gaussian density distribution.

        Args:
            n (int): Number of positions to sample

        Returns:
            list: List of [x, y, z] coordinates
        """
        
        # Get empirical parameters for the Gaussian distribution
        μ, _, _, _, _ = empirical_parameters(self.fragSize)
        # Use initial spatial dispersion (t=0)
        σ = self.spatial_dispersion(t=0)  # σ = σ0 * L^(-α)
        Rc = self.radius

        # The density is proportional to exp(-0.5 * ((r - μ*R)/(σ*R))^2)
        # In spherical coordinates, the radial probability includes a r^2 factor
        # We sample r from a distribution proportional to r^2 * exp(-0.5 * ((r - μ*R)/(σ*R))^2)

        # Define the mean and standard deviation in radial distance
        mean_r = μ * Rc
        std_r = σ * Rc

        # To sample from the non-standard distribution, use rejection sampling or approximate with a Gaussian
        # For simplicity, we'll approximate by sampling from a Gaussian and adjust for spherical symmetry
        # Note: This is an approximation; for high precision, implement rejection sampling

        # Sample radial distances (approximate Gaussian, clipped to avoid negative r)
        r = np.random.normal(loc=mean_r, scale=std_r, size=n)
        r = np.maximum(r, 0)  # Ensure non-negative radial distances

        # Sample angular coordinates for isotropic distribution
        theta = np.arccos(2 * np.random.uniform(0, 1, n) - 1)  # Uniform in cos(θ)
        phi = np.random.uniform(0, 2 * np.pi, n)               # Uniform in φ

        # Convert to Cartesian coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Return list of [x, y, z] coordinates
        return [[x[i], y[i], z[i]] for i in range(n)]




    def cloud_radius(self, t: float, inplace: bool = False) -> float:
        """
            Cloud radius at time t; Equation (3.1) of gdmpidc.md.
            
            Args:
                t (float): Time since impact [s]
                cloud_radius (float, optional): If provided, updates self.radius with this value [m]
            
            Returns:
                float: Cloud radius at time t [m]
        """
        updated_radius = self.radius + t*expansion_velocity(parent_mass=1000, L_min=0.001, L_max=15.0)
        if inplace:
            self.radius = updated_radius
        return updated_radius
    
    def density(self, pos: list, t: float) -> float:
        """Density of the cloud at time 𝑡 and radial position 𝑟."""
        r = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
        μ, ρ0, _, _, _ = empirical_parameters(self.fragSize)
        σ = self.spatial_dispersion(t)
        Rc = self.cloud_radius(t)
        return ρ0*np.exp(-0.5 * ((r - μ*Rc)/(σ*Rc))**2)

    def spatial_dispersion(self, t: float) -> float:
        """Standard deviation of the fragments of the cloud at time 𝑡."""
        _, _, σ0, α, γ = empirical_parameters(self.fragSize)
        return self.fragSize**(-α) * (σ0 + γ*t)


if __name__ == "__main__":
    cloud = Cloud(10.0, 10)  # Initial radius

    # Calculate radius at t=2 without updating self.radius
    result = cloud.cloud_radius(t=2.0)
    print(result)  # Uses self.radius = 10.0
    print(cloud.radius)  # Still 10.0

    # Calculate radius at t=2 and update self.radius
    result = cloud.cloud_radius(t=2.0, inplace=True)
    print(result)  # Uses updated self.radius = 20.0
    print(cloud.radius)  # Now 20.0
    
