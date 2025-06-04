#!usr/bin/env python3

__author__ = "Claude 3.7 Sonnet, Kamyar Modjtahedzadeh, Claude 3.5 Sonnet V2"
__date__ = "May 1, 2025 - June 2, 2025"

import warnings     # Must be imported before scipy and numpy.
from numpy import log10, random
from scipy import integrate

warnings.filterwarnings("ignore", category=integrate.IntegrationWarning)


# NASA Standard Breakup Model
def get_AM_value(log_Lc, fragment_type="upper_stage"):
    """
    Sample implementation from the NASA breakup model's area-to-mass ratio distribution.
    
    Args:
        log_Lc: log10 of the characteristic length in meters
        fragment_type: "upper_stage" or "spacecraft"
    
    Returns:
        A sampled A/M value in m²/kg
    """

    # Determine which distribution to use based on size
    if log_Lc > log10(0.11):  # Larger than 11 cm
        return sample_large_AM_distribution(log_Lc, fragment_type)
    elif log_Lc < log10(0.08):  # Smaller than 8 cm
        return sample_small_AM_distribution(log_Lc)
    else:  # Between 8 and 11 cm - transition region
        # Linearly interpolate between small and large fragment distributions
        weight = (log_Lc - log10(0.08)) / (log10(0.11) - log10(0.08))
        large_sample = sample_large_AM_distribution(log_Lc, fragment_type)
        small_sample = sample_small_AM_distribution(log_Lc)
        return small_sample * (1 - weight) + large_sample * weight


def sample_large_AM_distribution(log_Lc, fragment_type):
    """Sample from the distribution for fragments larger than 11 cm."""
    if fragment_type == "upper_stage":
        # Upper stage/Rocket body parameters
        # Weight parameter alpha
        if log_Lc <= -1.4:
            alpha = 1.0
        elif log_Lc < 0:
            alpha = 1.0 - 0.3571 * (log_Lc + 1.4)
        else:
            alpha = 0.5

        # First distribution mean
        if log_Lc <= -0.5:
            mu1 = -0.45
        elif log_Lc < 0:
            mu1 = -0.45 - 0.9 * (log_Lc + 0.5)
        else:
            mu1 = -0.9

        # First distribution standard deviation
        sigma1 = 0.55

        # Second distribution mean
        mu2 = -0.9

        # Second distribution standard deviation
        if log_Lc <= -1.0:
            sigma2 = 0.28
        elif log_Lc < 0.1:
            sigma2 = 0.28 - 0.1636 * (log_Lc + 1)
        else:
            sigma2 = 0.1
    else:
        # Spacecraft parameters
        # Weight parameter alpha
        if log_Lc <= -1.95:
            alpha = 0.0
        elif log_Lc < 0.55:
            alpha = 0.3 + 0.4 * (log_Lc + 1.2)
        else:
            alpha = 1.0

        # First distribution mean
        if log_Lc <= -1.1:
            mu1 = -0.6
        elif log_Lc < 0:
            mu1 = -0.6 - 0.318 * (log_Lc + 1.1)
        else:
            mu1 = -0.95

        # First distribution standard deviation
        if log_Lc <= -1.3:
            sigma1 = 0.1
        elif log_Lc < -0.3:
            sigma1 = 0.1 + 0.2 * (log_Lc + 1.3)
        else:
            sigma1 = 0.3

        # Second distribution mean
        if log_Lc <= -0.7:
            mu2 = -1.2
        elif log_Lc < -0.1:
            mu2 = -1.2 - 1.333 * (log_Lc + 0.7)
        else:
            mu2 = -2.0

        # Second distribution standard deviation
        if log_Lc <= -0.5:
            sigma2 = 0.5
        elif log_Lc < -0.3:
            sigma2 = 0.5 - (log_Lc + 0.5)
        else:
            sigma2 = 0.3

    # Sample from the bimodal distribution
    if random.random() < alpha:
        # Sample from first distribution
        log_AM = random.normal(mu1, sigma1)
    else:
        # Sample from second distribution
        log_AM = random.normal(mu2, sigma2)

    # Convert from log to linear
    return 10**log_AM


def sample_small_AM_distribution(log_Lc):
    """Sample from the distribution for fragments smaller than 8 cm."""
    # SOC = Standard Objects and Components (applies to both spacecraft and upper stages)

    # Mean
    if log_Lc <= -1.75:
        mu = -0.3
    elif log_Lc < -1.25:
        mu = -0.3 - 1.4 * (log_Lc + 1.75)
    else:
        mu = -1.0

    # Standard deviation
    if log_Lc <= -3.5:
        sigma = 0.2
    else:
        sigma = 0.2 + 0.1333 * (log_Lc + 3.5)

    # Sample from normal distribution
    log_AM = random.normal(mu, sigma)

    # Convert from log to linear
    return 10**log_AM


def cross_sectional_area(Lc):
    """Calculate cross-sectional area based on characteristic length."""
    if Lc < 0.00167:
        area = 0.540424 * (Lc**2)
    else:
        area = 0.556945 * (Lc**2.0047077)
    return area


def calculate_mass(Lc, fragment_type="upper_stage"):
    """Calculate fragment mass based on characteristic length using NASA model."""
    area = cross_sectional_area(Lc)
    log_Lc = log10(Lc)
    AM_ratio = get_AM_value(log_Lc, fragment_type)
    mass = area / AM_ratio
    return mass


def expansion_velocity(parent_mass=1000, L_min=0.001, L_max=1.0, breakup_type="collision"):
    """
    Calculate the average expansion velocity of the debris cloud based on Eq. (3.3) of gdmpidc.md.
    
    Args:
        parent_mass (float): Mass of the parent object in kg
        L_min (float): Minimum characteristic length in meters
        L_max (float): Maximum characteristic length in meters
        breakup_type (str): "collision" or "explosion"
        
    Returns:
        float: Average expansion velocity in m/s
    """
    # Define differential size distribution function n(L_c) based on Eq. (2.8), (2.9), (2.10)
    def differential_size_distribution(L_c):
        if breakup_type == "explosion":
            # n(L_c) = -d/dL_c[N(L_c)] = -d/dL_c[6*L_c^(-1.6)]
            return 6 * 1.6 * L_c**(-2.6)
        else:  # collision
            # n(L_c) = -d/dL_c[N(L_c)] = -d/dL_c[0.1*L_c^(-1.71)*M_parent^0.75]
            return 0.1 * 1.71 * L_c**(-2.71) * parent_mass**0.75

    # Define average velocity for fragments of size L_c based on Eq. (2.4), (2.5)
    def avg_velocity(L_c):
        # Calculate A/M ratio for this characteristic length
        log_Lc = log10(L_c)
        AM_ratio = get_AM_value(log_Lc)

        if breakup_type == "explosion":
            return 0.2 * log10(AM_ratio) + 1.85
        else:  # collision
            return 0.9 * log10(AM_ratio) + 2.9

    # Numerator: ∫ v̄(L_c) × n(L_c) dL_c
    def integrand_numerator(L_c):
        return avg_velocity(L_c) * differential_size_distribution(L_c)

    # Denominator: ∫ n(L_c) dL_c
    def integrand_denominator(L_c):
        return differential_size_distribution(L_c)

    # Compute the integrals
    numerator, _ = integrate.quad(integrand_numerator, L_min, L_max)
    denominator, _ = integrate.quad(integrand_denominator, L_min, L_max)

    # Return the average expansion velocity
    if denominator != 0:
        return numerator / denominator
    else:
        return 0.0


def empirical_parameters(Lc):
    """
        Define empirical parameters for the fragment based on its characteristic length.
        These parameters have been optimized to achieve the following coverage:
          - Small fragments (𝐿c ≈ 5 cm): ≳75% inside cloud radius
          - Medium fragments (𝐿c ≈ 9 cm): ≳90% inside cloud radius
          - Large fragments (𝐿c ≈ 15 cm): ≳99% inside cloud radius
        
        Returns:
            tuple: (μ, ρ0, σ0, α, γ)
    """

    # Small fragments
    if Lc < 0.08:
        μ, ρ0, σ0, α, γ = 0.64, 3.0, 0.12, 0.5, 0.005
    # Medium fragments
    elif 0.08 <= Lc <= 0.11:
        μ, ρ0, σ0, α, γ = 0.60, 4.0, 0.062, 0.65, 0.003
    # Large fragments
    else:
        μ, ρ0, σ0, α, γ = 0.60, 3.0, 0.047, 0.69, 0.001

    return μ, ρ0, σ0, α, γ


def packing_density(Lc: float) -> float:
    """
    Initial packing density of a cloud subsphere with fragment size 𝐿c.
    Please note that these parameters are Claude 3.5 Sonnet V2's best guesses.
    """
    if Lc < 0.08:             # Small fragments
        return 0.55           ## More dispersed due to higher ejection velocities
    elif 0.08 <= Lc <= 0.11:  # Medium fragments
        return 0.65           ## Moderate packing
    else:                     # Large fragments
        return 0.75           ## More tightly packed due to lower dispersion


def cumulative_distribution(Lc: float, M_parent: float = None, breakup_type: str ="collision") -> float:
    """Cumulative Distibution function (𝑁); represents the number of particles greater than 𝐿c."""
    if breakup_type == "collision":
        return 0.1*Lc**(-1.71)*M_parent**(0.75)
    elif breakup_type == "explosion":
        return 6*Lc**(-1.71)


def subrange_count(L1: float, L2: float, parent_mass: float, breakup_type: str = "collision") -> int:
    """Number of fragments between lengths 𝐿1 and 𝐿2 where 𝐿2 > 𝐿1."""
    if breakup_type == "collision":
        number_of_fragments = cumulative_distribution(L1, parent_mass) - cumulative_distribution(L2, parent_mass)
    elif breakup_type == "explosion":
        number_of_fragments = cumulative_distribution(L1, parent_mass, "explosion") - cumulative_distribution(L2, parent_mass, "explosion")
    return int(number_of_fragments)


def point_count(Lc: float, parent_mass: float, epsilon: float = 0.00001, breakup_type: str = "collision") -> int:
    if breakup_type == "collision":
        number_of_fragments = cumulative_distribution(Lc, parent_mass) - cumulative_distribution(Lc + epsilon, parent_mass)
    elif breakup_type == "explosion":
        number_of_fragments = cumulative_distribution(Lc, parent_mass, "explosion") - cumulative_distribution(Lc + epsilon, parent_mass, "explosion")
    return int(number_of_fragments)
