import numpy as np
import math

def sample_AM_distribution(log_Lc, fragment_type="upper_stage"):
    """
    Sample from the NASA breakup model's area-to-mass ratio distribution.
    
    Args:
        log_Lc: log10 of the characteristic length in meters
        fragment_type: "upper_stage" or "spacecraft"
    
    Returns:
        A sampled A/M value in m²/kg
    """
    
    # Determine which distribution to use based on size
    if log_Lc > math.log10(0.11):  # Larger than 11 cm
        return sample_large_AM_distribution(log_Lc, fragment_type)
    elif log_Lc < math.log10(0.08):  # Smaller than 8 cm
        return sample_small_AM_distribution(log_Lc)
    else:  # Between 8 and 11 cm - transition region
        # Linearly interpolate between small and large fragment distributions
        weight = (log_Lc - math.log10(0.08)) / (math.log10(0.11) - math.log10(0.08))
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
    if np.random.random() < alpha:
        # Sample from first distribution
        log_AM = np.random.normal(mu1, sigma1)
    else:
        # Sample from second distribution
        log_AM = np.random.normal(mu2, sigma2)
    
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
    log_AM = np.random.normal(mu, sigma)
    
    # Convert from log to linear
    return 10**log_AM



def calculate_mass(Lc, fragment_type="upper_stage"):
    """Calculate fragment mass based on characteristic length using NASA model."""
    
    # Step 1: Calculate cross-sectional area
    if Lc < 0.00167:
        area = 0.540424 * (Lc ** 2)
    else:
        area = 0.556945 * (Lc ** 2.0047077)
    
    # Step 2: Sample A/M ratio from the appropriate distribution
    log_Lc = math.log10(Lc)
    AM_ratio = sample_AM_distribution(log_Lc, fragment_type)
    
    # Step 3: Calculate mass
    mass = area / AM_ratio
    
    return mass