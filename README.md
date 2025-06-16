
# Gaussian Distribution Model of Post-Impact Debris Cloud (GDMPIDC)

A comprehensive computational model for simulating and analyzing satellite debris clouds following fragmentation events in orbit. This project implements the theoretical framework described in the accompanying research papers for modeling the spatial density distribution of fragments using Gaussian functions.

## Overview

The GDMPIDC project provides tools for:
- **Debris Cloud Generation**: Creating realistic post-impact debris clouds with size-dependent fragment distributions
- **Geometric Analysis**: Computing collision probabilities and trajectory intersections
- **Monte Carlo Simulation**: Statistical analysis of impact probabilities at event zero
- **3D Visualization**: Plotting and analyzing fragment spatial distributions

## Key Features

### 🎯 **Physical Modeling**
- Implementation of NASA Standard Breakup Model for fragment size distributions
- Gaussian spatial density functions with size-dependent parameters
- Time-evolution modeling for short-term debris cloud expansion
- Area-to-mass ratio calculations for realistic fragment properties

### 📊 **Statistical Analysis**
- Monte Carlo estimation of collision probabilities
- Importance sampling for enhanced computational efficiency
- Adaptive sampling with sequential refinement
- Confidence interval calculations using Wilson score method

### 🔧 **Computational Tools**
- Efficient 3D geometric algorithms for line-sphere intersections
- Parametric line generation and point-to-line distance calculations
- Fragment position sampling from non-uniform density distributions
- Optimized numerical integration and statistical sampling

## Project Structure

```
src/
├── gdmpidc.py              # Core debris cloud classes and modeling
├── gdmpidc_tools.py        # NASA breakup model and utility functions
├── geometric_analysis.py   # 3D geometry and collision analysis
├── docs/
│   └── theoretical/
│       ├── gdmpidc.md      # Theoretical foundation document
│       └── cissdcm.md      # Computational implementation details
└── subscripts/
    ├── misc/
    │   └── simple_sim.py   # Example simulation script
    └── plotters/           # Visualization tools
```

## Quick Start

### Basic Debris Cloud Creation

```python
from src.gdmpidc import Cloud
from src.geometric_analysis import get_entry_exit, count_points_near_line

# Create a debris cloud from a 100 kg parent object with 10m radius
cloud = Cloud(parent_mass=100e3, parent_radius=10)

print(f"Total fragments: {len(cloud.all_points):,}")
print(f"Cloud radius: {cloud.radius:.2f} m")
```

### Collision Probability Analysis

```python
from src.geometric_analysis import importance_sample_entry_exit, line_parametric_3d

# Generate entry/exit points on cloud sphere
entry, exit = importance_sample_entry_exit(cloud.radius)

# Create parametric line through cloud
line = line_parametric_3d(entry, exit)

# Count fragments within hit distance
hit_distance = 1.0  # meters
hits = count_points_near_line(line, cloud.all_points, hit_distance)

print(f"Fragments within {hit_distance}m: {hits}")
print(f"Hit percentage: {hits/len(cloud.all_points)*100:.2f}%")
```

### Monte Carlo Impact Probability

```python
# Run the main analysis (as shown in main.py)
if __name__ == "__main__":
    from main import main
    main(parent_mass=10000, parent_radius=1000)
```

## Example Results

Based on a 10,000 kg parent object with 1,000 m radius:

```
Cloud created in 14.08 seconds
Total fragments: 804,105
Cloud radius: 1100.64 m 

Fragment distribution (inside cloud radius only):
  Small fragments (< 8 cm):  804,077 (99.997%)
  Medium fragments (8-11 cm): 28 (0.003%)
  Large fragments (> 11 cm):  0 (0.000%)

Monte Carlo Results:
  Impact Probability: 0.012941
  Hits: 1,307 out of 101,000 trials
  95% Confidence Interval: [0.012262, 0.013656]
```

## Mathematical Foundation

### Fragment Size Distribution

The model implements the NASA Standard Breakup Model cumulative distributions:

**Collision Events:**
```
N(L_c) = 0.1 × L_c^(-1.71) × M_parent^(0.75)
```

**Explosion Events:**
```  
N(L_c) = 6 × L_c^(-1.6)
```

### Spatial Density Function

Fragments follow a Gaussian radial distribution:

```
ρ(r, L_c, t) = ρ₀(L_c, t) × exp[-½((r - μR_c(t))/(σ(L_c,t)R_c(t)))²]
```

Where:
- `μ`: Peak density location (typically 0.6-0.64)
- `σ`: Size-dependent dispersion parameter  
- `R_c(t)`: Time-dependent cloud radius

### Impact Probability Calculation

Using non-homogeneous Poisson process theory:

```
P_impact = 1 - exp[-∫ Λ(s) ds]
```

Where `Λ(s)` is the collision rate along trajectory path `s`.

## Algorithm Features

### Size-Dependent Parameters

| Fragment Size | Peak Location (μ) | Dispersion (σ₀) | Alpha (α) |
|---------------|-------------------|-----------------|-----------|
| Small (< 8cm) | 0.64             | 0.12            | 0.5       |
| Medium (8-11cm)| 0.60             | 0.062           | 0.65      |
| Large (> 11cm)| 0.60             | 0.047           | 0.69      |

### Sampling Methods

1. **Uniform Sampling**: Standard random point generation on sphere surface
2. **Importance Sampling**: Biased toward high-density regions for efficiency
3. **Rejection Sampling**: For complex probability distributions
4. **Adaptive Sampling**: Sequential refinement with convergence criteria

## Performance Characteristics

- **Cloud Generation**: ~14 seconds for 800k+ fragments
- **Monte Carlo Analysis**: ~5 seconds for 100k trials with importance sampling
- **Memory Usage**: Optimized for large fragment populations
- **Scalability**: Handles parent masses from 1 kg to 100+ tons

## Research Applications

This implementation supports research in:
- **Space Situational Awareness**: Debris cloud tracking and prediction
- **Collision Risk Assessment**: Spacecraft-debris encounter probabilities  
- **Mission Planning**: Safe trajectory design through debris fields
- **Debris Mitigation**: Understanding fragmentation event consequences

## Dependencies

- `numpy`: Numerical computations and array operations
- `scipy`: Statistical distributions and numerical integration
- `matplotlib`: Visualization and plotting (for plotter scripts)

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd debris-cloud

# Install dependencies (handled automatically by Replit)
# Dependencies are managed through replit.nix and pyproject.toml
```

## Usage Examples

### Fragment Analysis

```python
# Analyze fragment distribution by size
for category in cloud.subclouds:
    for Lc, subcloud in cloud.subclouds[category].items():
        print(f"Size {Lc:.3f}m: {subcloud.nFrag} fragments")
        print(f"Subcloud radius: {subcloud.radius:.2f}m")
```

### Custom Trajectory Analysis

```python
# Define custom entry/exit points
entry_point = (0, 0, cloud.radius)
exit_point = (0, 0, -cloud.radius)
line = line_parametric_3d(entry_point, exit_point)

# Analyze hits with different distance thresholds
for distance in [0.5, 1.0, 2.0, 5.0]:
    hits = count_points_near_line(line, cloud.all_points, distance)
    print(f"Distance {distance}m: {hits} hits ({hits/len(cloud.all_points)*100:.3f}%)")
```

## Contributing

This project implements peer-reviewed research in orbital debris modeling. Contributions should maintain scientific accuracy and computational efficiency. Key areas for development:

- **Temporal Evolution**: Extended time-scale modeling beyond 120 seconds
- **Perturbation Forces**: Enhanced solar radiation pressure and drag models  
- **Visualization**: Advanced 3D plotting and animation capabilities
- **Optimization**: Performance improvements for larger-scale simulations

## References

1. Johnson, N. L., et al. (2001). NASA's New Breakup Model of EVOLVE 4.0. *Advances in Space Research*, **28**(9), 1377-1384.
2. Modjtahedzadeh, K. (2025). Gaussian Distribution Model of Post-Impact Debris Cloud. Boeing Intelligence & Analytics.
3. Modjtahedzadeh, K. (2025). Computational Implementation of Spherically Symmetric Debris Cloud Models. Boeing Intelligence & Analytics.

## License

This project implements scientific research methodologies for orbital debris analysis. Please cite the accompanying research papers when using this code in academic or commercial applications.

---

**Authors**: Kamyar Modjtahedzadeh, Claude 3.7 Sonnet, Grok 3, Claude 3.5 Sonnet V2  
**Institution**: Boeing Intelligence & Analytics  
**Date**: May 2025 - June 2025
