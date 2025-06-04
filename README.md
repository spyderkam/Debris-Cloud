
# Gaussian Distribution Model of Post-Impact Debris Cloud (GDMPIDC)

A Python-based simulation framework for modeling debris clouds generated from space object breakup events using Gaussian distribution models and NASA Standard Breakup Model parameters.

## Overview

This project implements a sophisticated debris cloud modeling system that simulates the spatial distribution and temporal evolution of fragments following collision or explosion events in space. The model uses empirically-derived parameters and follows established NASA breakup models to provide realistic debris cloud simulations.

## Core Components

### 1. Fragment Modeling (`src/gdmpidc.py`)

The core module defines three main classes:

#### `Fragment`
Represents individual debris fragments with properties including:
- **Size**: Characteristic length (meters)
- **Position**: 3D coordinates (x, y, z)
- **Mass**: Calculated using NASA Standard Breakup Model
- **Velocity**: Ejection velocity based on breakup type (collision/explosion)
- **Physical Properties**: Cross-sectional area, empirical parameters

#### `SubCloud`
Models fragments of a specific characteristic length as a sub-population within the larger debris cloud:
- **Gaussian Distribution**: Fragments follow a spherically symmetric Gaussian density distribution
- **Temporal Evolution**: Cloud radius expands over time with configurable expansion velocity
- **Density Calculation**: Provides spatial density at any point and time
- **Position Sampling**: Generates realistic fragment positions using Monte Carlo methods

#### `Cloud`
Manages the complete debris cloud containing multiple `SubCloud` instances:
- **Size Resolution**: Different resolution levels for small (< 8cm), medium (8-11cm), and large (> 11cm) fragments
- **Adaptive Sampling**: Variable resolution based on fragment size
- **Spatial Filtering**: Automatically filters fragments within the cloud boundary

### 2. Supporting Tools (`src/gdmpidc_tools.py`)

Provides essential calculation functions:

#### NASA Standard Breakup Model Implementation
- **Area-to-Mass Ratio**: `get_AM_value()` - Samples A/M ratios for fragments
- **Size Distributions**: Handles both collision and explosion breakup types
- **Fragment Counting**: `point_count()`, `subrange_count()` for fragment population calculations

#### Physical Calculations
- **Mass Calculation**: `calculate_mass()` - Determines fragment mass from characteristic length
- **Cross-sectional Area**: `cross_sectional_area()` - Calculates fragment area
- **Expansion Velocity**: `expansion_velocity()` - Computes average cloud expansion rate

#### Empirical Parameters
- **Fragment Dispersion**: `empirical_parameters()` - Returns optimized dispersion parameters
- **Packing Density**: `packing_density()` - Initial spatial packing based on fragment size

### 3. Geometric Analysis (`src/geometric_analysis.py`)

Implements geometric operations for trajectory analysis:

#### Point Generation
- **`get_entry_exit()`**: Generates random entry/exit points on sphere surface
- **Diameter Control**: Optional parameter to prevent/allow diametrically opposite points

#### Line Operations  
- **`line_parametric_3d()`**: Creates parametric line equations in 3D space
- **`count_points_near_line()`**: Counts fragments within specified distance of a trajectory line

## Usage Examples

### Basic Cloud Creation

```python
from src.gdmpidc import Cloud
from src.gdmpidc_tools import *

# Create debris cloud from 100,000 kg parent object with 10m radius
parent_mass = 100e3  # kg
parent_radius = 10   # meters
cloud = Cloud(parent_mass, parent_radius)

print(f"Total fragments: {len(cloud.all_points)}")
print(f"Cloud radius: {cloud.radius:.2f} meters")
```

### Trajectory Analysis

```python
from src.geometric_analysis import get_entry_exit, line_parametric_3d, count_points_near_line

# Generate random trajectory through cloud
p1, p2 = get_entry_exit(cloud.radius, diameter=False)
trajectory = line_parametric_3d(p1, p2)

# Count fragments within 1 meter of trajectory
hit_distance = 1.0
hits = count_points_near_line(trajectory, cloud.all_points, hit_distance)
hit_percentage = (hits / len(cloud.all_points)) * 100

print(f"Fragments within {hit_distance}m of trajectory: {hits} ({hit_percentage:.2f}%)")
```

### Fragment Analysis by Size

```python
# Analyze specific fragment size
from src.gdmpidc import SubCloud

# Create sub-cloud for 5cm fragments
subcloud = SubCloud(
    characteristic_length=0.05,  # 5cm fragments
    parent_rad=10,              # parent radius
    num_fragments=1000,         # number of fragments
    breakup_type="collision"    # or "explosion"
)

# Sample fragment positions
positions = subcloud.sample_positions(100)
print(f"Generated {len(positions)} fragment positions")

# Calculate density at specific location
density = subcloud.density([0, 0, 0], t=0)  # density at origin at t=0
print(f"Density at origin: {density:.6f}")
```

## Simple Simulation Example

The `simple_sim.py` demonstrates a complete workflow:

1. **Cloud Generation**: Creates debris cloud from parent object breakup
2. **Trajectory Definition**: Generates random entry/exit points through cloud
3. **Impact Analysis**: Counts fragments intersecting with trajectory path
4. **Performance Metrics**: Reports computation times and hit statistics

```bash
python simple_sim.py
```

Sample output:
```
Computational PID cloud creation time: 2.34 seconds

• A mass of 93800 kg has been broken up into 15432 fragments between 0.001 and 1.0 meters
• The radius of the debris cloud at collision is 18.45 meters
• The number of fragments within 1.0 meters of trajectory: 127 (0.82% of total)
• Computation time for point count: 0.03 seconds
```

## Key Features

### Realistic Physics
- **NASA Standard Breakup Model**: Implements official breakup model parameters
- **Empirical Validation**: Parameters optimized for realistic fragment distributions
- **Temporal Evolution**: Models cloud expansion and dispersion over time

### Computational Efficiency
- **Adaptive Resolution**: Variable sampling resolution based on fragment size
- **Monte Carlo Methods**: Efficient position sampling using Gaussian distributions  
- **Linear Complexity**: O(n) algorithms for geometric calculations

### Flexibility
- **Multiple Breakup Types**: Supports both collision and explosion scenarios
- **Configurable Parameters**: Adjustable cloud properties and simulation parameters
- **Modular Design**: Independent components for easy extension and modification

## Scientific Applications

This framework is designed for:
- **Space Debris Risk Assessment**: Analyzing collision probabilities and debris evolution
- **Mission Planning**: Evaluating spacecraft trajectory risks through debris fields
- **Debris Mitigation Studies**: Understanding fragment dispersion patterns
- **Monte Carlo Simulations**: Statistical analysis of debris cloud behavior

## Installation & Dependencies

The project uses standard scientific Python libraries:
- `numpy` - Numerical computations and array operations
- `scipy` - Scientific computing and integration functions

```bash
# Dependencies are automatically managed by Replit
# No manual installation required
```

## Performance Considerations

- **Cloud Creation**: Typically 1-5 seconds for 10,000-50,000 fragments
- **Geometric Analysis**: Linear time complexity O(n) for point-line distance calculations
- **Memory Usage**: Efficient storage of fragment positions and properties
- **Scalability**: Handles fragment populations from hundreds to hundreds of thousands

## File Structure

```
src/
├── gdmpidc.py           # Core debris cloud classes
├── gdmpidc_tools.py     # NASA breakup model tools
├── geometric_analysis.py # Trajectory analysis functions
└── docs/                # Technical documentation
    ├── theoretical/     # Mathematical background
    └── log_notes/       # Development notes

simple_sim.py           # Example simulation
main.py                 # Main entry point
```

## Authors

- **Kamyar Modjtahedzadeh** - Primary developer and researcher
- **Claude AI** - Development assistance and optimization

## License

This project is developed for scientific and educational purposes. Please refer to your organization's licensing requirements for commercial use.

## Contributing

This is a research-focused project. For questions or collaboration opportunities, please contact the primary author through the repository.
