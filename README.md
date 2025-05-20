
# Gaussian Distribution Model of Post-Impact Debris Cloud

A computational model describing the evolution of satellite debris clouds following fragmentation events in orbit. The approach models the spatial density distribution of fragments using a Gaussian function that incorporates the physical characteristics of debris particles.

## Features

- Models post-impact debris (PID) cloud evolution
- Implements NASA Standard Breakup Model
- Calculates fragment distributions based on characteristic length
- Visualizes debris cloud spatial distribution
- Provides statistical analysis of fragment dispersion

## Requirements

- Python 3.11+
- NumPy ≥ 2.2.5
- SciPy ≥ 1.15.2
- Matplotlib (for visualization)
- Markdown

## Project Structure

```
├── src/
│   ├── docs/              # Documentation
│   │   ├── theoretical/   # Theoretical foundations
│   │   └── log_notes/     # Development notes
│   ├── gdmpidc.py        # Main implementation
│   ├── gdmpidc_tools.py  # Utility functions
│   ├── percent_inside.py # Statistical analysis
│   ├── plotter.py        # Visualization
│   └── plotter_inside.py # Inner fragment plotting
└── main.py               # Entry point
```

## Usage

1. Set up the fragment parameters:
```python
cloud = Cloud(characteristic_length=0.05, num_fragments=1000, breakup_type="collision")
```

2. Generate fragment positions:
```python
positions = cloud.sample_positions()
```

3. Visualize the debris cloud:
```python
# The visualization code will automatically execute when running the main script
```

## Model Parameters

- Small fragments (≤ 8 cm)
- Medium fragments (8-11 cm)
- Large fragments (> 11 cm)

Each size category has specific empirical parameters optimized for:
- Small fragments: ≳75% inside cloud radius
- Medium fragments: ≳90% inside cloud radius
- Large fragments: ≳99% inside cloud radius

## Documentation

Detailed documentation is available in the `src/docs/` directory:
- `gdmpidc.md`: Comprehensive model description
- `cissdcm.md`: Computational implementation details

## Author

Kamyar Modjtahedzadeh  
Boeing Intelligence & Analytics  
