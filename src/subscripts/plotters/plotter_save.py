__author__ = "Claude 4.0 Sonnet"
__date__ = "June 2, 2025"

# Interactive 3D scatter plot using Plotly

import plotly.graph_objects as go
import plotly.offline as pyo
import numpy as np

## Assuming vec_r is a list of lists like [[x0, y0, z0], [x1, y1, z1], ...]
vec_r = np.array(vec_r)  # Convert to numpy array for easier handling

## Extract x, y, z coordinates
x = vec_r[:, 0]
y = vec_r[:, 1]
z = vec_r[:, 2]

# Create scatter plot trace
scatter_trace = go.Scatter3d(
    x=x, y=y, z=z,
    mode='markers',
    marker=dict(
        size=3,
        color='blue',
        opacity=0.8
    ),
    name='Debris Fragments'
)

# Generate spherical shell coordinates for wireframe
u = np.linspace(0, 2 * np.pi, 20)  # Azimuthal angle
v = np.linspace(0, np.pi, 20)      # Polar angle
x_sphere = cloud.radius * np.outer(np.cos(u), np.sin(v))
y_sphere = cloud.radius * np.outer(np.sin(u), np.sin(v))
z_sphere = cloud.radius * np.outer(np.ones(np.size(u)), np.cos(v))

# Create wireframe surface for the spherical shell
wireframe_trace = go.Surface(
    x=x_sphere, y=y_sphere, z=z_sphere,
    opacity=0.3,
    colorscale=[[0, 'red'], [1, 'red']],
    showscale=False,
    name='Cloud Boundary'
)

# Create the figure
fig = go.Figure(data=[scatter_trace, wireframe_trace])

# Update layout
fig.update_layout(
    title=f'Interactive Debris Cloud Visualization<br>{len(vec_r)} fragments within {cloud.radius:.1f}m radius',
    scene=dict(
        xaxis_title='x [m]',
        yaxis_title='y [m]',
        zaxis_title='z [m]',
        aspectmode='cube'  # Equal aspect ratio
    ),
    width=800,
    height=600
)

# Save as interactive HTML file to 3d_htmls directory
import os
output_dir = os.path.join(os.path.dirname(__file__), "subscripts", "3d_htmls")
#output_dir = os.path.join(os.path.dirname(__file__), "3d_htmls")
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "interactive_cloud.html")
pyo.plot(fig, filename=output_path, auto_open=False)
print(f"Interactive plot saved as '{output_path}' - open this file in your browser!")
