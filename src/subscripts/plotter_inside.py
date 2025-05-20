
__date__ = "May 18, 2025"

# 3D scatter plot for fragments inside cloud radius only
import matplotlib.pyplot as plt
import numpy as np

# Filter points inside cloud radius
inside_points = [point for point in vec_r if np.linalg.norm(point) <= cloud.radius]
inside_points = np.array(inside_points)

# Extract x, y, z coordinates of inside points
x = inside_points[:, 0]
y = inside_points[:, 1]
z = inside_points[:, 2]

# Create a new figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create scatter plot for inside points
ax.scatter(x, y, z, c='blue', marker='o', s=20)

# Set labels for axes
ax.set_xlabel(r'$x$', fontsize=20)
ax.set_ylabel(r'$y$', fontsize=20)
ax.set_zlabel(r'$z$', fontsize=20)

# Turn off axis numbers (tick labels)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# Generate spherical shell coordinates
u = np.linspace(0, 2 * np.pi, 20)  # Azimuthal angle
v = np.linspace(0, np.pi, 20)      # Polar angle
x_sphere = cloud.radius * np.outer(np.cos(u), np.sin(v))
y_sphere = cloud.radius * np.outer(np.sin(u), np.sin(v))
z_sphere = cloud.radius * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the spherical shell as a wireframe
ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='r', alpha=0.5, linewidth=1.5)

# Ensure equal aspect ratio for proper sphere shape
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio for x, y, z axes

# Show plot
plt.show()
