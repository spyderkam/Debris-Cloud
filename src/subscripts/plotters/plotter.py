__date__ = "May 14, 2025"

# 3D scatter plot

# import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

## Assuming vec_r is a list of lists like [[x0, y0, z0], [x1, y1, z1], ...]
vec_r = np.array(vec_r)  # Convert to numpy array for easier handling

## Extract x, y, z coordinates
x = vec_r[:, 0]
y = vec_r[:, 1]
z = vec_r[:, 2]

## Create a new figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

## Create scatter plot
ax.scatter(x, y, z, c='blue', marker='o', s=20)

## Set labels for axes
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
#max_range = max(np.max(np.abs(vec_r)), cloud.radius) * 1.1  # Include data and sphere, with padding
#ax.set_xlim(-max_range, max_range)
#ax.set_ylim(-max_range, max_range)
#ax.set_zlim(-max_range, max_range)

# Ensure equal aspect ratio for proper sphere shape
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio for x, y, z axes

## Show plot
plt.savefig("cloud.png", bbox_inches="tight")
plt.show()
