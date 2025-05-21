# Fragment Proximity to a Euclidean Line at Event Zero

Kam Modjtahedzadeh  
Boeing Intelligence & Analytics  
May 20, 2025 -- May 21, 2025

> I have a spherical shell with radius $R_\mathrm c(t)$ centered at $(x, y, z) = (0, 0, 0)$ at $t=0$. I want to find the probability that a particle in the shell will be within $r$ meters of a straight line passing through that shell.

## Creating a Sphere Element in Entry/Exit Point Generator

**<u>_Kam_:</u>**

I want a Python function that generates two random Cartesian points on the surface of a sphere with radius `radius` (entry and exit). This function has three inputs: a *required* input, the `radius`, an *optional* input, the `center=(0,0,0)` another *optional* input called `diameter=False`. If `True`, then the endpoints have the possibility of creating a diameter as the line that passes through `entry_point` and `exit_point` also passes through the center of the sphere, i.e., a diameter of the spherical shell. If it is set to its default `False`, then `entry_point`  and `exit_point` *cannot* be the endpoints of a diameter. Needless to say by now, but the output of this function is to be a tuple: `(entry_point, exit_point)`; name the function appropriately.

**<u>_Claude Sonnet 3.7_:</u>**

The ***`get_entry_exit`*** function uses the Gaussian sampling method to generate uniformly distributed points on a sphere, similar to the approach mentioned in your document on spherically symmetric debris cloud models.

Key features:

  1.  It generates the first point using a normalized vector from a 3D Gaussian distribution
  2.  For `diameter=True`, it creates the second point exactly opposite to the first
  3.  For `diameter=False`, it ensures the points aren't diametrically opposite by checking their dot product
  4.  Both points are scaled by the radius and shifted according to the center coordinates
  5. The function uses the dot product threshold of $-0.98$ to ensure the points aren't nearly diametrically opposite. This corresponds to avoiding angles greater than about $168°$ between the points, since $\cos(168°) ≈ -0.98$.

The approach aligns with the computational methods detailed in Section 2 of the "Computational Implementation of Spherically Symmetric Debris Cloud Models" document.

### Mathematical Formulation of `get_entry_exit`

$$\text{inputs:} ~ R, \mathbf{c}, d$$

where:

  - $R \in \mathbb{R}^+$ is the radius of the sphere
  - $\mathbf{c} \in \mathbb{R}^3$ is the center of the sphere, default $(0,0,0)$
  - $d \in \{\texttt{True}, \,\texttt{False}\}$ is the diameter flag, default False

$$\text{output:} ~ (\mathbf{p}_{\text{entry}}, \mathbf{p}_{\text{exit}}) \in \mathbb{R}^3 \times \mathbb{R}^3$$

#### Algorithm

1. Generate a random unit vector $\mathbf u$ on the unit sphere and normalize it to be $\mathbf{v}_1$:

$$
\begin{align*}
	&\mathbf{u} \sim \mathrm{Normal}(0, 1)^3  \\
	&\mathbf{v}_1 = \frac{\mathbf{u}}{\|\mathbf{u}\|_2}
\end{align*}
$$

2. Compute the second vector $\mathbf{v}_2$ based on the diameter flag: 
$$
\mathbf{v}_2 = \begin{cases}
	-\mathbf{v}_1, & \text{if } d = \texttt{True} \\
	\frac{\mathbf{w}}{\|\mathbf{w}\|_2}, & \text{if } d = \texttt{False}\text{ and } \mathbf{v}_1 \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|_2} > -0.98
\end{cases}
$$ where $\mathbf{w} \sim \mathrm{Normal}(0, 1)^3$ is regenerated until the dot product condition is satisfied.

4. Compute the entry and exit points:

$$
\begin{align*}
	\mathbf{p}_{\text{entry}} = \mathbf{c} + R\mathbf{v}_1 \\
	\mathbf{p}_{\text{exit}} = \mathbf{c} + R\mathbf{v}_2
\end{align*}
$$

4. Return the tuple $(\mathbf{p}_{\text{entry}}, \mathbf{p}_{\text{exit}})$

#### Computational Details

- The function uses the property that normalizing vectors sampled from a multivariate normal distribution $\mathrm{Normal}(0, 1)^3$ creates a uniform distribution on the unit sphere.

- For the non-diameter case, the constraint $\mathbf{v}_1 \cdot \mathbf{v}_2 > -0.98$ ensures the vectors are not nearly antipodal, as: $$\mathbf{v}_1 \cdot \mathbf{v}_2 = \cos(\vartheta)$$ where $\vartheta$ is the angle between the vectors. The value $-0.98$ corresponds to approximately $\vartheta < 168^{\circ}$.​​​​​​​​​​​​​​​​


## Parametric 3D Line

**<u>_Kam_:</u>**

I have two points $\mathbf{p}_1 = (x_1, y_1, z_1)$ and $\mathbf{p}_2 = (x_2, y_2, z_2)$ in 3D Euclidean space. I want the equation of the line that passes through $\mathbf{p}_1$ and $\mathbf{p}_2$.

**<u>_Grok 3_:</u>**

Given two points $\mathbf{p}_1 = (x_1, y_1, z_1)$ and $\mathbf{p}_2 = (x_2, y_2, z_2)$ in 3D Euclidean space, the line passing through them can be expressing in *parametric form*:

$$r(\lambda) = (x_1, y_1, z_1) + \lambda (x_1 - x_2,\, y_1 - y_2,\, z_1 - z_2)$$

or, in component form:

$$
\begin{align*}
	&X = x_1 + \lambda(x_2 - x_1)  \\
	&Y = y_1 + \lambda(y_2 - y_1)  \\
	&Z = z_1 + \lambda(z_2 - z_1)
\end{align*}
$$

where $\lambda \in \mathbb{R}$ is the parameter. So, if $\lambda = 0$ then $(X, Y, Z) = (x_1, y_1. y_1)$; likewise, if $\lambda = 1 \Rightarrow X = x_2$ and $\lambda = 0.5 \Rightarrow X = (x_1 + x_2)/2$.

Alternatively, if $x_2 \neq x_1$, $y_2 \neq y_1$, and $z_2 \neq z_1$, the line can be written in *symmetric form*:

$$
\frac{X-x_1}{x_2-x_1} = \frac{Y-y_1}{y_2-y_1} = \frac{Z-z_1}{z_2-z_1}
$$

If any denominator is zero (e.g., $x_2 = x_1$), the corresponding coordinate is constant (e.g., $X = x_1$), and that term is omitted from the symmetric form.

***

**What I have:**
1. 
