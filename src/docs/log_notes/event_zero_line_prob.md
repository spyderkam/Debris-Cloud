# Fragment Proximity to a Euclidean Line at Event Zero

Kam Modjtahedzadeh  
Boeing Intelligence & Analytics  
May 20, 2025 -- June 16, 2025

> I have a spherical shell with radius $R_\mathrm c(t)$ centered at $(x, y, z) = (0, 0, 0)$ at $t=0$. I want to find the probability that a particle in the shell will be within a distance $\ell$ away from a straight line passing through that shell.

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

$$\text{inputs:} ~ R, \mathbf{c}, \texttt{d}$$

where:

  - $R \in \mathbb{R}^+$ is the radius of the sphere
  - $\mathbf{c} \in \mathbb{R}^3$ is the center of the sphere, default $(0,0,0)$
  - $\texttt{d} \in \{\texttt{True}, \,\texttt{False}\}$ is the diameter flag, default $\texttt{False}$

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
	-\mathbf{v}_1, & \text{if } \texttt{d} = \texttt{True} \\
	\frac{\mathbf{w}}{\|\mathbf{w}\|_2}, & \text{if } \texttt{d} = \texttt{False}\text{ and } \mathbf{v}_1 \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|_2} > -0.98
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

$$r(\lambda) = (x_1, y_1, z_1) + \lambda \cdot \underbrace{(x_1 - x_2,\, y_1 - y_2,\, z_1 - z_2)}_{\text{direction vector (} \mathbf{k}\mathrm{)}}$$

or, in (parametric) component form:

$$
\begin{align*}
	&X = x_1 + \lambda(x_2 - x_1)  \\
	&Y = y_1 + \lambda(y_2 - y_1)  \\
	&Z = z_1 + \lambda(z_2 - z_1)
\end{align*}
$$

where $\lambda \in \mathbb{R}$ is the scaling parameter. So, if $\lambda = 0$ then $(X, Y, Z) = (x_1, y_1, z_1)$; likewise, if $\lambda = 1 \Rightarrow X = x_2$ and $\lambda = 0.5 \Rightarrow X = (x_1 + x_2)/2$.

Alternatively, if $x_2 \neq x_1$, $y_2 \neq y_1$, and $z_2 \neq z_1$, the line can be written in *symmetric form*:

$$
\frac{X-x_1}{x_2-x_1} = \frac{Y-y_1}{y_2-y_1} = \frac{Z-z_1}{z_2-z_1}
$$

If any denominator is zero (e.g., $x_2 = x_1$), the corresponding coordinate is constant (e.g., $X = x_1$), and that term is omitted from the symmetric form.

## Counting Hits

**<u>_Kam_:</u>**

**What I have:**

1. _`get_entry_exit`_: Generates entry point ($\mathbf p_1 = \mathbf p_{\mathrm{entry}}$) and exit point ($\mathbf p_2 = \mathbf p_{\mathrm{exit}}$).
2. _`line_parametric_3d`_: Computes a parametric equation of the line passing through points $\mathbf p_1$and $\mathbf p_2$ and returns $r(\lambda)$.

**What I want:**

First, I will generate debris particles inside the cloud at time $t = 0$.
 
```python
cloud = Cloud(characteristic_length=0.05, num_fragments=500)
all_points = cloud.sample_positions()
inside_points = [point for point in all_points if np.linalg.norm(point) <= cloud.radius]
```

Next, I will generate entry and exit points and compute the parametric three-dimensional line equation that passes through those points.

```python
p1, p2 = get_entry_exit(cloud.radius)
line = line_parametric_3d(p1, p2)
```

Then I will use the `give_it_a_name` function to calculate and return the integer number of how many points in `inside_points` are within `l` distance from `line`. (This is what I want.)

**<u>_Claude 3.5 Sonnet V2_:</u>**

Let's call this `give_it_a_name`  function `count_points_near_line`. Its implementation:

1.  Takes the parametric line function, array of points, and maximum distance
2.  For each point, finds the closest point on the line using vector projection
3.  Calculates the minimum distance from point to line
4.  Counts points within the specified distance

### Mathematical Foundation of Line-Point Distance Calculation

The `count_points_near_line` function implements a fundamental geometric algorithm for calculating the shortest distance between points and a line in three-dimensional space. This analysis breaks down the mathematical principles underlying the implementation.

#### Orthogonal Projection Method <!--Vector Projection (onto a Line)-->

Recall the parametric line representation,

$$\mathbf{r}(\lambda) = \mathbf{p}_1 + \lambda(\mathbf{p}_2 - \mathbf{p}_1)$$

where the direction vector of the line is $\mathbf{k} = \mathbf{p}_2 - \mathbf{p}_1$. For any point $\mathbf{p}$ in space, finding the closest point on the line requires determining the optimal parameter $\lambda^\star$ that minimizes the distance $|\mathbf{p} - \mathbf{r}(\lambda)|$. 

##### Mathematical Derivation of $\lambda^\star$

To find the optimal parameter, instead of minimizing the distance, minimize the _squared_ distance between $\mathbf{p}$ and the line point $\mathbf{r}(\lambda)$. The condition for minimization requires that the vector from $\mathbf{p}$ to the closest line point be perpendicular to the direction vector $\mathbf{k}$.

Setting up the _perpendicularity condition_:

$$(\mathbf{p} - \mathbf{r}(\lambda^\star)) \cdot \mathbf{k} = 0$$

substituting the parametric equation:

$$(\mathbf{p} - \mathbf{p}_1 - \lambda^\star\mathbf{k}) \cdot \mathbf{k} = 0$$

expanding the dot product:

$$(\mathbf{p} - \mathbf{p}_1) \cdot \mathbf{k} - \lambda^\star(\mathbf{k} \cdot \mathbf{k}) = 0$$

The vector from $\mathbf{p}_1$ to the query point $\mathbf{p}$  (displacement vector) is defined as $\mathbf{q} = \mathbf{p} - \mathbf{p}_1$;

$$\mathbf{q} \cdot \mathbf{k} - \lambda^\star|\mathbf{k}|^2 = 0$$

solving for $\lambda^\star$:

$$\lambda^\star = \frac{\mathbf{q} \cdot\mathbf{k}}{|\mathbf{k}|^2} = \frac{(\mathbf{p} - \mathbf{p}_1) \cdot (\mathbf{p}_2 - \mathbf{p}_1)}{|\mathbf{p}_2 - \mathbf{p}_1|^2}$$

This formula represents the scalar projection of $\mathbf{q}$ onto $\mathbf{k}$, normalized by the squared magnitude of $\mathbf{k}$; in other words,

$$\lambda^\star = \frac{\text{scalar projection of } \mathbf{q} \text{ onto } \mathbf{k}}{|\mathbf{k}|} $$

#### Distance Calculation

Once the optimal parameter is determined, the closest point on the line is:

$$\mathbf{r}(\lambda^\star) = \mathbf{p}_1 + \lambda^\star(\mathbf{p}_2 - \mathbf{p}_1)$$

The minimum distance from point $\mathbf{p}$ to the line is the Euclidean distance between $\mathbf{p}$ and this closest point:

$$l_{\mathrm{min}} = |\mathbf{p} - \mathbf{r}(\lambda^\star)| = \sqrt{\sum_{i=1}^{3}(p_i - r_i(\lambda^\star))^2}$$

### Computational Implementation

#### Implementation Analysis

The function implementation follows this mathematical framework precisely. The code extracts $\mathbf{p}_1$ and $\mathbf{p}_2$ by evaluating the parametric line function at $\lambda = 0$ and $\lambda = 1$ respectively. It then computes the direction vector $\mathbf{k} = \mathbf{p}_2 - \mathbf{p}_1$ and the displacement vector $\mathbf{q} = \mathbf{p} - \mathbf{p}_1$.

The projection coefficient calculation implements the dot product formula $\lambda^\star = \mathbf{q} \cdot \mathbf{k} / |\mathbf{k}|^2$ using component-wise operations. The closest point determination and final distance calculation follow the standard Euclidean distance formula. If $l_\mathrm{min} \leq \ell$, then that point is counted as a hit. 

#### Computational Complexity

The algorithm exhibits linear time complexity $\mathcal O(n)$ with respect to the number of points, as each point requires a constant number of operations regardless of the problem size. The space complexity is $\mathcal O(1)$ since the function processes points individually without storing intermediate results.

This mathematical approach provides an exact solution to the point-to-line distance problem, making it suitable for precise proximity analysis in debris cloud modeling applications where accurate geometric calculations are essential for risk assessment and trajectory planning.​​​​​​​​​​​​​​​​

## Analytical Framework for Impact Probability Calculation

Given the parent mass, initial cloud radius, cumulative distribution function, $\rho(r,\, L_{\mathrm{c}})$, and hit distance threshold, the impact probability will be calculated such that a linear trajectory passing through random entry and exit points on the cloud sphere at time $t = 0$ will encounter at least one fragment within distance $\ell$. A fundamental result in stochastic processes from *Poisson process theory* is the application of a *non-homogeneous Poisson process*; since encounters occur randomly at rate (per unit length) $\Lambda(\mathcal{L}')$ along path $\mathcal{L}'$, then the number of encounters follows a Poisson distribution with parameter $\int\Lambda d\mathcal{L}'$. The probability of at least one encounter is $1 - \mathbb{P}_{\text{avoid}} = 1 - \exp[-\int\Lambda\, d\mathcal{L}']$. The total impact probability is computed as the weighted sum over all fragment size categories: $$ 
\begin{align*}
    \mathbb{P}_\text{impact} \,=\, 
    %%%%%  
    \sum_{L_\mathrm c \,=\, L_\mathrm{min} }^{L_{\mathrm{max}}} f \cdot \mathbb{P}_\text{impact}(L_{\mathrm{c}})
\end{align*}
$$ where $f$ represents the fraction of fragments given from percentage distributions per size category. For a specific characteristic length, the probability of impact along a trajectory is: $$
\begin{align*}
    \mathbb{P}_{\text{impact}}(L_\mathrm c,\, \mathcal{L}) \,=\, 1 - \exp\left[-\int_\mathcal{L}\Lambda (L_\mathrm c,\, \mathcal{L}^\prime) \,d\mathcal{L}^\prime\right]
\end{align*}
$$ The collision rate at position $s$ along trajectory $\mathcal{L}$ is: $$\Lambda(s,\, \mathcal{L}) \,=\, \rho_N(r(s),\, \mathcal{L}) \cdot A_{\mathrm{eff}}$$ where $A_{\mathrm{eff}} = \pi\ell^2$ is the effective impact cross-sectional area and $\rho_N(r)$ is the the number density of fragments at radial $r$. The number density is related to the mass density through, $$
\begin{align*}
    \rho(r, L_\mathrm{c}) \,=\, \frac{\rho(r, L_\mathrm{c})}{\overline{M}} \frac{dN}{dL_\mathrm{c}} 
\end{align*}
$$ where $\overline{M}$ is the average mass of fragments with characteristic length $L_\mathrm{c}$.

### Monte Carlo Implementation

#### Basic Monte Carlo Estimator

The impact probability is estimated through repeated sampling:

$$ \hat{\mathbb{P}}_{\text{impact}} \,=\, \frac{1}{\Upsilon_{\mathrm{trials}}} \sum_{j=1}^{\Upsilon_{\mathrm{trials}}} I_j $$

where $I_j$ is an indicator function:

$$ I_j \,=\, \begin{cases} 1 & \text{if trajectory } j \text{ encounters any fragment} \\ 0 & \text{otherwise} \end{cases} $$

#### Variance and Confidence Intervals

The variance of the Monte Carlo estimator is:

$$ \mathrm{Var}\!\left(\hat{\mathbb{P}}\right) \,=\, \frac{\mathbb{P}(1-\mathbb{P})}{\Upsilon_{\mathrm{trials}}} $$

The confidence interval using the Wilson score method is:

$$ \mathbb{P} \,\in\, \left[\frac{\hat{\mathbb{P}} + \frac{\zeta^2}{2\Upsilon} - \zeta\sqrt{\frac{\hat{\mathbb{P}}(1-\hat{\mathbb{P}})}{\Upsilon} + \frac{\zeta^2}{4\Upsilon^2}}}{1 + \frac{\zeta^2}{\Upsilon}},\ \frac{\hat{\mathbb{P}} + \frac{\zeta^2}{2\Upsilon} + \zeta\sqrt{\frac{\hat{\mathbb{P}}(1-\hat{\mathbb{P}})}{\Upsilon} + \frac{\zeta^2}{4\Upsilon^2}}}{1 + \frac{\zeta^2}{\Upsilon}}\right] $$

where $\zeta$ is the z-score (standard normal quantile).

#### Adaptive Sample Size with Sequential Refinement

Since $\mathbb{P}$ is unknown a priori, an adaptive sampling approach is employed. Starting with an initial batch of $\Upsilon_0 = 10^4$ trials, the required sample size is updated after each batch:

$$ \Upsilon_{\mathrm{next}} \,=\, \frac{\zeta^2 (1-\hat{\mathbb{P}}_{\mathrm{current}})}{\hat{\mathbb{P}}_{\mathrm{current}} \epsilon^2} $$

The sampling continues until convergence is achieved:

$$ \frac{\mathrm{CI}_{\mathrm{upper}} - \mathrm{CI}_{\mathrm{lower}}}{\hat{\mathbb{P}}_{\mathrm{current}}} \,<\, \epsilon $$

where $\mathrm{CI}$ is confidence interval given from the Wilson score method. As more trials are added, the $\sqrt{\Upsilon}$ in the denominator of the Wilson formula makes the interval tighter, reducing the relative width until it's less than the tolerance $\epsilon$.

### Importance Sampling Enhancement

#### Biased Sampling Distribution

Instead of uniform sampling on the sphere, use importance sampling where samples are drawn from a biased distribution (called the importance distribution) that prioritizes "important" regions of the sample space—those more likely to contribute to the quantity being estimated, and then adjust for this bias by weighting the samples. To sample entry and exit points this way, use the following probability density function (PDF):

$$ \mathcal{P}(\mathbf{p}_{\mathrm{entry}}, \mathbf{p}_{\mathrm{exit}}) \,\propto\, \exp\left[-\frac{\ell_{\mathrm{min}}^2(\mathbf{p}_{\mathrm{entry}}, \mathbf{p}_{\mathrm{exit}})}{2\sigma_{\mathrm{IS}}^2}\right] $$

where $\ell_{\mathrm{min}}$ is a function that calculates the minimum distance between the trajectory defined by these points and the peak density sphere at radius $\mu R_{\mathrm{c}}$ and $\sigma_{\mathrm{IS}}$ is the *importance sampling* standard deviation, a tuning parameter that controls how strongly the sampling is biased toward trajectories passing near the peak density radius.[^1]


[^1]: A small $\sigma_{\mathrm{IS}}$ means we heavily favor trajectories close to $\mu R_{\mathrm{c}}$, a large $\sigma_{\mathrm{IS}}$ means the bias is weaker, approaching uniform sampling. Typically, $\sigma_{\mathrm{IS}} \approx 0.2 R_\mathrm c$ provides a good balance between bias strength and sampling efficiency, though this can be adjusted based on the specific debris distribution.


##### Why This Shape

1. **When $\ell_{\mathrm{min}} = 0$** (trajectory passes exactly through the peak density shell):
    - The exponent is 0
    - $e^0 = 1$ (maximum value)
    - These trajectories get the highest sampling probability
2. **As $\ell_{\mathrm{min}}$ increases** (trajectory gets farther from peak density):
    - The exponent becomes more negative
    - The value decreases toward 0
    - These trajectories get lower sampling probability
3. **The role of $\sigma_{\mathrm{IS}}$**:
    - Small $\sigma_{\mathrm{IS}}$: Sharp peak, strong bias (only trajectories very close to peak get sampled)
    - Large $\sigma_{\mathrm{IS}}$: Gentle peak, weak bias (more like uniform sampling)

##### Why Gaussian

We choose this Gaussian form because:

  1. **It matches your debris distribution** - your fragments follow a Gaussian density distribution, so it makes sense to use a Gaussian importance function
  2. **Smooth falloff** - trajectories slightly farther from the peak still get sampled, just less frequently
  3. **Well-understood mathematics** - Gaussian functions have nice properties for analysis
  4. **Tunable** - the $\sigma_{\mathrm{IS}}$ parameter lets you control how aggressive the biasing is

#### Weighted Estimator

The importance-sampled estimator becomes:

$$ \hat{\mathbb{P}}_{\mathrm{IS}} \,=\, \frac{1}{\Upsilon_{\mathrm{trials}}} \sum_{j=1}^{\Upsilon_{\mathrm{trials}}} I_j w_j $$

where the weight is:

$$ w_j = \frac{h(\mathbf{p}_{\mathrm{entry}}, \mathbf{p}_{\mathrm{exit}})}{\mathcal{P}(\mathbf{p}_{\mathrm{entry}}, \mathbf{p}_{\mathrm{exit}})} $$

where and $h$ represents the uniform distribution on the sphere.

## Sample Results for Large Asteroid

```
############################################################
IMPACT PROBABILITY AT EVENT ZERO
Monte Carlo Implementation of Collision Probability Analysis
############################################################

Creating debris cloud...
Parent mass: 10000 kg
Parent radius: 1000 m
Hit distance threshold: 1.0 m

Cloud created in 14.08 seconds
Total fragments: 804,105
Cloud radius: 1100.64 m 

Fragment distribution (inside cloud radius only):
  Small fragments (< 8 cm):  804,077 (99.997%)
  Medium fragments (8-11 cm): 28 (0.003%)
  Large fragments (> 11 cm):  0 (0.000%)

############################################################
MONTE CARLO IMPACT PROBABILITY CALCULATION
############################################################

1. Standard Monte Carlo Estimation:
Starting Monte Carlo simulation with 10000 trials...        
Using 40,205 sampled fragments (5.0% of total)
Trial 0/10000 (0.0%) - ETA: 0.0s
Trial 1000/10000 (10.0%) - ETA: 25.3s
Trial 2000/10000 (20.0%) - ETA: 22.3s
Trial 3000/10000 (30.0%) - ETA: 19.0s
Trial 4000/10000 (40.0%) - ETA: 16.1s
Trial 5000/10000 (50.0%) - ETA: 13.4s
Trial 6000/10000 (60.0%) - ETA: 10.7s
Trial 7000/10000 (70.0%) - ETA: 8.0s
Trial 8000/10000 (80.0%) - ETA: 5.3s
Trial 9000/10000 (90.0%) - ETA: 2.7s

Results:
  Impact Probability: 0.005200
  Hits: 52 out of 10,000 trials
  Fragments Used: 40,205 (5.0% of total)
  95% Confidence Interval: [0.003968, 0.006812]
  Standard Error: 0.000719
  Computation Time: 26.48 seconds

2. Adaptive Monte Carlo with Sequential Refinement:        
Starting adaptive Monte Carlo with target precision 5.0%...
Starting Monte Carlo simulation with 1000 trials...        
Using 80,410 sampled fragments (10.0% of total)
Trial 0/1000 (0.0%) - ETA: 0.0s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.3s
Trial 2000/5000 (40.0%) - ETA: 15.1s
Trial 3000/5000 (60.0%) - ETA: 10.3s
Trial 4000/5000 (80.0%) - ETA: 5.9s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 21.0s
Trial 2000/5000 (40.0%) - ETA: 15.5s
Trial 3000/5000 (60.0%) - ETA: 10.6s
Trial 4000/5000 (80.0%) - ETA: 5.8s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 21.3s
Trial 2000/5000 (40.0%) - ETA: 16.8s
Trial 3000/5000 (60.0%) - ETA: 11.6s
Trial 4000/5000 (80.0%) - ETA: 6.1s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.7s
Trial 2000/5000 (40.0%) - ETA: 15.2s
Trial 3000/5000 (60.0%) - ETA: 10.1s
Trial 4000/5000 (80.0%) - ETA: 5.0s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.6s
Trial 2000/5000 (40.0%) - ETA: 15.2s
Trial 3000/5000 (60.0%) - ETA: 10.0s
Trial 4000/5000 (80.0%) - ETA: 4.9s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 25.5s
Trial 2000/5000 (40.0%) - ETA: 17.4s
Trial 3000/5000 (60.0%) - ETA: 11.0s
Trial 4000/5000 (80.0%) - ETA: 5.3s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.6s
Trial 2000/5000 (40.0%) - ETA: 13.8s
Trial 3000/5000 (60.0%) - ETA: 9.2s
Trial 4000/5000 (80.0%) - ETA: 4.6s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.2s
Trial 2000/5000 (40.0%) - ETA: 13.5s
Trial 3000/5000 (60.0%) - ETA: 9.0s
Trial 4000/5000 (80.0%) - ETA: 4.5s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.3s
Trial 2000/5000 (40.0%) - ETA: 13.6s
Trial 3000/5000 (60.0%) - ETA: 9.1s
Trial 4000/5000 (80.0%) - ETA: 4.7s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 19.5s
Trial 2000/5000 (40.0%) - ETA: 14.1s
Trial 3000/5000 (60.0%) - ETA: 9.6s
Trial 4000/5000 (80.0%) - ETA: 5.0s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.1s
Trial 2000/5000 (40.0%) - ETA: 15.1s
Trial 3000/5000 (60.0%) - ETA: 10.0s
Trial 4000/5000 (80.0%) - ETA: 5.0s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.1s
Trial 2000/5000 (40.0%) - ETA: 14.4s
Trial 3000/5000 (60.0%) - ETA: 9.5s
Trial 4000/5000 (80.0%) - ETA: 4.8s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.7s
Trial 2000/5000 (40.0%) - ETA: 14.5s
Trial 3000/5000 (60.0%) - ETA: 9.7s
Trial 4000/5000 (80.0%) - ETA: 4.8s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.3s
Trial 2000/5000 (40.0%) - ETA: 13.5s
Trial 3000/5000 (60.0%) - ETA: 9.3s
Trial 4000/5000 (80.0%) - ETA: 4.7s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.2s
Trial 2000/5000 (40.0%) - ETA: 14.3s
Trial 3000/5000 (60.0%) - ETA: 9.7s
Trial 4000/5000 (80.0%) - ETA: 4.9s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 21.2s
Trial 2000/5000 (40.0%) - ETA: 16.9s
Trial 3000/5000 (60.0%) - ETA: 10.9s
Trial 4000/5000 (80.0%) - ETA: 5.3s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.7s
Trial 2000/5000 (40.0%) - ETA: 15.0s
Trial 3000/5000 (60.0%) - ETA: 9.7s
Trial 4000/5000 (80.0%) - ETA: 4.8s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.7s
Trial 2000/5000 (40.0%) - ETA: 14.3s
Trial 3000/5000 (60.0%) - ETA: 9.7s
Trial 4000/5000 (80.0%) - ETA: 4.9s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 19.3s
Trial 2000/5000 (40.0%) - ETA: 14.5s
Trial 3000/5000 (60.0%) - ETA: 9.6s
Trial 4000/5000 (80.0%) - ETA: 4.8s
Adding 5000 more trials...
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.4s
Trial 2000/5000 (40.0%) - ETA: 14.9s
Trial 3000/5000 (60.0%) - ETA: 9.8s
Trial 4000/5000 (80.0%) - ETA: 4.8s

Adaptive Results:
  Impact Probability: 0.012941
  Hits: 1,307 out of 101,000 trials
  95% Confidence Interval: [0.012262, 0.013656]
  Total Computation Time: 5.10 seconds

3. Sensitivity Analysis - Different Hit Distances:
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.3s
Trial 2000/5000 (40.0%) - ETA: 13.6s
Trial 3000/5000 (60.0%) - ETA: 9.0s
Trial 4000/5000 (80.0%) - ETA: 4.5s
  Hit distance 0.5 m: ℙ = 0.002400 ± 0.000692
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 20.2s
Trial 2000/5000 (40.0%) - ETA: 15.1s
Trial 3000/5000 (60.0%) - ETA: 10.1s
Trial 4000/5000 (80.0%) - ETA: 5.0s
  Hit distance 1.0 m: ℙ = 0.014000 ± 0.001662
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 17.8s
Trial 2000/5000 (40.0%) - ETA: 13.4s
Trial 3000/5000 (60.0%) - ETA: 8.9s
Trial 4000/5000 (80.0%) - ETA: 4.4s
  Hit distance 2.0 m: ℙ = 0.053600 ± 0.003185
Starting Monte Carlo simulation with 5000 trials...
Using 80,410 sampled fragments (10.0% of total)
Trial 0/5000 (0.0%) - ETA: 0.0s
Trial 1000/5000 (20.0%) - ETA: 18.3s
Trial 2000/5000 (40.0%) - ETA: 13.6s
Trial 3000/5000 (60.0%) - ETA: 9.0s
Trial 4000/5000 (80.0%) - ETA: 4.5s
  Hit distance 5.0 m: ℙ = 0.258400 ± 0.006191

############################################################
ANALYSIS COMPLETE
############################################################
```

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTI4OTYzMzE5Ml19
-->