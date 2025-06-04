# Computational Implementation of Spherically Symmetric Debris Cloud Models

Kamyar Modjtahedzadeh <br>
Boeing Intelligence & Analytics <br>
May 13, 2025

## 1 Distribution of Area-to-Mass Ratios

Based on the NASA breakup model's area-to-mass ratio distribution [3], a piecewise statistical function has been structured for the fragment's logarithm characteristic length ($L_c$) input and its cross-sectional area-to-mass ratio ($A/M$) output. The function partitions the domain into three regions based on fragment size:

$$
\xi (x) \,=\,
\begin{cases}
\xi_{\text{small}}(\phi) & \text{if } \phi < \log(0.08) \\
\xi_{\text{inter}}(\phi) & \text{if } \log(0.08) \leq \phi \leq \log(0.11) \\
\xi_{\text{large}}(\phi) & \text{if } \phi > \log(0.11)
\end{cases} \tag{1.1}
$$

where $\phi \equiv \log L_c$, $\xi_{\text{small}}$ is the small fragment distribution and $\xi_{\text{large}}$ is the large fragment distribution.

### 1.1 Interpolation of Middle Region

For the middle region (8 cm ≤ $L_c$ ≤ 11 cm), the $\xi$ function uses linear interpolation;

$$\xi_{\text{inter}}(\phi) = \xi_{\text{small}}(\phi) \cdot (1-w(\phi)) + \xi_{\text{large}}(\phi) \cdot w(\phi) \tag{1.2}$$

where the weight function is given by,

$$w(\phi) = \frac{\phi - \log(0.08)}{\log(0.11) - \log(0.08)} \tag{1.3}$$

This creates a smooth transition between the statistical distributions for small and large fragments [3].[^1]

### 1.2 Statistical Sampling Functions

Each of the base functions $\xi_{\text{small}}$ and $\xi_{\text{large}}$ perform statistical sampling from normal distributions in logarithmic space. The output value is ultimately transformed from logarithmic space back to linear space. This mathematical structure creates a continuous statistical model for $A/M$ values that accommodates different physical regimes across the size spectrum of orbital debris fragments [3].

#### 1.2.1 Small Fragment Distribution

The small fragment distribution function $\xi_{\text{small}}(\phi)$ implements a unimodal log-normal sampling model:

$$\xi_{\text{small}}(\phi) \,=\, 10^{N(\bar{\chi},\, \sigma_{\chi})}  \tag{1.4}$$

where $N(\bar{\chi}, \sigma_{\chi})$ represents a random sample from a normal distribution, $\bar{\chi}$ is the mean of $\chi \equiv \log \frac{A}{M}$, and $\sigma_{\chi}$ is the standard deviation of $\chi$. The aforementioned normal distribution is:

$$\frac{1}{\sigma_{\chi}\sqrt{2\pi}} \exp\left[- \frac{(\chi - \bar{\chi})^2}{2\sigma_{\chi}^2}\right] \tag{1.5}$$

in Johnson et al. (2001), this normal distribution models the logarithm of area-to-mass ratios [3]. In other words, $\log(A/M)$ peaks and is symmetric around that peak value of $\bar{\chi}$.

Moreover, $\bar{\chi}$ and $\sigma_{\chi}$ follow piecewise functions. The expression for $\bar{\chi}(\phi)$ is:

$$
\bar{\chi}(\phi) \,=\,
\begin{cases}
-0.3 & \text{if } \phi \leq -1.75 \\
-0.3 - 1.4(\phi + 1.75) & \text{if } -\!
1.75 < \phi < -1.25 \\
-1.0 & \text{if } \phi \geq -1.25,
\end{cases} \tag{1.6}
$$

and for the standard deviation of $\chi$:

$$
\sigma_{\chi}(\phi) \,=\,
\begin{cases}
0.2 & \text{if } \phi \leq -3.5 \\
0.2 + 0.1333(\phi + 3.5) & \text{if } \phi > -3.5.
\end{cases} \tag{1.7}
$$

#### 1.2.2 Large Fragment Distribution

The large fragment distribution implements a bimodal log-normal sampling model [3];

$$
\xi_{\text{large}}(\phi) \,=\,
\begin{cases}
10^{N(\bar{\chi}_1,\, \sigma_{\chi_1})} & \text{with probability } \omega(\phi) \\
10^{N(\bar{\chi}_2,\, \sigma_{\chi_2})} & \text{with probability } 1 - \omega(\phi)
\end{cases} \tag{1.8}
$$

where the weighting parameter, $\omega$, determines how the bimodal distribution splits fragments between two different normal distributions. The term bimodal captures the fact that larger fragments tend to come from two different populations with distinct physical properties: Fragments with higher and lower $A/M$ values. Particles with higher area-to-mass ratios—represented by $N(\bar{\chi}_1, \sigma_{\chi_1})$—are typically less dense and from planar material components while those with lower ratios—$N(\bar{\chi}_2, \sigma_{\chi_2})$—usually more dense and from compact material. This approach is consistent with the physical reality that parent bodies are composed of materials with varying densities and structural properties.

All parameters are piecewise functions of $\phi$. For upper stage fragments, the weight parameter:

$$
\omega(\phi) \,=\,
\begin{cases}
1.0 & \text{if } \phi \leq -1.4 \\
1.0 - 0.3571(\phi + 1.4) & \text{if } -\!1.4 < \phi < 0 \\
0.5 & \text{if } \phi \geq 0
\end{cases} \tag{1.9}
$$

the first distribution mean:

$$
\bar{\chi}_1(\phi) \,=\,
\begin{cases}
-0.45 & \text{if } \phi \leq -0.5 \\
-0.45 - 0.9 \cdot (\phi + 0.5) & \text{if } -\!0.5 < \phi < 0 \\
-0.9 & \text{if } \phi \geq 0,
\end{cases} \tag{1.10}
$$

the first distribution standard deviation:

$$\sigma_{\chi_1}(\phi) \,=\,0.55 \tag{1.11}$$

the second distribution mean:

$$\bar{\chi}_2(\phi) \,=\,-0.9 \tag{1.12}$$

and the second distribution standard deviation:

$$
\sigma_{\chi_2}(\phi) \,=\,
\begin{cases}
0.28 & \text{if } \phi \leq -1.0 \\
0.28 - 0.1636 \cdot (\phi + 1) & \text{if } -\!1.0 < \phi < 0.1 \\
0.1 & \text{if } \phi \geq 0.1
\end{cases} \tag{1.13}
$$

## 2 Sampling in 3D Space with a Spherically Symmetric Density

Given the spherically symmetric volumetric mass density ($\tilde{\rho}$) of a debris cloud, fragment positions can be generated that must follow that distribution. The statistical sampling here requires consideration of the spherical geometry to accurately represent the spatial distribution of fragments.

### 2.1 Spherical Volume Element in Probability Sampling

The spherical volume element is:

$$dV \,=\, r^2 \sin\theta\, dr\,d\theta\,d\varphi \tag{2.1}$$

The normalized probability density function (PDF)—i.e., probability per unit volume—is [1]:

$$\mathcal{P}(r) \,=\, \frac{\tilde{\rho}(r)}{M_{\text{total}}} \tag{2.2}$$

where $M_{\text{total}}$ is the total mass of all the fragments. The probability of finding a fragment in a specific volume element is then proportional to $\tilde{\rho}\,dV$;

$$\mathbb{P}_{r,\theta,\varphi} \,\propto\, \tilde{\rho}(r) \cdot r^2 \sin\theta \tag{2.3}$$

The statistical approach presented here reflects techniques developed in the NASA EVOLVE 4.0 framework for modeling fragmentation debris clouds in orbital environments [3].

### 2.2 Marginal and Conditional PDFs

To sample a three-dimensional point in spherical coordinates, the marginal PDFs for the radial and angular components.[^2]

#### 2.2.1 Radial PDF

Following (2.3), normalizing the radial probability—i.e., the radial component of $\mathbb{P}_{r,\theta,\varphi}$—requires dividing it by the radial segment of the total mass [1, 2];


$$\int_0^{\infty} \mathbb{P}_r(r)\,dr \,=\,\int_0^{\infty} \frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\, dr'}\,dr \,=\,1 \tag{2.4}$$

The marginal PDF for $r$ is obtained by integrating $\mathbb{P}$ over $\theta$ and $\varphi$:

$$\mathcal{P}_r(r) = \int_0^{2\pi} \int_0^{\pi} \mathbb{P}_{r,\theta,\varphi}\, d\theta\, d\varphi \tag{2.5}$$

$$= \int_0^{2\pi} \int_0^{\pi} \frac{\tilde{\rho}(r) \cdot r^2 \sin\theta}{4\pi \int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\, dr'}\, d\theta\, d\varphi \tag{2.6}$$

$$= \frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\, dr'} \tag{2.7}$$

The radial PDF, which shows that the probability density for the radial distance is proportional to $\tilde{\rho} \times r^2$, reflecting the increasing volume of spherical shells at larger $r$.

#### 2.2.2 Angular PDFs

The joint angular PDF, conditioned on $r$, is:

$$\mathcal{P}_{\theta,\varphi}(r) = \frac{\mathbb{P}_{r,\theta,\varphi}}{\mathbb{P}_r(r)} \tag{2.8}$$

$$= \left(\frac{\tilde{\rho}(r) \cdot r^2 \sin\theta}{4\pi \int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\, dr'}\right) \bigg/ \left(\frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\, dr'}\right) \tag{2.9}$$

$$= \frac{\sin\theta}{4\pi} \tag{2.10}$$

Line (2.9) can be factored using separation of variables;

$$\mathcal{P}_{\theta,\varphi} = \mathcal{P}_{\theta}(\theta) \cdot \mathcal{P}_{\varphi}(\varphi) \tag{2.11}$$

where $\mathcal{P}_{\theta}$ and $\mathcal{P}_{\varphi}$ are given by,

$$\mathcal{P}_{\theta}(\theta) = \frac{\sin\theta}{2} \text{ for } \theta \in [0, \pi] \tag{2.12}$$

$$\mathcal{P}_{\varphi}(\varphi) = \frac{1}{2\pi} \text{ for } \varphi \in [0, 2\pi) \tag{2.13}$$

The sine term in (2.12) indicates that $\theta$ is not uniformly distributed; instead, the polar angle distribution is weighted by the solid angle element.

### 2.3 Sampling Procedure

To sample a position in spherical coordinates, each coordinate must be generated according to its respective PDF.

#### 2.3.1 Azimuthal Angle

The azimuthal angle $\varphi$ follows a uniform distribution from (2.13);

$$\varphi \sim \text{Uniform}[0, 2\pi) \tag{2.14}$$

$$U_1 \sim \text{Uniform}[0, 1] \tag{2.15}$$

where Uniform refers to the uniform probability distribution.[^3] $\varphi$ can then be written in terms of $U_1$:

$$\varphi \,=\, 2\pi U_1 \tag{2.16}$$

#### 2.3.2 Polar Angle

The polar angle's PDF is given from (2.12). To sample $\theta$, use the cumulative distribution function (CDF):

$$F_{\theta}(\theta) \,=\, \int_0^{\theta} \mathcal{P}_{\theta}(\theta')\,d\theta' = \frac{1 - \cos\theta}{2}, \tag{2.17}$$  

then set $F_{\theta}(\theta) = U_2$ where $U_2 \sim \text{Uniform}[0, 1]$,

$$\theta \,=\, \arccos(1 - 2U_2). \tag{2.18}$$

Alternatively, since $\cos\theta \in [-1, 1]$ and $U_2$ follows a uniform distribution on $[0, 1]$, a standard linear transformation for converting a uniform random variable from one range to another can be used. The formula to be used is:

$$X \,=\, a + (b - a)U, \tag{2.19}$$

where $X$ is uniformly distributed on $[a, b]$ [4, 6].[^4] Since $a$ and $b$ are respectively equal to $-1$ and $+1$, $X \to \cos\theta = 2U_2 - 1$; and so,

$$\theta \,=\, \arccos(2U_2 - 1) \tag{2.20}$$

This ensures the angular distribution is isotropic over the sphere.

#### 2.3.3 Radial Distance

The PDF for $r$ is:

$$\mathcal{P}_r(r) \,=\,\frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2\,dr'} \text{ for } r \in [0, \infty) \tag{2.21}$$

sampling from this PDF depends on the form of $\tilde{\rho}(r)$. Sampling from this PDF depends on the form of $\tilde{\rho}(r)$. While several specialized methods exist, two fundamental approaches are particularly relevant: the inverse CDF method and rejection sampling.

##### First compute the CDF:

$$F_r(r) \,=\, \int_0^r \mathcal{P}_r(r')\,dr' \,=\, \frac{\int_0^r \tilde{\rho}(r') \cdot r'^2\,dr'}{\int_0^{\infty} \tilde{\rho}(r'') \cdot r''^2\,dr''} \tag{2.22}$$

then, solve $F_r(r) = U_3$ for $r$. This requires an analytical or numerical inverse, which may be complex for arbitrary $\tilde{\rho}(r)$ [5, 7].

##### Rejection Sampling 

If the inverse CDF is intractable, use rejection sampling. First, choose a proposal distribution $q(r)$; e.g., a Gaussian or exponential, that is easy to sample from and bounds the target PDF. Next, sample $r \sim q(r)$ and $U \sim \text{Uniform}[0, 1]$. Then, accept $r$ if,

$$U \,\leq\, \frac{\mathcal{P}_r(r)}{\beta \cdot q(r)} \tag{2.23}$$

where $\beta$ is a constant such that $\mathcal{P}_r(r) \leq \beta \cdot q(r)$ for all $r$ [6, 7].

###### Example for a Gaussian Dispersion

For a Gaussian volumetric mass density, such as the one proposed in Modjtahedzadeh (2025) [1],

$$\rho(r) \,=\, \rho_0 \exp\left[-0.5 \left(\frac{r - \mu_R}{\sigma_R}\right)^2\right] \tag{2.24}$$

rejection sampling is often used. Rejection sampling is optimal for this radial PDF due to the mathematical intractability of its CDF inverse [5]. The product structure of the geometric factor ($r^2$) and Gaussian term creates a distribution ideally suited for rejection techniques [7], avoiding computationally expensive numerical integration while maintaining accuracy [6]. This approach is particularly efficient for three-dimensional physical simulations [4], and the bounded nature of the cloud's radius facilitates the construction of appropriate envelope functions [7], making rejection sampling both mathematically sound and computationally practical for debris cloud modeling. The radial PDF subsequently becomes:

$$\mathcal{P}_r(r) \propto r^2 \exp\left[-0.5 \left(\frac{r - \mu_R}{\sigma_R}\right)^2\right] \tag{2.25}$$

In computational implementation, two approaches can be used based on required accuracy and performance constraints. A pure rejection sampling implementation requires defining an envelope distribution $q(r)$ from which samples are easier to draw, typically a Gaussian centered at $\mu_R$ with standard deviation $\sigma_R$. The algorithm then iteratively samples $r_{\text{proposed}} \sim q(r)$ and accepts this sample with acceptance probability:

$$\frac{r^2_{\text{proposed}} \exp\left[-0.5 \left(\frac{r_{\text{proposed}} - \mu_R}{\sigma_R}\right)^2\right]}{\beta \cdot q(r_{\text{proposed}})} \tag{2.26}$$

where $\beta$ is chosen such that this ratio is always $\leq 1$. When an independently generated quantity $U$ is less than this ratio, accept the sample; otherwise, try again. This ensures that the accepted samples follow the desired target distribution.

For applications where computational efficiency is prioritized, rectification can be employed by directly sampling from a normal distribution with mean $\mu_R$ and standard deviation $\sigma_R$, then adjusting for non-negativity constraints by truncating negative values;

$$N \sim \text{Normal}(\mu_R, \sigma_R) \tag{2.27}$$

$$r = \max(0, N) \tag{2.28}$$

Rectification is a necessary physical correction as it will set $r \to 0$ if $r$ is computed to be negative. From a mathematical perspective, this creates a modified distribution that has the same probability density as the normal distribution for all positive values and has a point mass (spike) at exactly $r = 0$ that accumulates all the probability that would have been spread across negative values. Since the normal distribution has support on $(-\infty, +\infty)$, directly sampling $N$ produces a non-zero probability of generating negative values, which would be physically meaningless for radial distances. The likelihood of negative samples depends on the ratio $\mu_R/\sigma_R$. While a proper truncated normal distribution would be mathematically more precise [7], the rectified Gaussian approximation represents a computationally efficient approximation that preserves the essential statistical properties of the distribution while ensuring physical validity [8].

Anyhow, this approximation preserves the Gaussian profile while ensuring physically meaningful non-negative radial distances. The azimuthal angle is sampled uniformly just as in (2.14) and the polar angle is sampled with proper weighting by the solid angle element:

$$\cos\theta \,\sim\, \text{Uniform}[-1, 1] \tag{2.29}$$

$$\theta \,=\, \arccos(\cos\theta) \tag{2.30}$$

the complete position vector in Cartesian coordinates is then given from the standard coordinate transformation,

$$\vec{r} \,=\, (r\sin\theta\cos\varphi, r\sin\theta\sin\varphi, r\cos\theta). \tag{2.31}$$

This rectified Gaussian approximation method produces statistically accurate fragment distributions while maintaining reasonable simulation performance.

## References

[1] Modjtahedzadeh, K. (2025). Gaussian Distribution Model of Post-Impact Debris Cloud. *Manuscript in preparation*. Boeing Intelligence & Analytics.

[2] Ellgen, P. (2020). Probability Density Functions for Velocity Components in Spherical Coordinates. *Thermodynamics and Chemical Equilibrium*, LibreTexts Chemistry.

[3] Johnson, N. L., Krisko, P.H., Liou, J.-C., & Anz-Meador, P.D. (2001). NASA's New Breakup Model of EVOLVE 4.0. *Advances in Space Research*, 28(9), 1377–1384.

[4] Davis-Ross, K., & Chunn, J. (2021). An Introduction to Probability and Simulation. *Chapman and Hall/CRC*.

[5] Ross, S. M. (2019). Introduction to Probability Models (12th ed.). *Academic Press*.

[6] Law, A. M. (2015). Simulation Modeling and Analysis (5th ed.). *McGraw-Hill*.

[7] Devroye, L. (1986). Non-Uniform Random Variate Generation. *Springer-Verlag*.

[8] Gentle, J. E. (2003). Random Number Generation and Monte Carlo Methods (2nd ed.). *Springer*.

[^1]: It is called $\xi_{\text{inter}}$ instead of $\xi_{\text{middle}}$ because this function does not represent a unique statistical distribution for medium-sized fragments. Instead, it performs an interpolation (weighted blending) between the small and large fragment distributions. The name describes its mathematical operation rather than just the size category. <br><br> Without interpolation, a fragment at 7.9 cm would use entirely different statistical parameters than one at 8.1 cm, creating an unrealistic discontinuity in physical properties. Linear interpolation creates a gradual transition: At exactly 8 cm almost entirely small fragment distribution is used, at 9.5 cm (midpoint) weighting of both distributions is used, and at exactly 11 cm almost the entirety of large fragment distribution is used. This represents the physical reality that fragmentation processes produce a continuous spectrum of fragment properties rather than discrete categories with sharp boundaries.

[^2]: A marginal PDF is the probability density function of one or more variables from a joint distribution, obtained by averaging out (integrating) the other variables. In the context of sampling fragment positions in a 3D debris cloud, where positions are defined by spherical coordinates $(r, \theta, \varphi)$, the marginal PDF for the radial distance $r$ describes the likelihood of a fragment being at a certain distance from the origin, ignoring the angular coordinates $\theta$ and $\phi$. Its derived by summing up the probabilities over all possible angles, focusing only on the radial part of the distribution.

[^3]: $U \sim \text{Uniform}[0, 1]$ means that $U$ has equal probability of taking any value between 0 and 1

[^4]: In general, a uniform random variable $\mathcal{U}$ that's distributed on any interval $[c, d]$ can be transformed into a uniform distribution on $[a, b]$ using a modified formula: $$\mathcal{X} = a + (b - a) \frac{\mathcal{U} - c}{d - c}$$ This formula first normalizes $\mathcal{U}$ to $[0, 1]$ and then scales it to $[a, b]$ [5, 7].​​​​​​​​​​​​​​​​
