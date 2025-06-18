# Computational Implementation of Spherically Symmetric Debris Cloud Models

Kamyar Modjtahedzadeh  
Boeing Intelligence & Analytics  
June 9, 2025

> This paper presents a comprehensive computational framework for implementing spherically symmetric debris cloud models based on the NASA Standard Breakup Model. The implementation addresses three critical components of orbital debris simulation: statistical sampling of fragment area-to-mass ratios, three-dimensional spatial positioning according to Gaussian density distributions, and numerical fragment count determination across size categories. A piecewise statistical function is developed to model area-to-mass ratio distributions, incorporating unimodal log-normal sampling for small fragments and bimodal distributions for larger debris with smooth interpolation across the transition region. For spatial sampling, marginal and conditional probability density functions are derived in spherical coordinates, with specialized techniques including inverse cumulative distribution function methods and rejection sampling for Gaussian volumetric mass densities. The framework accounts for the steep power-law size distribution characteristic of fragmentation events, with adaptive resolution strategies that allocate computational resources efficiently across small, medium, and large fragment populations. Practical algorithms are provided for both high-accuracy rejection sampling and computationally efficient rectified Gaussian approximations, enabling flexible implementation based on performance requirements. This computational approach provides the foundation for simulating post-impact debris clouds in orbital environments, supporting risk assessment and collision avoidance applications in space situational awareness.

**Keywords:** Orbital debris, Fragmentation simulation, NASA EVOLVE 4.0, Statistical sampling

## 1 Distribution of Area-to-Mass Ratios

Based on the NASA breakup model's area-to-mass ratio distribution [3], a piecewise statistical function has been structured for the fragment's logarithm characteristic length ($L_{\mathrm{c}}$) input and its cross-sectional area-to-mass ratio ($A/M$) output. The function partitions the domain into three regions based on fragment size:

$$ \xi (x) = \begin{cases} \xi_{\mathrm{small}}(\phi) & \text{if } \phi < \log(0.08) \ \xi_{\mathrm{inter}}(\phi) & \text{if } \log(0.08) \leq \phi \leq \log(0.11) \ \xi_{\mathrm{large}}(\phi) & \text{if } \phi > \log(0.11) \end{cases} \tag{1.1} $$

where $\phi \equiv \log L_{\mathrm{c}}$, $\xi_{\mathrm{small}}$ is the small fragment distribution and $\xi_{\mathrm{large}}$ is the large fragment distribution.

### 1.1 Interpolation of Middle Region

For the middle region (8 cm ≤ $L_{\mathrm{c}}$ ≤ 11 cm), the $\xi$ function uses linear interpolation;

$$\xi_{\mathrm{inter}}(\phi) = \xi_{\mathrm{small}}(\phi) \cdot (1-w(\phi)) + \xi_{\mathrm{large}}(\phi) \cdot w(\phi) \tag{1.2}$$

where the weight function is given by,

$$w(\phi) = \frac{\phi - \log(0.08)}{\log(0.11) - \log(0.08)} \tag{1.3}$$

This creates a smooth transition between the statistical distributions for small and large fragments [3].[^1]

### 1.2 Statistical Sampling Functions

Each of the base functions $\xi_{\mathrm{small}}$ and $\xi_{\mathrm{large}}$ perform statistical sampling from normal distributions in logarithmic space. The output value is ultimately transformed from logarithmic space back to linear space. This mathematical structure creates a continuous statistical model for $A/M$ values that accommodates different physical regimes across the size spectrum of orbital debris fragments [3].

#### 1.2.1 Small Fragment Distribution

The small fragment distribution function $\xi_{\mathrm{small}}(\phi)$ implements a unimodal log-normal sampling model:

$$\xi_{\mathrm{small}}(\phi) = 10^{N_{\chi}} \tag{1.4}$$

where $N_{\chi} \sim \text{Normal}(\bar{\chi}, \sigma_{\chi})$ represents a random sample from a normal distribution, $\bar{\chi}$ is the mean of $\chi \equiv \log \frac{A}{M}$, and $\sigma_{\chi}$ is the standard deviation of $\chi$. The aforementioned normal distribution is:

$$\frac{1}{\sigma_{\chi}\sqrt{2\pi}} \exp\left[- \frac{(\chi - \bar{\chi})^2}{2\sigma_{\chi}^2}\right] \tag{1.5}$$

in Johnson et al. (2001), this normal distribution models the logarithm of area-to-mass ratios [3]. In other words, $\log(A/M)$ peaks and is symmetric around that peak value of $\bar{\chi}$.

Moreover, $\bar{\chi}$ and $\sigma_{\chi}$ follow piecewise functions. The expression for $\bar{\chi}(\phi)$ is:

$$ \bar{\chi}(\phi) = \begin{cases} -0.3 & \text{if } \phi \leq -1.75 \ -0.3 - 1.4(\phi + 1.75) & \text{if } -1.75 < \phi < -1.25 \ -1.0 & \text{if } \phi \geq -1.25, \end{cases} \tag{1.6} $$

and for the standard deviation of $\chi$:

$$ \sigma_{\chi}(\phi) = \begin{cases} 0.2 & \text{if } \phi \leq -3.5 \ 0.2 + 0.1333(\phi + 3.5) & \text{if } \phi > -3.5. \end{cases} \tag{1.7} $$

#### 1.2.2 Large Fragment Distribution

The large fragment distribution implements a bimodal log-normal sampling model [3];

$$ \xi_{\mathrm{large}}(\phi) = \begin{cases} 10^{N_{\chi_1}} & \text{with probability } \omega(\phi) \ 10^{N_{\chi_2}} & \text{with probability } 1 - \omega(\phi) \end{cases} \tag{1.8} $$

where the weighting parameter, $\omega$, determines how the bimodal distribution splits fragments between two different normal distributions. The term bimodal captures the fact that larger fragments tend to come from two different populations with distinct physical properties: Fragments with higher and lower $A/M$ values. Particles with higher area-to-mass ratios—represented by $10^{N_{\chi_1}}$—are typically less dense and from planar material components while those with lower ratios—$10^{N_{\chi_2}}$—usually more dense and from compact material. This approach is consistent with the physical reality that parent bodies are composed of materials with varying densities and structural properties.

All parameters are piecewise functions of $\phi$. For upper stage fragments, the weight parameter:

$$ \omega(\phi) = \begin{cases} 1.0 & \text{if } \phi \leq -1.4 \ 1.0 - 0.3571(\phi + 1.4) & \text{if } -1.4 < \phi < 0 \ 0.5 & \text{if } \phi \geq 0 \end{cases} \tag{1.9} $$

the first distribution mean:

$$ \bar{\chi}_1(\phi) = \begin{cases} -0.45 & \text{if } \phi \leq -0.5 \ -0.45 - 0.9 \cdot (\phi + 0.5) & \text{if } -0.5 < \phi < 0 \ -0.9 & \text{if } \phi \geq 0, \end{cases} \tag{1.10} $$

the first distribution standard deviation:

$$\sigma_{\chi_1}(\phi) = 0.55 \tag{1.11}$$

the second distribution mean:

$$\bar{\chi}_2(\phi) = -0.9 \tag{1.12}$$

and the second distribution standard deviation:

$$ \sigma_{\chi_2}(\phi) = \begin{cases} 0.28 & \text{if } \phi \leq -1.0 \ 0.28 - 0.1636 \cdot (\phi + 1) & \text{if } -1.0 < \phi < 0.1 \ 0.1 & \text{if } \phi \geq 0.1 \end{cases} \tag{1.13} $$

## 2 Sampling in 3D Space with a Spherically Symmetric Density

Given the spherically symmetric volumetric mass density ($\tilde{\rho}$) of a debris cloud, fragment positions can be generated that must follow that distribution. The statistical sampling here requires consideration of the spherical geometry to accurately represent the spatial distribution of fragments.

### 2.1 Spherical Volume Element in Probability Sampling

The spherical volume element is:

$$dV = r^2 \sin\theta, dr,d\theta,d\varphi \tag{2.1}$$

The normalized probability density function (PDF)—i.e., probability per unit volume—is [1]:

$$
\begin{align}
    \mathcal{P}(r) \,=\, \frac{\tilde{\rho}(r)}{M_{\mathrm{total}}}  \ ,
\end{align}
$$

where $M_{\mathrm{total}}$ is the total mass of all the fragments. The probability of finding a fragment in a specific volume element is then proportional to $\tilde{\rho},dV$;

$$\mathbb{P}_{r,\theta,\varphi} \propto \tilde{\rho}(r) \cdot r^2 \sin\theta \tag{2.3}$$

The statistical approach presented here reflects techniques developed in the NASA EVOLVE 4.0 framework for modeling fragmentation debris clouds in orbital environments [3].

### 2.2 Marginal and Conditional PDFs

To sample a three-dimensional point in spherical coordinates, the marginal PDFs for the radial and angular components.[^2]

#### 2.2.1 Radial PDF

Following (2.3), normalizing the radial probability—i.e., the radial component of $\mathbb{P}_{r,\theta,\varphi}$—requires dividing it by the radial segment of the total mass [1, 2];

$$\int_0^{\infty} \mathbb{P}_r(r),dr = \int_0^{\infty} \frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2, dr'},dr = 1 \tag{2.4}$$

The marginal PDF for $r$ is obtained by integrating $\mathbb{P}$ over $\theta$ and $\varphi$:

$$\mathcal{P}_r(r) = \int_0^{2\pi} \int_0^{\pi} \mathbb{P}_{r,\theta,\varphi}, d\theta, d\varphi \tag{2.5}$$

$$= \int_0^{2\pi} \int_0^{\pi} \frac{\tilde{\rho}(r) \cdot r^2 \sin\theta}{4\pi \int_0^{\infty} \tilde{\rho}(r') \cdot r'^2, dr'}, d\theta, d\varphi \tag{2.6}$$

$$= \frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2, dr'} \tag{2.7}$$

The radial PDF, which shows that the probability density for the radial distance is proportional to $\tilde{\rho} \times r^2$, reflecting the increasing volume of spherical shells at larger $r$.

#### 2.2.2 Angular PDFs

The joint angular PDF, conditioned on $r$, is:

$$\mathcal{P}_{\theta,\varphi}(r) = \frac{\mathbb{P}_{r,\theta,\varphi}}{\mathbb{P}_r(r)} \tag{2.8}$$

$$= \left(\frac{\tilde{\rho}(r) \cdot r^2 \sin\theta}{4\pi \int_0^{\infty} \tilde{\rho}(r') \cdot r'^2, dr'}\right) \bigg/ \left(\frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2, dr'}\right) \tag{2.9}$$

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

$$\varphi = 2\pi U_1 \tag{2.16}$$

#### 2.3.2 Polar Angle

The polar angle's PDF is given from (2.12). To sample $\theta$, use the cumulative distribution function (CDF):

$$F_{\theta}(\theta) = \int_0^{\theta} \mathcal{P}_{\theta}(\theta'),d\theta' = \frac{1 - \cos\theta}{2}, \tag{2.17}$$

then set $F_{\theta}(\theta) = U_2$ where $U_2 \sim \text{Uniform}[0, 1]$,

$$\theta = \arccos(1 - 2U_2). \tag{2.18}$$

Alternatively, since $\cos\theta \in [-1, 1]$ and $U_2$ follows a uniform distribution on $[0, 1]$, a standard linear transformation for converting a uniform random variable from one range to another can be used. The formula to be used is:

$$X = a + (b - a)U, \tag{2.19}$$

where $X$ is uniformly distributed on $[a, b]$ [4, 6].[^4] Since $a$ and $b$ are respectively equal to $-1$ and $+1$, $X \to \cos\theta = 2U_2 - 1$; and so,

$$\theta = \arccos(2U_2 - 1) \tag{2.20}$$

This ensures the angular distribution is isotropic over the sphere.

#### 2.3.3 Radial Distance

The PDF for $r$ is:

$$\mathcal{P}_r(r) = \frac{\tilde{\rho}(r) \cdot r^2}{\int_0^{\infty} \tilde{\rho}(r') \cdot r'^2,dr'} \text{ for } r \in [0, \infty) \tag{2.21}$$

sampling from this PDF depends on the form of $\tilde{\rho}(r)$. While several specialized methods exist, two fundamental approaches are particularly relevant: the inverse CDF method and rejection sampling.

##### Inverse CDF Method

First compute the CDF:

$$F_r(r) = \int_0^r \mathcal{P}_r(r'),dr' = \frac{\int_0^r \tilde{\rho}(r') \cdot r'^2,dr'}{\int_0^{\infty} \tilde{\rho}(r'') \cdot r''^2,dr''} \tag{2.22}$$

then, solve $F_r(r) = U_3$ for $r$. This requires an analytical or numerical inverse, which may be complex for arbitrary $\tilde{\rho}(r)$ [5, 7].

##### Rejection Sampling

If the inverse CDF is intractable, use rejection sampling. First, choose a proposal distribution $q(r)$; e.g., a Gaussian or exponential, that is easy to sample from and bounds the target PDF. Next, sample $r \sim q(r)$ and $U \sim \text{Uniform}[0, 1]$. Then, accept $r$ if,

$$U \leq \frac{\mathcal{P}_r(r)}{\beta \cdot q(r)} \tag{2.23}$$

where $\beta$ is a constant such that $\mathcal{P}_r(r) \leq \beta \cdot q(r)$ for all $r$ [6, 7].

**Example for a Gaussian Dispersion** For a Gaussian volumetric mass density, such as the one proposed in Modjtahedzadeh (2025) [1],

$$\rho(r) = \rho_0 \exp\left[-0.5 \left(\frac{r - \mu_{\mathrm{R}}}{\sigma_{\mathrm{R}}}\right)^2\right] \tag{2.24}$$

rejection sampling is often used. Rejection sampling is optimal for this radial PDF due to the mathematical intractability of its CDF inverse [5]. The product structure of the geometric factor ($r^2$) and Gaussian term creates a distribution ideally suited for rejection techniques [7], avoiding computationally expensive numerical integration while maintaining accuracy [6]. This approach is particularly efficient for three-dimensional physical simulations [4], and the bounded nature of the cloud's radius facilitates the construction of appropriate envelope functions [7], making rejection sampling both mathematically sound and computationally practical for debris cloud modeling. The radial PDF subsequently becomes:

$$\mathcal{P}_r(r) \propto r^2 \exp\left[-0.5 \left(\frac{r - \mu_{\mathrm{R}}}{\sigma_{\mathrm{R}}}\right)^2\right] \tag{2.25}$$

In computational implementation, two approaches can be used based on required accuracy and performance constraints. A pure rejection sampling implementation requires defining an envelope distribution $q(r)$ from which samples are easier to draw, typically a Gaussian centered at $\mu_{\mathrm{R}}$ with standard deviation $\sigma_{\mathrm{R}}$. The algorithm then iteratively samples $r_{\mathrm{proposed}} \sim q(r)$ and accepts this sample with acceptance probability:

$$\frac{r^2_{\mathrm{proposed}} \exp\left[-0.5 \left(\frac{r_{\mathrm{proposed}} - \mu_{\mathrm{R}}}{\sigma_{\mathrm{R}}}\right)^2\right]}{\beta \cdot q(r_{\mathrm{proposed}})} \tag{2.26}$$

where $\beta$ is chosen such that this ratio is always $\leq 1$. When an independently generated quantity $U$ is less than this ratio, accept the sample; otherwise, try again. This ensures that the accepted samples follow the desired target distribution.

For applications where computational efficiency is prioritized, rectification can be employed by directly sampling from a normal distribution with mean $\mu_{\mathrm{R}}$ and standard deviation $\sigma_{\mathrm{R}}$, then adjusting for non-negativity constraints by truncating negative values;

$$N_r \sim \text{Normal}(\mu_{\mathrm{R}}, \sigma_{\mathrm{R}}) \tag{2.27}$$

$$r = \max(0, N_r) \tag{2.28}$$

Rectification is a necessary physical correction as it will set $r \to 0$ if $r$ is computed to be negative. From a mathematical perspective, this creates a modified distribution that has the same probability density as the normal distribution for all positive values and has a point mass (spike) at exactly $r = 0$ that accumulates all the probability that would have been spread across negative values. Since the normal distribution has support on $(-\infty, +\infty)$, directly sampling $N_r$ produces a non-zero probability of generating negative values, which would be physically meaningless for radial distances. The likelihood of negative samples depends on the ratio $\mu_{\mathrm{R}}/\sigma_{\mathrm{R}}$. While a proper truncated normal distribution would be mathematically more precise [7], the rectified Gaussian approximation represents a computationally efficient approximation that preserves the essential statistical properties of the distribution while ensuring physical validity [8].

This approximation preserves the Gaussian profile while ensuring physically meaningful non-negative radial distances. The azimuthal angle is sampled uniformly just as in (2.14) and the polar angle is sampled with proper weighting by the solid angle element:

$$\cos\theta \sim \text{Uniform}[-1, 1] \tag{2.29}$$

$$\theta = \arccos(\cos\theta) \tag{2.30}$$

the complete position vector in Cartesian coordinates is then given from the standard coordinate transformation,

$$\vec{r} = (r\sin\theta\cos\varphi, r\sin\theta\sin\varphi, r\cos\theta). \tag{2.31}$$

This rectified Gaussian approximation method produces statistically accurate fragment distributions while maintaining reasonable simulation performance.

## 3 Numerical Size Distribution

### 3.1 Counts by Size Category

To determine the number of fragments within specific size ranges, use the cumulative distribution function ($N$). The number of fragments between characteristic lengths $L_1$ and $L_2$, where $L_1 < L_2$, is $N(L_1) - N(L_2)$.

The practical bounds for all size categories are set to $L_{\mathrm{min}} = 0.001$ m and $L_{\mathrm{max}} = 1.0$ m. Using the cumulative distribution function for collisions,

$$\Upsilon_{\mathrm{small}} = 13482.1M_{\mathrm{parent}}^{0.75} \tag{3.1}$$

$$\Upsilon_{\mathrm{medium}} = 3.15M_{\mathrm{parent}}^{0.75} \tag{3.2}$$

$$\Upsilon_{\mathrm{large}} = 4.26M_{\mathrm{parent}}^{0.75} \tag{3.3}$$

where $\Upsilon$ is the symbol for number of.

#### 3.1.1 Percentage Per Category (Collision)

The total number of fragments is:

$$\Upsilon_{\mathrm{total}} = \Upsilon_{\mathrm{small}} + \Upsilon_{\mathrm{medium}} + \Upsilon_{\mathrm{large}} = 13489.5M_{\mathrm{parent}}^{0.75} \tag{3.4}$$

The percentage distributions are: 
&nbsp;&nbsp;&nbsp;&nbsp;• Small fragments: $99.95\%$  
&nbsp;&nbsp;&nbsp;&nbsp;% • Medium fragments: $0.023\%$ 
&nbsp;&nbsp;&nbsp;&nbsp;% • Large fragments: $0.032\%$  %

#### 3.1.2 Intra-Category Distribution

It has been established how the fragments are distributed based on their classification; however, the three size categories also have internal distributions. The small, medium, and large fragments themselves follow the same respective power law distribution pattern. Within any subrange size category, the fragments are not uniformly distributed as smaller fragments within any classification are more numerous than larger ones. For any subrange within any category, from $_{\mathrm{a}}$ to $_{\mathrm{b}}$, where $0.001$ m $\leq L_{\mathrm{a}} < L_{\mathrm{b}} \leq 1.0$ m, the number of fragments is:

$$\Upsilon_{\mathrm{subrange}} = N(L_{\mathrm{a}}) - N(L_{\mathrm{b}}) = 0.1M_{\mathrm{parent}}^{0.75} \left( L_{\mathrm{a}}^{-1.71} - L_{\mathrm{b}}^{-1.71} \right) \tag{3.5}$$

### 3.2 Discrete Fragment Count

#### 3.2.1 Size Step Configuration

The power law exponent ($-1.71$) means fragment counts change much more rapidly at smaller sizes. Small fragments need fine resolution to capture this steep variation, while longer fragments change more gradually. Allocate computational resources where matter most; high resolution in where counts change rapidly (small), moderate resolution in the transition zone (medium), and standard resolution where changes are gradual (large).

#### 3.2.2 Interval Sampling

To estimate the number of fragments for any particular characteristic length, the following formula is used:

$$\Upsilon(L_{\mathrm{c}}) = N(L_{\mathrm{c}}) - N(L_{\mathrm{c}} + \varepsilon), \tag{3.6}$$

where $\varepsilon$ is an infinitesimal increment. With this formulation, $\Upsilon(L_{\mathrm{c}}) \approx \varepsilon \cdot _{\mathrm{c}})$. Loop through the discrete $L_{\mathrm{c}}$ values with the assigned $\Delta L_{\mathrm{c}}$ size steps and calculate the number of pieces for each characteristic length with the above formula for $\Upsilon(L_{\mathrm{c}})$.

## 4 Impact Probability at Event Zero

Given the parent mass, initial cloud radius, cumulative distribution function, $\rho(r, L_{\mathrm{c}})$, and hit distance threshold, the impact probability will be calculated such that a linear trajectory passing through random entry and exit points on the cloud sphere at time $t = 0$ will encounter at least one fragment within distance $\ell$. A fundamental result in stochastic processes from Poisson process theory is the application of a non-homogeneous Poisson process [12]; since encounters occur randomly at rate (per unit length) $\Lambda(\mathcal{L})$ along path $\mathcal{L}$, then the number of encounters follows a Poisson distribution with parameter $\int\Lambda d\mathcal{L}'$. The probability of at least one encounter is $1-\mathbb{P}P_{\mathrm{avoid}} = 1-\exp[-\int\Lambda d\mathcal{L}']$ [10, 12]. The total impact probability is computed as the weighted sum over all fragment size categories:

$$\mathbb{P}P_{\mathrm{impact}} = \sum_{L_{\mathrm{c}} =L_{\mathrm{min}}}^{L_{\mathrm{max}}} f \cdot \mathbb{P}P_{\mathrm{impact}}(L_{\mathrm{c}}) , \tag{4.1}$$

where $f$ represents the fraction of fragments given from percentage distributions per size category, as in Section 3.1.1.

For a specific characteristic length, the probability of impact along a trajectory is [11]:

$$\mathbb{P}_{\mathrm{impact}}(L_{\mathrm{c}}, \mathcal{L}) = 1-\exp\left[-\int_L \Lambda(L_{\mathrm{c}}, \mathcal{L}'),d\mathcal{L}'\right] . \tag{4.2}$$

The collision rate at position $r$ along trajectory $\mathcal{L}$ is [9, 13–16]:

$$\Lambda(r, \mathcal{L}) = \rho_N(r, \mathcal{L})\cdot \pi\ell^2 , \tag{4.3}$$

$\rho_N(r)$ is the number density of fragments at radial distance $r$ and $\pi\ell^2$ represents the effective impact cross-sectional area. The number density is related to the mass density through,

$$\rho_N(r, L_{\mathrm{c}}) = \frac{\rho(r, L_{\mathrm{c}})}{M} \frac{dN}{d_{\mathrm{c}}} , \tag{4.4}$$

where $M$ is the average mass of fragments with characteristic length $_{\mathrm{c}}$.

### 4.1 Monte Carlo Implementation

#### 4.1.1 Estimator Design and Variance Analysis

The impact probability is estimated through repeated sampling [4–6]:

$$\hat{\mathbb{P}P}_{\mathrm{impact}} = \frac{1}{\Upsilon_{\mathrm{trials}}} \sum_{j=1}^{\Upsilon_{\mathrm{trials}}} I_j, \tag{4.5}$$

where $I_j$ is an indicator function:

$$I_j = \begin{cases} 1 & \text{if trajectory } j \text{ encounters any fragment} \ 0 & \text{otherwise.} \end{cases} \tag{4.6}$$

This indicator function captures whether a given trajectory $j$ results in a collision, effectively turning each trial into a binary outcome. The reliability of this approximation depends on the number of trials and the inherent variability of the sampling process, which is quantified through the variance of the estimator.

The variance of the Monte Carlo estimator is [5, 6]:

$$\text{Var}(\mathbb{P}) = \frac{\mathbb{P}(1-\mathbb{P})}{\Upsilon_{\mathrm{trials}}} . \tag{4.7}$$

The confidence interval using the Wilson score method is [17]:

$$\mathbb{P} \in \left[\frac{\hat{\mathbb{P}P} + \frac{\zeta^2}{2\Upsilon} - \zeta\sqrt{\frac{\hat{\mathbb{P}}(1-\hat{\mathbb{P}P})}{\Upsilon} + \frac{\zeta^2}{4\Upsilon^2}}}{1+\frac{\zeta^2}{\Upsilon}}, \frac{\hat{\mathbb{P}P} + \frac{\zeta^2}{2\Upsilon} + \zeta\sqrt{\frac{\hat{\mathbb{P}}(1-\hat{\mathbb{})}{\Upsilon} + \frac{\zeta^2}{4\Upsilon^2}}}{1+\frac{\zeta^2}{\Upsilon}}\right] , \tag{4.8}$$

where $\zeta$ is the z-score (standard normal quantile) [4–6].

#### 4.1.2 Adaptive Sample Size with Sequential Refinement

Since $\mathbb{P}$ is unknown a priori, an adaptive sampling approach is employed [5, 6]. Starting with an initial batch of $\Upsilon_0 = 10^4$ trials, the required sample size is updated after each batch;

$$\Upsilon_{\mathrm{next}} = \frac{\zeta^2(1-\hat{\mathbb{P}P}_{\mathrm{current}})}{\hat{\mathbb{P}P}_{\mathrm{current}}\epsilon^2} . \tag{4.9}$$

where $\epsilon$ is the tolerance. The sampling continues until convergence is achieved [4–6]:

$$\frac{\mathbb{P}P_{\mathrm{upper}} - \mathbb{P}P_{\mathrm{lower}}}{\hat{\mathbb{P}P}_{\mathrm{current}}} < \epsilon , \tag{4.10}$$

where $\mathbb{P}P_{\mathrm{lower}}$ and $\mathbb{P}P_{\mathrm{upper}}$ come from the confidence interval formula in (4.8). As more trials are added, the $\sqrt{\Upsilon}$ in the denominator of the Wilson formula makes the interval tighter, reducing the relative width until it's less than the tolerance [4–6].

### 4.2 Importance Sampling

Instead of uniform sampling on the sphere, use importance sampling where samples are drawn from a biased distribution (called the importance distribution) that prioritizes "important" regions of the sample space—those more likely to contribute to the quantity being estimated, and then adjust for this bias by weighting the samples [5, 6, 16]. To sample entry and exit points this way, use the following non-normalized PDF:



$$\mathcal{P}(p_{\mathrm{entry}}, p_{\mathrm{exit}}) \propto \exp\left[-\frac{\ell_{\mathrm{min}}^2(p_{\mathrm{entry}}, p_{\mathrm{exit}})}{2\sigma_{\mathrm{IS}}^2}\right] , \tag{4.11}$$

where $\ell_{\mathrm{min}}$ is a function that calculates the minimum distance between the trajectory defined by the entry and exit points and the peak density sphere at radius $\mu R_{\mathrm{c}}$ and $\sigma_{\mathrm{IS}}$ is the importance sampling standard deviation, a tuning parameter that controls how strongly the sampling is biased toward trajectories passing near the peak density radius [5, 6].

Furthermore, the importance-sampled estimator becomes,

$$\hat{P}_{\mathrm{IS}} = \frac{1}{\Upsilon_{\mathrm{trials}}} \sum_{j=1}^{\Upsilon_{\mathrm{trials}}} I_j w_j, \tag{4.12}$$

where the weight is [5, 6]:

$$w_j = \frac{h_j(p_{\mathrm{entry}}, p_{\mathrm{exit}})}{\mathcal{P}_j(p_{\mathrm{entry}}, p_{\mathrm{exit}})} , \tag{4.13}$$

with $h$ representing the uniform distribution on the sphere. Without weights, estimates would be biased high because trajectories near the peak (where fragments are dense) are oversampled. The weights correct this by giving oversampled trajectories proportionally lower contributions [5,18]. This approach works because it gives an unbiased estimate of the true probability, even though sampling is done non-uniformly [5, 6, 18].[^5]

## References

[1]  Modjtahedzadeh, K. (2025). _Gaussian distribution model of post-impact debris cloud_. Manuscript in preparation. Boeing Intelligence & Analytics.  
  
[2]  Ellgen, P. (2020). _Probability density functions for velocity components in spherical coordinates_. In _Thermodynamics and chemical equilibrium_. LibreTexts Chemistry.  
  
[3]  Johnson, N. L., Krisko, P. H., Liou, J.-C., & Anz-Meador, P. D. (2001). NASA's new breakup model of EVOLVE 4.0. _Advances in Space Research,_ **28**(9), 1377–1384.  
  
[4]  Davis-Ross, K., & Chunn, J. (2021). _An introduction to probability and simulation_. Chapman and Hall/CRC.  
  
[5]  Ross, S. M. (2019). _Introduction to probability models_ (12th ed.). Academic Press.  
  
[6]  Law, A. M. (2015). _Simulation modeling and analysis_ (5th ed.). McGraw-Hill.  
  
[7]  Devroye, L. (1986). _Non-uniform random variate generation_. Springer-Verlag.  
  
[8]  Gentle, J. E. (2003). _Random number generation and Monte Carlo methods_ (2nd ed.). Springer.  
  
[9]  Chapman, S., & Cowling, T. G. (1970). _The mathematical theory of non-uniform gases_ (3rd ed.). Cambridge University Press.  
  
[10] Duderstadt, J. J., & Hamilton, L. J. (1976). _Nuclear reactor analysis_. John Wiley & Sons.  
  
[11] Letizia, F., Colombo, C., & Lewis, H. G. (2016). Collision probability due to space debris clouds through a continuum approach. _Journal of Guidance, Control, and Dynamics,_ **39**(10), 2240–2249.  
  
[12] Haight, F. A. (1967). _Handbook of the Poisson distribution_. John Wiley & Sons.  
  
[13] McQuarrie, D. A. (2000). _Statistical mechanics_. University Science Books.  
  
[14] Perkins, D. H. (2000). _Introduction to high energy physics_ (4th ed.). Cambridge University Press.  
  
[15] Grün, E., Zook, H. A., Fechtig, H., & Giese, R. H. (1985). Collisional balance of the meteoritic complex. _Icarus,_ **62**(2), 244–272.  
  
[16] Bird, G. A. (1994). _Molecular gas dynamics and the direct simulation of gas flows_. Oxford University Press.  
  
 [17] Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. _Journal of the American Statistical Association,_ **22**(158), 209–212.  
  
[18] Kloek, T., & van Dijk, H. K. (1978). Efficient estimation of income distribution parameters. _Journal of Econometrics,_ **8**(1), 61–74.  


[^1]: It is called $\xi_{\mathrm{inter}}$ instead of $\xi_{\mathrm{middle}}$ because this function does not represent a unique statistical distribution for medium-sized fragments. Instead, it performs an interpolation (weighted blending) between the small and large fragment distributions. The name describes its mathematical operation rather than just the size category. Without interpolation, a fragment at 7.9 cm would use entirely different statistical parameters than one at 8.1 cm, creating an unrealistic discontinuity in physical properties. Linear interpolation creates a gradual transition: At exactly 8 cm almost entirely small fragment distribution is used, at 9.5 cm (midpoint) weighting of both distributions is used, and at exactly 11 cm almost the entirety of large fragment distribution is used. This represents the physical reality that fragmentation processes produce a continuous spectrum of fragment properties rather than discrete categories with sharp boundaries.

[^2]: A marginal PDF is the probability density function of one or more variables from a joint distribution, obtained by averaging out (integrating) the other variables. In the context of sampling fragment positions in a 3D debris cloud, where positions are defined by spherical coordinates $(r, \theta, \varphi)$, the marginal PDF for the radial distance $r$ describes the likelihood of a fragment being at a certain distance from the origin, ignoring the angular coordinates $\theta$ and $\varphi$. Its derived by summing up the probabilities over all possible angles, focusing only on the radial part of the distribution.

[^3]: $\mathcal{U} \sim \text{Uniform}[0, 1]$ means that $U$ has equal probability of taking any value between 0 and 1.

[^4]: In general, a uniform random variable $\mathcal{U}$ that's distributed on any interval $[c, d]$ can be transformed into a uniform distribution on $[a, b]$ using a modified formula: $\mathcal{X} = a + (b - a) \frac{\mathcal{U} - c}{d - c}$ This formula first normalizes $\mathcal{U}$ to $[0, 1]$ and then scales it to $[a, b]$ [5, 7].

[^5]: The weights could also be normalized first such that the normalized weight is equal to $w_j/\sum_{i=1}^{\Upsilon} w_i$; however, the non-normalized weight is numerically more stable [5, 6].
