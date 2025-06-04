# Gaussian Distribution Model of Post-Impact Debris Cloud

Kamyar Modjtahedzadeh  
Boeing Intelligence & Analytics  
May 5, 2025

> A comprehensive model describing the evolution of satellite debris clouds following fragmentation events in orbit is presented. The approach models the spatial density distribution of fragments using a Gaussian function that incorporates the physical characteristics of debris particles. The model accounts for size-dependent stratification through an inverse relationship between the applied perturbational forces and characteristic length. These prominent non-gravitational forces cause the initial spherical cloud with differently sized fragments to evolve along distinct trajectories, transforming the geometry of the original structure. It follows that this representation provides the groundwork for a computational model to simulate and predict the system for short time intervals post-fragmentation, i.e., $\lesssim{120}$ seconds. 

**Keywords:** Space debris, Density distribution, Cloud propagation, Orbital mechanics

## 1 Basic Modeling of PID Cloud Radius using Packing Density

### 1.1 Fundamental Model of Post-Impact Fragmentation

A uniform solid sphere with radius $R_\mathrm{s}$ exists in space when it either explodes or collides with another body such as an interceptor and causes the sphere to fragment into equivalent pieces of post-impact debris (PID) with all the PID fragments being evenly spaced. The PID forms a cloud resembling a sphere with radius $R_\mathrm{c}$.

To model the relationship between the radius of the cloud sphere and the original solid sphere, the problem will be analyzed using principles of mass conservation and packing density. Since mass is conserved and the original sphere was uniform, the total volume of all PID pieces must equal the volume of the original sphere; however, the PID pieces are now distributed throughout a larger space with voids between them.

### 1.2 Packing Density

While the actual material volume remains the same, the PID particles no longer complete their respective sphere, creating empty space between fragments. The packing density, $\eta$, represents the fraction of the volume that is actually occupied by material (with the rest being empty space).[^1]

For the solid sphere, the packing fraction $\eta_\mathrm{s}$ is virtually 1 since it is uniform and fully packed. By conservation of mass: $\eta_\mathrm{s} V_\mathrm{s} = \eta_\mathrm{c} V_\mathrm{c}$, and since $\eta_\mathrm{s} = 1$, $V_\mathrm{s} = \eta_\mathrm{c} V_\mathrm{c}$. Substituting the formulas for spherical volumes: $\frac{4}{3}\pi R_\mathrm{s}^3 = \eta_\mathrm{c} \cdot \frac{4}{3}\pi R_\mathrm{c}^3$.

Solving for the ratio of radii:

$$ \frac{R_\mathrm{c}}{R_\mathrm{s}} = \left(\frac{1}{\eta_\mathrm{c}}\right)^{1/3} \tag{1.1} $$

This formula shows that the ratio of the cloud radius to the solid sphere radius depends on the cube root of the reciprocal of the packing density. Since $\eta_\mathrm{c} < 1$, $R_\mathrm{s} < R_\mathrm{c}$ as expected.

## 2 Spacial Dispersion and Fragment Size

### 2.1 Proposed Model of Fragment Normal Distribution

Realistically, the distribution of particles is not uniform throughout the cloud. Now that the basic relationship between $R_\mathrm{c}$ and $R_\mathrm{s}$ has been established in Equation (1.1), the model will be expanded such that the post-intercept debris particles effectively form a normal distribution.

#### 2.1.1 Radial Density Function

The volumetric mass density of the debris is modeled as a function of radial distance from the center using a Gaussian distribution,

$$ \rho(r) = \rho_0 \exp\left[-\frac{(r - \mu R_\mathrm{c})^2}{2\sigma^2 R_\mathrm{c}^2}\right] \tag{2.1} $$

where $\rho_0$ is a normalization constant, $\mu$ represents the peak density radial position,[^2] $\sigma$ is the standard deviation of the debris distribution, and $r$ is the radial distance from the center. Moreover, it follows that the mass-integral formula is given by:

$$ M_\mathrm{c} = 4\pi \int_{0}^{R_\mathrm{c}} \rho r^2 dr \tag{2.2} $$

where $M_\mathrm{c}$ is the total mass of the cloud sphere and hence the original solid sphere.[^3]

#### 2.1.2 Radial Packing Density Function

With this non-uniform distribution, a position-dependent packing density can be defined:

$$ \eta_\mathrm{c}(r) = \frac{\rho(r)}{\rho_{\mathrm{max}}} \cdot \eta_{\mathrm{max}} \tag{2.3} $$

$\rho_{\mathrm{max}}$ is the maximum value of $\rho(r)$ and $\eta_{\mathrm{max}}$ is the maximum packing density (which is at the peak of the distribution).

To put Equation (2.3) in literal terms: The effective packing density of the spherical cloud at radial position $r$ is equivalent to the packing density at the peak of the distribution × the ratio of the cloud's volume mass density at $r$ to its maximum volume mass density.

### 2.2 Size Dependency[^4]

The model presented in Sections 1 and 2.1 operated under the assumption that fragmentation produces equivalent particles. While this premise helped establish the Gaussian distribution framework, a single impact actually produces a wide spectrum of fragment sizes [1,3]. In this section, the model is expanded to account for heterogeneous particles within a single event.

To implement this expansion, Equation (2.1) is extended to incorporate the fragment characteristic length ($L_\mathrm{c}$), which is defined as the average of the three orthogonal dimensions of a fragment according to the NASA Standard Breakup Model [1]. The parameters that now depend on $L_\mathrm{c}$ are: the normalization constant $\rho_0(L_\mathrm{c})$, $\mu(L_\mathrm{c})$ which now defines how the peak density location varies with fragment size, and $\sigma(L_\mathrm{c})$ which is the distribution of fragments with characteristic length $L_\mathrm{c}$; consequently, $\rho(r) \rightarrow \rho(r, L_\mathrm{c})$. This formulation allows for mixed-sized particles created from the same event to display recognizably different spatial distributions within the PID cloud, consistent with observations [2, 6].

#### 2.2.1 Characteristic Length and Spacial Distribution Relationship

The NASA Standard Breakup Model links ejection velocity to $L_\mathrm{c}$, standard deviation can then be connected to characteristic length based on the average values of both explosion and collision velocity distributions [8]. Johnson et al. (2001) describe these values as:

$$ \bar{v}_{\mathrm{explosion}} = 0.2\log \frac{A}{M} + 1.85 \tag{2.4} $$

$$ \bar{v}_{\mathrm{collision}} = 0.9\log \frac{A}{M} + 2.9 \tag{2.5} $$

where $A$ is the cross-sectional area of the fragment and $M$ is its mass.[^5] Fragment velocity depends on $A/M$, with smaller fragments typically having higher ratios and therefore higher velocities [1]. Since spatial dispersion is proportional to velocity and velocity is inversely related to fragment size [3]:

$$ \sigma(L_\mathrm{c}) \simeq \sigma_0 L_\mathrm{c}^{-\alpha} \tag{2.6} $$

where $\sigma_0$ and $\alpha$ are empirical parameters.

#### 2.2.2 Mass Conservation and Fragment Size Distribution Functions

For the model to conserve the total mass of the original object, it must now account for the size of the fragments. The total mass (of all the PID particles) is now:

$$ M_{\mathrm{total}} = 4\pi \int_{L_{\mathrm{min}}}^{L_{\mathrm{max}}} \int_{0}^{R_\mathrm{c}} r^2 \rho(r, L_\mathrm{c})dr n(L_\mathrm{c})dL_\mathrm{c} \tag{2.7} $$

where $n(L_\mathrm{c})$ is the differential size distribution function. It is given from the cumulative size distribution function via,

$$ n(L_\mathrm{c}) = -\frac{d}{dL_\mathrm{c}}N(L_\mathrm{c}) \tag{2.8} $$

$N$ represents the number of particles with characteristic length larger than $L_\mathrm{c}$ whereas $n$ denotes the number of particles at exactly $L_\mathrm{c}$. For explosions and collisions respectively [1, 2]:

$$ N_{\mathrm{explosion}}(L_\mathrm{c}) = 6L_\mathrm{c}^{-1.6} \tag{2.9} $$

$$ N_{\mathrm{collision}}(L_\mathrm{c}) = 0.1L_\mathrm{c}^{-1.71} M_{\mathrm{parent}}^{0.75} \tag{2.10} $$

The parent object is the original intact object(s) before fragmentation, which is just the solid sphere in this model.

## 3 Temporal Evolution of the PID Cloud Density

### 3.1 Transition from Initial State

A static model of the PID cloud at the moment of impact has been set up to eventually model the system over time. With orbital mechanics, the cloud sphere with Gaussian density will be transforms into more complex shapes. For a debris cloud generated in orbit, this morphing emulates the below sequence [4, 7]:

$$
	\mathrm{Spherical} \,\xrightarrow{\hphantom{l}t \,\approx\, 0.25 T_{\mathrm{orb}}\hphantom{l}} \,\mathrm{Ellipsoidal} \,\xrightarrow{\hphantom{l}t \,\approx\, T_{\mathrm{orb}}\hphantom{l}} \,\mathrm{Toroidal}
$$

where $T_{\mathrm{orb}}$ represents the orbital period of the parent object.

As fragments follow their individual orbital paths, the variations in velocity vectors between fragments (differential velocities) cause the cloud to elongate along the orbital trajectory while largely maintaining its cohesion in directions perpendicular to the orbital path. After approximately one orbital period, the cloud forms a characteristic toroidal structure as fragments complete their first revolution around Earth [4, 7].

### 3.2 Time-Dependent Density Function for Short Timescales

To extend the model up to short time periods after impact ($\approx 120 \, \mathrm{s}$), $\rho(r, L_\mathrm{c})$ must be modified to account for the PID cloud's evolution with time. The parameters of the radial density are: $\rho_0(L_\mathrm{c})$, $\mu(L_\mathrm{c})$, $\sigma(L_\mathrm{c})$, and $R_\mathrm{c}$.

Making the cloud radius time-dependent does make sense because differential velocities cause the cloud to naturally expand over time due to the different velocities of fragments. In the initial moments after impact, expansion velocity is roughly constant [5]; for a short simulation, a simple linear expansion model should suffice,

$$ R_\mathrm{c}(t) = R_\mathrm{c}(0) + v_{\mathrm{expansion}} \cdot t \tag{3.1} $$

The standard deviation, on the other hand, is a dimensionless parameter that represents the spread of the distribution as a fraction of $R_\mathrm{c}$. Without perturbations, it will be assumed to remain constant during the early stages of the temporal evolution of the debris cloud; $\sigma(L_\mathrm{c}) \not\rightarrow \sigma(L_\mathrm{c}, t)$. As for the parameter representing the peak density position, for a short timescale simulation, keeping $\mu$ constant is a valid zeroth-order approximation as Vallado & Oltrogge (2017) examine the temporal evolution of debris clouds and establish the distinct evolutionary phases laid out in Section 3.1. Their analysis shows that the initial spherical phase lasts significantly longer than 120 seconds before transitioning to an ellipsoidal shape due to differential orbital mechanics; $\mu(L_\mathrm{c}) \not\rightarrow \mu(L_\mathrm{c}, t)$. Last but not least, the normalization constant must evolve to preserve total mass as the cloud expands.

$$ \rho_0(L_\mathrm{c}, t) = \rho_0(L_\mathrm{c}, 0)\left[\frac{R_\mathrm{c}(0)}{R_\mathrm{c}(t)}\right]^3 \tag{3.2} $$

this time-dependent normalization factor $\rho_0(t)$ ensures that the total mass of the parent object is conserved across time. To summarize, in Equation (2.1), $\rho_0$ and $\mu$ are still not functions of time, whereas $\sigma$ and $R_\mathrm{c}$ are. Withal, $\rho(r, L_\mathrm{c}) \rightarrow \rho(r, L_\mathrm{c}, t)$.

#### 3.2.1 Expansion Velocity

The expansion velocity parameter in Equation (3.1) can be estimated from the average fragment velocity distribution. For short timescales, this parameter can be derived using the NASA Standard Breakup Model [5]. The average expansion velocity is determined by integrating over the fragment size distribution:

$$ \bar{v}_{\mathrm{expansion}} = \frac{\int_{L_\mathrm{c}} \bar{v}n dL_\mathrm{c}}{\int_{L_\mathrm{c}} n dL_\mathrm{c}} \tag{3.3} $$

Where $\bar{v}(L_\mathrm{c})$ represents the average velocity for fragments of characteristic length $L_\mathrm{c}$, which can be estimated using Equations (2.4) and (2.5) for explosions and collisions, respectively.

## 4 Orbital Mechanics

### 4.1 Two-Body Propagation of Fragment Orbits

Now that a time dependence for the PID cloud's density has been established, the focus will be shifted to modeling particles for individual fragment trajectories and the orbital mechanics of individual fragments.

#### 4.1.1 Individualistic Fragments

The cloud density model in Section 3 describes the statistical distribution of fragments, but individual fragments follow deterministic orbital paths. The transition from the statistical cloud model to individual fragment trajectories requires mapping the density distribution to specific initial state vectors. The probability is directly related to the mass density, $\mathcal P(\vec{r}|L_\mathrm{c}, t) \propto \rho(|\vec{r}|, L_\mathrm{c}, t)$. For a fragment with characteristic length $L_\mathrm{c}$ at time $t$ after impact, its position vector $\vec{r}$ corresponds to the probability density function:

$$ \mathcal P(\vec{r}|L_\mathrm{c}, t) = \frac{\rho(|\vec{r}|, L_\mathrm{c}, t)}{4\pi\int_0^{R_\mathrm{c}(t)} \rho(r, L_\mathrm{c}, t)r^2 dr} \tag{4.1} $$

This normalized probability density function relates the mass density distribution from Equation (2.1) to the spatial distribution of fragments across three-dimensional space. In this formulation, the radial distribution follows the Gaussian model while the angular components of $\vec{r}$ may vary according to impact dynamics.[^6]

The initial velocity vector for each fragment combines the parent object's velocity with the ejection velocity relative to the parent:

$$ \vec{v}_{\mathrm{fragment}} = \vec{v}_{\mathrm{parent}} + \vec{v}_{\mathrm{ejection}} \tag{4.2} $$

where $\vec{v}_{\mathrm{ejection}}$ is derived from the velocity distribution in the NASA Standard Breakup Model [1, 5]. For short timescales, this is the same as explosion/collision velocity given in Equations (2.4) and (2.5) respectively.

#### 4.1.2 Two-Body Orbital Dynamics

Once initial state vectors are established,[^7] each fragment's trajectory is governed by two-body orbital dynamics [4, 7]. The acceleration of a fragment in Earth's Newtonian gravitational field is given by,

$$ \vec{\ddot{r}}_{\mathrm{grav}} = -\frac{GM_\oplus}{r^2}\hat{r} \tag{4.3} $$

where $M_\oplus$ is Earth's mass and $|\vec{r}|$ is the distance of the particle from its center. For numerical propagation over relatively short timescales, the state vectors can be updated using a numerical integration scheme such as the fourth-order Runge-Kutta method.

#### 4.1.3 Size Effects on Orbital Evolution

The characteristic length influences orbital evolution through area-to-mass ratio effects [3, 8]. For larger fragments, gravitational forces dominate, while smaller fragments are increasingly affected by non-gravitational perturbations such as solar radiation pressure and atmospheric drag:

$$ \vec{a}_{\mathrm{total}} = \vec{\ddot{r}}_{\mathrm{grav}} + \vec{a}_{\mathrm{pert}} \tag{4.4} $$

The non-gravitational acceleration term, $\vec{a}_{\mathrm{pert}}$, is proportional to $A/M$ as these forces act on the cross-sectional area of the particle and the resulting acceleration depends inversely on mass due to Newton's second law:

$$ \vec{a}_{\mathrm{pert}} \simeq \frac{k_{\mathrm{pert}}}{L_\mathrm{c}} \tag{4.5} $$

where $k_{\mathrm{pert}}$ is a shape-dependent proportionality constant [3, 8].[^8]

### 4.2 Perturbation Forces

Building upon the relationship between fragment size and orbital evolution discussed in Section 4.1.3, this section focuses on non-gravitational perturbations, which are particularly relevant for modeling size-dependent effects in debris cloud evolution. While gravitational perturbations (Earth's non-spherical gravity field and third-body effects) do influence fragment orbits, they affect all fragments equally regardless of mass. For modeling the differential evolution of the Gaussian debris cloud, non-gravitational forces that scale with fragment size are of primary importance [3, 8].

#### 4.2.1 Non-Gravitational Perturbations Overview

The non-gravitational perturbation acceleration can be expressed as [3, 7]:

$$ \vec{a}_{\mathrm{pert}} = \vec{a}_{\mathrm{SRP}} + \vec{a}_d + \vec{a}_{\mathrm{EM}} \tag{4.6} $$

where $\vec{a}_{\mathrm{SRP}}$ represents solar radiation pressure, $\vec{a}_d$ represents atmospheric drag, and $\vec{a}_{\mathrm{EM}}$ accounts for electromagnetic forces. These forces are particularly significant because, unlike gravitational perturbations, their effect scales with the fragment's area-to-mass ratio [3, 8]. That being said, electromagnetic forces will be neglected because small debris fragments dominate large ones and non-gravitational perturbations like solar radiation pressure and atmospheric drag have a much stronger effect on smaller debris particles than electromagnetic forces [3, 6].

#### 4.2.2 Solar Radiation Pressure

Solar radiation pressure (SRP) becomes increasingly significant for fragments with high area-to-mass ratios, making it particularly important for smaller debris particles [3, 8]. The acceleration due to SRP can be modeled as,

$$ \vec{a}_{\mathrm{SRP}} = -P_{\mathrm{SRP}}\frac{A}{M}C_R\hat{r}_\odot \tag{4.7} $$

where $P_{\mathrm{SRP}}$ is the solar radiation pressure at the fragment's distance from the Sun (≈ 4.56 N/m$^2$ at Earth's orbit), $C_R$ is the dimensionless radiation pressure coefficient accounts for the optical properties of the fragments surface (between 1 and 2, typically 1.44 for debris fragments),[^9] and $\hat{r}_\odot$ is the unit vector from the fragment to the Sun [9, 10].

Following the relationship established in Equation 4.5, the solar radiation pressure's acceleration scales inversely with characteristic length:

$$ \vec{a}_{\mathrm{SRP}} \simeq -P_{\mathrm{SRP}}\frac{k}{L_\mathrm{c}}C_R\hat{r}_\odot \tag{4.8} $$

Please note that $k$ is numerically different than $k_{\mathrm{pert}}$. Anyhow, this inverse relationship explains why smaller fragments ($L_\mathrm{c} < 1$ cm) experience much stronger perturbations from SRP than larger fragments [8]. For these small fragments, SRP can cause significant orbital changes over relatively short timescales, affecting eccentricity and right ascension of the ascending node [3].

#### 4.2.3 Atmospheric Drag

For fragments in low Earth orbit (≲ 800 km), atmospheric drag represents the dominant non-gravitational perturbation [3, 7]. The acceleration due to drag follows:

$$ \vec{a}_d = -\frac{1}{2}\rho_{\mathrm{atm}}\frac{A}{M}C_d v_{\mathrm{rel}}^2\hat{v}_{\mathrm{rel}} \tag{4.9} $$

$$ \simeq -\frac{1}{2}\rho_{\mathrm{atm}}\frac{k}{L_\mathrm{c}}C_d v_{\mathrm{rel}}^2\hat{v}_{\mathrm{rel}} \tag{4.10} $$

where $\rho_{\mathrm{atm}}$ is the atmospheric density at the fragment's altitude, $C_d$ is the drag coefficient (usually between 2 and 2.3 for debris fragments [11]), and $v_{\mathrm{rel}}$ is the fragment's velocity relative to the atmosphere.

The atmospheric density varies significantly with altitude and solar activity, typically decreasing exponentially with increasing altitude [3, 7]. For debris cloud evolution modeling, empirical atmospheric models such as NRLMSISE-00 or Jacchia-Roberts are commonly used to estimate density at a given altitude, latitude, longitude, time, and solar conditions [7].

#### 4.2.4 Impact on Gaussian Distribution Parameters

A refined evolutionary model incorporating these effects would update $\sigma$ from Equation (2.6) to be time-dependent as it is now affected by SRP:

$$ \sigma(L_\mathrm{c}, t) = \sigma(L_\mathrm{c}, 0) + \gamma t L_\mathrm{c}^{-\alpha} \tag{4.11} $$

where $\gamma$ is an empirical parameter. The peak radial position, on the other hand, is unaffected by perturbations for short timescales. For large fragments, the SRP force is negligible as their mass-to-area ratio is very high [1]. For medium and small fragments, the initial deviation of fragments is dominated by the impulse from the fragmentation event itself [5]; over and above that, the ejection velocities from the fragmentation event far exceed the velocity changes induced by SRP within the first 120 seconds (typically ≪ 1 m·s$^{-1}$ for small fragments) [4].

Now that perturbations have been introduced causing $\sigma(L_\mathrm{c}) \rightarrow \sigma(L_\mathrm{c}, t)$ and $\mu(L_\mathrm{c}) \not\rightarrow \mu(L_\mathrm{c}, t)$, the fragment normal distribution with all its functionality is:

$$ \rho(r, L_\mathrm{c}, t) = \rho_0(L_\mathrm{c}, t) \exp\left[-\frac{1}{2}\left(\frac{r - \mu(L_\mathrm{c})R_\mathrm{c}(t)}{\sigma(L_\mathrm{c}, t)R_\mathrm{c}(t)}\right)^2\right] \tag{4.12} $$

with $R_\mathrm{c}(t)$, $\rho_0(L_\mathrm{c}, t)$, and $\sigma(L_\mathrm{c}, t)$ given in Equations (3.1), (3.2), and (4.11) respectively.

## Acknowledgments

Thank you to my colleague and mentor, Robert A. Amborski, for his invaluable guidance and support throughout this research. His generous provision of key references significantly strengthened the theoretical foundation of this work. Robert's previous research into finding the probability of 
impact given an interceptor flying through a cloud of debris is what directly and ultimately created the postulate that fragment dispersion post-impact is 
symmetric about the peak of the distribution. I would also like to acknowledge the technical discussions and constructive feedback that helped refine the mathematical framework presented in this paper.

## References

[1] Johnson, N. L., Krisko, P. H., Liou, J.-C., & Anz-Meador, P. D. (2001). NASA’s New Breakup Model of EVOLVE 4.0. _Advances in Space Research_, **28**(9), 1377, 1384.

[2] Liou, J.-C., Johnson, N. L., & Hill, N. (2001). Risks in Space from Orbiting Debris. _\mathrm{s}cience_, **311**(5759), 340–341.

[3] Wright, D. (2007). Space Debris. _Physics Today_, **60**(7), 36–42.

[4] Chobotov, V. A., & Spencer, D. B. (1990). A Review of Orbital Debris Modeling Techniques, _Journal of Spacecraft and Rockets_, **27**(3), 261–269.

[5] Nagl, D. E., Reynolds, R. C., & McKay, G. A. (1992). Fragmentation Models: Techniques and Validation. _Advances in Space Research_, **12**(8), 107–116.

[6] Braun, V., Oltrogge, D., & Alfano, S. (2017). Statistical Analysis of Orbital Breakup Events. _Acta Astronautica_, **139**, 56–68.

[7] Vallado, D. A., & Oltrogge, D. L. (2017). Three-Dimensional Volumetric Risk Assessment for Fragmentation Event Debris Fields. _Journal of Geophysical Research: Space Physics_, **122**(11), 11,213–11,228.

[8] Cimmino, N., Isoletta, G., Piergentili, F., Santoni, F., Cardona, T., Graziani, F., Fuentes, J. L., & Cecchini, A. (2021). Tuning of NASA Standard Breakup Model for Fragmentation Events Modelling. _Aerospace_, **8**(7), 185.

[9] Padole, P. (2024). How Solar Radiation Pressure Affects Satellites: Calculating Instantaneous Force and Acceleration. _Imperial College London Blogs_. https://blogs.imperial.ac.uk/parikshit/2024/09/23/solar-radiation-pressure-satellites/

[10] List, M., Bremer, S., Rievers, B., & Lämmerzahl, C. (2015). Modelling of Solar Radiation Pressure Effects: Parameter Analysis for the MICROSCOPE Mission. _International Journal of Aerospace Engineering_, **2015**, 928206.

[11] Anz-Meador, P. D., & Liou, J.-C. (2001). An assessment of the NASA explosion fragmentation model to 1 mm characteristic sizes. _Advances in Space Research_, **28**(9), 1383-1388.

## Appendix A: Reasoning Behind $\mathcal P \propto \rho$

Equation (4.1) represents the relationship between the Gaussian mass density function $\rho(r, L_\mathrm{c}, t)$ and a probability density function for fragment positions in 3D space. The mass density function describes how fragment mass is distributed radially within the cloud; to create a probability density function in space, we need to ensure it integrates to 1 over all space, $\iiint_\mathcal V \mathcal P d\mathcal V = 1$. In spherical coordinates this becomes:

$$ 4\pi \int_0^{R_\mathrm{c}(t)} \mathcal P(r|L_\mathrm{c}, t)r^2 dr = 1 \tag{A.1} $$

For the probability to be proportional to the mass density, dividing $\rho$ by the total mass normalizes this distribution:

$$ \mathcal P(\vec{r}|L_\mathrm{c}, t) = \frac{\rho(|\vec{r}|, L_\mathrm{c}, t)}{M_{\mathrm{total}}(t)} \tag{A.2} $$

where,

$$ M_{\mathrm{total}}(t) = 4\pi \int_0^{R_\mathrm{c}(t)} \rho(r, L_\mathrm{c}, t)r^2 dr \tag{A.3} $$

This normalized probability distribution provides the mathematical foundation for accurately positioning fragments according to their physical distribution while maintaining statistical validity through proper conversion of mass density into sampling weights. A sampling weight is a numerical value assigned to each potential position in space that determines its probability of being selected during random sampling. In this context, positions with higher mass density receive higher weights, making them more likely to be chosen when generating fragment positions.

## Appendix B: Size Classification

Size classifications in this model are as the NASA EVOLVE 4.0 model [1]:

-   Small ($L_\mathrm{c} \lesssim 8~\text{cm}$); experience rapid dispersion due to dominant SRP and atmospheric drag effects [3, 8]
-   Medium ($8~\text{cm}\lesssim L_\mathrm{c} \lesssim 11~\text{cm}$); experience moderate deviation primarily due to SRP [3]
-   Large ($L_\mathrm{c} \gtrsim 11$ cm); remain close to their initial ballistic trajectories [1]

This size-dependent stratification alters the Gaussian density distribution over time, with different size classes evolving at different rates according to their respective area-to-mass ratios.

## Appendix C: Relationship Between $L_\mathrm{c}$ and $A/M$

The relationship between characteristic length and area-to-mass ratio in the NASA Breakup Model is not a simple direct formula, but rather a complex set of statistical distributions. First, $L_\mathrm{c}$ is related to cross-sectional area by:

$$ 
	A = \begin{cases} 
		0.540424 L_\mathrm{c}^2 & \text{if } L_\mathrm{c} < 0.00167 \text{ m} \\
		0.556945 L_\mathrm{c}^{2.0047077} & \text{otherwise} \end{cases} 
\tag{C.1}
$$

Then, $A/M$ follows statistical distributions based on $L_\mathrm{c}$ [1].

[^1]: Packing density, packing fraction, and packing coefficient are all correct terminology.

[^2]: If the peak density position is at two-thirds of the cloud radius, then the mean value of the distribution is $\frac{2}{3}R_\mathrm{c}$.

[^3]: Equation (2.2) is obtained by breaking down the cloud sphere into infinitesimally thin shells with radii between 0 and $R_\mathrm{c}$.

[^4]: Refer to Appendix B for fragment size classifications.

[^5]: Note that fragments with the same characteristic length do not all have the same area or mass. In earlier models of the NASA Standard Breakup Model, debris was assumed to be spherical, so a given $L_\mathrm{c}$ would correspond to a specific area and mass. So, $\log \frac{A}{M}$ is related to $L_\mathrm{c}$ through a statistical distribution rather than a direct function; hence, it is not the same for identical $L_\mathrm{c}$.

[^6]: Refer to Appendix A for more details.

[^7]: The position vector ($\vec{r}$) and the velocity vector ($\vec{v}_{\mathrm{fragment}}$) define the initial orbital state of each particle at the moment of separation from the parent object.

[^8]: The constant $k_{\mathrm{pert}}$ captures shape and material properties that affect this relationship. For example, for a spherical fragment $k_{\mathrm{pert}} = \frac{3}{\rho}$, for a cube fragment $k_{\mathrm{pert}} = \frac{6}{\rho}$, and for irregular shaped fragments, $k_{\mathrm{pert}}$ would be determined empirically.

[^9]: For perfect absorption, $C_R = 1$ (momentum is transferred only once); for perfect reflection, $C_R = 2$ (momentum is transferred twice, once upon absorption and once upon reflection).
