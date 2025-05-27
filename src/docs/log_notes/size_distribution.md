# Numerical Size Distribution

Kamyar Modjtahedzadeh  
Boeing Intelligence & Analytics  
May 22, 2025 -- May 27, 2025

> What is the quantitative distribution of fragment size (small, medium, large) immediately post impact? What about distribution of specific $L_\mathrm{c}$? How _many_ total fragments for each $L_\mathrm{c}$ to _almost_ satisfy mass conservation (doesn’t and _shouldn’t_ have to be exact!)?

## Framework from NASA EVOLVE 4.0

$n(L_\mathrm c)$ is the differential number density of fragments at _exactly_ size $L_\mathrm c$ and $N(L_\mathrm c)$ represents the cumulative number of fragments *larger* than $L_{\mathrm{c}}$.

$$n(L_\mathrm c) = -\frac{d}{dL_{\mathrm{c}}} N(L_{\mathrm{c}})$$
  
### Fragment Counts by Size Category

To determine the number of fragments within specific size ranges, use the cumulative distribution function, 
$$N_{\mathrm{collision}}(L_{\mathrm{c}}) \,=\, 0.1 L_{\mathrm{c}}^{-1.71} M_{\mathrm{parent}}^{0.75}$$

where $N_{\mathrm{collision}}$ will just be written as $N$ for short. The number of fragments between two characteristic lengths $L_1$ and $L_2$ (where $L_1 < L_2$) is calculated as:

$$\text{number of fragments between } L_1 \text{ and } L_2 \,=\, N(L_1) - N(L_2)$$

The practical bounds for all size categories are set to $L_{\mathrm{min}} = 0.001 \, \mathrm{m}$ and $L_{\mathrm{max}} = 1.0 \, \mathrm{m}$.

  - Small fragments ($0.001 \, \mathrm{m} \leq L_{\mathrm{c}} \leq 0.08 \, \mathrm{m}$):
	$$
	\begin{align*}
		N_{\mathrm{small}} \,=\, N(0.001) - N(0.08) \,=\, 13482.1 M_{\mathrm{parent}}^{0.75}
	\end{align*}
	$$
  - Medium fragments ($0.08 \, \mathrm{m} \leq L_{\mathrm{c}} \leq 0.11 \, \mathrm{m}$)
	$$
	\begin{align*}
		N_{\mathrm{medium}} \,=\, N(0.08) - N(0.11) \,=\, 3.15  M_{\mathrm{parent}}^{0.75}
	\end{align*}
	$$
  - Large fragments ($0.11 \, \mathrm{m} \leq L_{\mathrm{c}} \leq 1.0 \, \mathrm{m}$)
	$$
	\begin{align*}
		N_{\mathrm{large}} \,=\, N(0.11) - N(1) \,=\, 4.26 M_{\mathrm{parent}}^{0.75}
	\end{align*}
	$$

### Percentage Distribution Calculation

The total number of fragments is:

$$N_{\mathrm{total}} \,=\, N_{\mathrm{small}} + N_{\mathrm{medium}} + N_{\mathrm{large}} \,=\, 13489.5 M_{\mathrm{parent}}^{0.75}$$

The percentage distributions are:

  - Small fragments: $99.95\%$
  - Medium fragments: $0.023\%$
  - Large fragments: $0.032\%$

## Internal Distribution

It has been established how the fragments are distributed based on their size stratifications; however, the three categories also have internal distributions.

### Distribution Within Medium Fragments

The medium fragments themselves follow the same power-law distribution pattern. Within the medium size category, the fragments are *not* uniformly distributed as smaller fragments within this range are much more numerous than larger ones.

#### Mathematical Framework for Sub-Distributions

For any subrange within the medium category, from $L_{\mathrm{a}}$ to $L_{\mathrm{b}}$, where $0.08 \leq L_{\mathrm{a}} < L_{\mathrm{b}} \leq 0.11$, the number of fragments is:

$$N_{\mathrm{subrange}} = N(L_{\mathrm{a}}) - N(L_{\mathrm{b}}) = 0.1 M_{\mathrm{parent}}^{0.75} \left(L_{\mathrm{a}}^{-1.71} - L_{\mathrm{b}}^{-1.71}\right)$$

#### Quantitative Distribution Within Medium Range

Here's how the $3.15 M_{\mathrm{parent}}^{0.75}$ medium fragments are distributed:

| Size Range [$\mathrm{cm}$] | Fragment Count | $\%$ of Medium | $\%$ of Total |
|----------------------------|----------------|----------------|---------------|
| $8.0$--$8.5$    | $0.74M_{\mathrm{parent}}^{0.75}$ | $23.5\%$       | $0.0055\%$    |
| $8.5$--$9.0$    | $0.63M_{\mathrm{parent}}^{0.75}$ | $20.0\%$       | $0.0047\%$    |
| $9.0$--$9.5$    | $0.54M_{\mathrm{parent}}^{0.75}$ | $17.2\%$       | $0.0040\%$    |
| $9.5$--$10.0$   | $0.47M_{\mathrm{parent}}^{0.75}$ | $14.9\%$       | $0.0035\%$    |
| $10.0$--$10.5$  | $0.41M_{\mathrm{parent}}^{0.75}$ | $13.0\%$       | $0.0030\%$    |
| $10.5$--$11.0$  | $0.36M_{\mathrm{parent}}^{0.75}$ | $11.4\%$       | $0.0027\%$    |

#### Key Observations

  - **Size Stratification Within Medium Range**: The $8.0$ to $8.5$ centimeter subrange contains nearly **twice as many fragments** as the $10.5$ to $11$ centimeter subrange, despite both being $0.5\, \mathrm{cm}$ wide intervals.

  - **Power-Law Behavior**: The distribution follows the same $L_{\mathrm{c}}^{-1.71}$ power-law within the medium range, meaning smaller fragments always dominate numerically.

  - **Relative Fragment Density**: At any specific characteristic length within the medium range, the fragment density can be approximated using the differential distribution: $$n(L_{\mathrm{c}}) = 0.171 M_{\mathrm{parent}}^{0.75} L_{\mathrm{c}}^{-2.71}$$ For example, at $L_{\mathrm{c}} = 8.0\,\mathrm{cm}$, there are approximately **twice as many fragments per unit length** compared to $L_{\mathrm{c}} = 10.0\,\mathrm{cm}$.

This same stratification pattern applies to small and large fragment categories as well - within each category, smaller fragments are always more numerous than larger ones due to the underlying power-law physics of fragmentation processes.

### Small and Large Fragments

### Small Fragments

| Size Range [$\mathrm{cm}$] | Fragment Count | $\%$ of Small | $\%$ of Total |
|----------------------------|----------------|---------------|---------------|
| $0.1$--$1.0$    | $13227M_{\mathrm{parent}}^{0.75}$ | $98.1\%$      | $98.1\%$      |
| $1.0$--$2.0$    | $183M_{\mathrm{parent}}^{0.75}$   | $1.4\%$       | $1.4\%$       |
| $2.0$--$3.0$    | $40M_{\mathrm{parent}}^{0.75}$    | $0.3\%$       | $0.3\%$       |
| $3.0$--$4.0$    | $16M_{\mathrm{parent}}^{0.75}$    | $0.1\%$       | $0.1\%$       |
| $4.0$--$5.0$    | $8M_{\mathrm{parent}}^{0.75}$     | $0.1\%$       | $0.1\%$       |
| $5.0$--$6.0$    | $5M_{\mathrm{parent}}^{0.75}$     | $0.0\%$       | $0.0\%$       |
| $6.0$--$7.0$    | $3M_{\mathrm{parent}}^{0.75}$     | $0.0\%$       | $0.0\%$       |
| $7.0$--$8.0$    | $2M_{\mathrm{parent}}^{0.75}$     | $0.0\%$       | $0.0\%$       |

The power-law distribution creates **extreme stratification** within small fragments. The smallest subrange ($0.1$ to $1.0$ centimeter) contains $\mathbf{98.1\%}$ of all small fragments, which represents **$\boldsymbol{98.1\%}$ of the entire debris population**. This means that fragments smaller than $1 \,\mathrm{cm}$ constitute nearly the entire debris cloud by count. The remaining small fragment subranges show rapid decline: the $1.0$ to $2.0$ centimeter range contains only $1.4\%$ of small fragments, and all subranges above $2.0 \,\mathrm{cm}$ each contain less than $0.5\%$ of the small fragment population. This creates a heavily skewed distribution where the tiniest debris pieces overwhelmingly dominate.

### Large Fragments

| Size Range [$\mathrm{cm}$] | Fragment Count | $\%$ of Large | $\%$ of Total |
|----------------------------|----------------|---------------|---------------|
| $11$--$20$      | $2.79M_{\mathrm{parent}}^{0.75}$ | $65.5\%$      | $0.0207\%$    |
| $20$--$30$      | $0.78M_{\mathrm{parent}}^{0.75}$ | $18.4\%$      | $0.0058\%$    |
| $30$--$40$      | $0.30M_{\mathrm{parent}}^{0.75}$ | $7.2\%$       | $0.0023\%$    |
| $40$--$50$      | $0.15M_{\mathrm{parent}}^{0.75}$ | $3.6\%$       | $0.0011\%$    |
| $50$--$60$      | $0.09M_{\mathrm{parent}}^{0.75}$ | $2.1\%$       | $0.0006\%$    |
| $60$--$70$      | $0.06M_{\mathrm{parent}}^{0.75}$ | $1.3\%$       | $0.0004\%$    |
| $70$--$80$      | $0.04M_{\mathrm{parent}}^{0.75}$ | $0.9\%$       | $0.0003\%$    |
| $80$--$100$     | $0.05M_{\mathrm{parent}}^{0.75}$ | $1.1\%$       | $0.0003\%$    |

The distribution across large fragment subranges is more spread out: the $20.0$ to $30.0$ centimeter range still contains a meaningful $\mathbf{18.4\%}$ of large fragments, and even the $30.0$ to $40.0$ centimeter range holds $\mathbf{7.2\%}$. This indicates that while smaller sizes still dominate within the large category, the distribution is more balanced across different size subranges compared to the extreme concentration observed in small fragments. The large fragments represent the most **trackable and observable** debris pieces, despite comprising only **$\boldsymbol{0.032\%}$ of the total fragment count**.

## Discrete Fragment Count

### Size Step Configuration

The computational formulation employs different increment sizes optimized for each fragment category:

**Small fragments** ($0.001$ to $0.08$ meters) use increments of $\Delta L_{\mathrm{c}} = 0.0001$ meters, creating $790$ discrete calculation intervals. This fine granularity captures the steep power-law distribution where fragment counts change rapidly over small size ranges.

**Medium fragments** ($0.08$ to $0.11$ meters) use increments of $\Delta L_{\mathrm{c}} = 0.005$ meters, creating $60$ discrete calculation intervals. The larger increment size reflects the more moderate rate of change in this transition region.

**Large fragments** ($0.11$ to $1.0$ meters) use increments of $\Delta L_{\mathrm{c}} = 0.001$ meters, creating $890$ discrete calculation intervals. This increment size provides adequate resolution across the extended size range while maintaining computational efficiency.

### Interval Sampling <!--Size Increment Processing-->

To estimate the number of fragments for any particular characteristic length, the following formula is used:

$$\Upsilon(L_{\mathrm{c}}) \,=\, N(L_{\mathrm{c}}) - N(L_{\mathrm{c}} +\epsilon)$$

where the size step $\epsilon$ is set to `0.00001` meters. With this formulation, $\Upsilon(L_{\mathrm{c}}) \approx {\epsilon} \cdot n(L_{\mathrm{c}})$. Loop through the discrete $L_{\mathrm{c}}$ values with the assigned $\Delta L_{\mathrm{c}}$ increments and calculate each `nFrag` with the above equation for $\Upsilon(L_{\mathrm{c}})$.

