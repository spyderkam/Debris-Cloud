# Size Distribution

Kamyar Modjtahedzadeh  
Boeing Intelligence & Analytics  
May 22, 2025

> What is the quantitative distribution of fragment size (small, medium, large) immediately post impact? What about distribution of specific $L_\mathrm{c}$? How _many_ total fragments for each $L_\mathrm{c}$ to _almost_ satisfy mass conservation (doesn’t and _shouldn’t_ have to be exact!)?

## Framework from NASA EVOLVE 4.0

$n(L_\mathrm c)$ is the differential number density of fragments at _exactly_ size $L_\mathrm c$ and $N(L_\mathrm c)$ represents the cumulative number of fragments *larger* than $L_{\mathrm{c}}$.

$$n(L_\mathrm c) = -\frac{d}{dL_{\mathrm{c}}} N(L_{\mathrm{c}})$$
  
### Fragment Counts by Size Category

To determine the number of fragments within specific size ranges, use the cumulative distribution function, 
$$N_{\mathrm{collision}}(L_{\mathrm{c}}) \,=\, 0.1 L_{\mathrm{c}}^{-1.71} M_{\mathrm{parent}}^{0.75}$$

where $N_{\mathrm{collision}}$ will just be written as $N$ for short. The number of fragments between two characteristic lengths $L_1$ and $L_2$ (where $L_1 < L_2$) is calculated as:

$$\text{fragments between } L_1 \text{ and } L_2 \,=\, N(L_1) - N(L_2)$$

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
