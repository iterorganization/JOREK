---
title: "Run with Diamagnetic Drift"
nav_order: 1
parent: "Activation of Physical Terms"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# Running with Diamagnetic Drift Terms

- Diamagnetic drift terms are currently implemented in models 303, 333, 400, 500, 555, and **600**.
- To switch the terms on, set the input parameter `tauIC` in the namelist input file to the correct value:

$$
\tau_\text{IC} = \frac{m_\text{ion}}{2e\,F_0\sqrt{\mu_0\rho_0}}
$$

- This expression assumes $T_\text{ion}=T_\text{el}$, which accounts for the factor of two in the denominator.
- By default, the diamagnetic terms are switched off: `tauIC` is preset to zero.
- According to the definition above, the `tauIC` input parameter is always defined as for a single-temperature model.

## Examples

For an ASDEX Upgrade deuterium plasma with $n_0=10^{20}\,\mathrm{m}^{-3}$ and $F_0=4.2$, the correct value of `tauIC` is $3.8\times10^{-3}$:

$$
\begin{aligned}
\mu_0 &= 4\pi\times10^{-7}\,\mathrm{N}\,\mathrm{A}^{-2}, \\
e &= 1.6022\times10^{-19}\,\mathrm{C}, \\
m_D &= 3.344\times10^{-27}\,\mathrm{kg}, \\
\rho_0 &= 3.344\times10^{-7}\,\mathrm{kg}\,\mathrm{m}^{-3}, \\
\sqrt{\mu_0\rho_0} &= 6.5\times10^{-7}.
\end{aligned}
$$

For a deuterium plasma, the expression can be written as:

$$
\tau_\text{IC} = \frac{1.61\times10^8}{F_0\sqrt{n_0}}.
$$

## Some estimates

- The **diamagnetic velocity** is $\mathbf{v}_\text{dia}=-(\nabla p\times\mathbf{B})/(qnB^2)$, with magnitude $v_\text{dia}=|\nabla p|/(qnB)$, where $q$ is the particle's electric charge and $n$ is the particle density.
- The **diamagnetic frequency** is given by $\tau_\text{dia}=v_D/C$, where $C\approx2\pi r$ is the poloidal circumference of the flux surface under consideration. At the edge of an ASDEX Upgrade plasma, $C\approx4\,\mathrm{m}$.
