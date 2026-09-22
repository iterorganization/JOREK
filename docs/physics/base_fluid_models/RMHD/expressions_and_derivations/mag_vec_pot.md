---
title: "Magnetic Vector Potential"
nav_order: 6
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Magnetic Vector Potential

The ansatz for the magnetic vector potential is ***(this ansatz is not relevant for the derivation of the JOREK equations; details about the Gauge missing)***:

$$\begin{equation*}
\mathbf{A}=\nabla\alpha\times\nabla\phi+\psi\nabla\phi
\end{equation*}$$

where $\partial_3\alpha\equiv 0$ is assumed. The curl of $\mathbf{A}$ is given by

$$\begin{equation*}\begin{split}
\nabla\times\mathbf{A}
 &= \nabla\times(\nabla\alpha\times\nabla\phi)+\nabla\times(\psi\nabla\phi) \\
 &= -\Delta^*\alpha\nabla\phi+\frac{1}{R^2}\nabla_\text{pol}\underbrace{\partial_3\alpha}_{\equiv 0} + \nabla\psi\times\nabla\phi + \psi\underbrace{\nabla\times\nabla\phi}_{\equiv 0} \\
 &= -\Delta^*\alpha\nabla\phi+\nabla\psi\times\nabla\phi
\end{split}\end{equation*}$$

As this has to be in agreement with the magnetic field expression $\mathbf{B}=F\nabla\phi+\nabla\psi\times\nabla\phi$, the following condition can be derived:

$$\begin{equation*}
  F\equiv-\Delta^*\alpha
\end{equation*}$$
