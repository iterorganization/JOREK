---
title: "The Continuity Equation Expression"
nav_order: 12
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the Continuity Equation

The continuity equation in JOREK is given by:

$$\begin{equation}
  \partial_t\rho=-\nabla\cdot(\rho\mathbf{v})+\nabla\cdot(D\nabla\rho) + S_\rho
\end{equation}$$

The expression for the velocity $\mathbf{v}$ in the reduced MHD model is:

$$\begin{equation*}
\mathbf{v} = -R^2\nabla u \times \nabla\phi + v_{||}\mathbf{B}
\end{equation*}$$

Substituting this into the first term of the right hand side of the continuity equation, we get:

$$\begin{equation*}\begin{split}
-\nabla\cdot(\rho\mathbf{v}) &= -\nabla\cdot\left[\rho\left(-R^2\nabla u \times \nabla\phi + v_{||}\mathbf{B}\right)\right] \\
&= \nabla\cdot\left[\rho R^2\nabla u \times \nabla\phi\right] - \nabla\cdot\left(\rho v_{||}\mathbf{B}\right)
\end{split}\end{equation*}$$

Using the vector identity $\nabla\cdot(\mathbf{A}\times\mathbf{C}) = \mathbf{C}\cdot(\nabla\times\mathbf{A}) - \mathbf{A}\cdot(\nabla\times\mathbf{C})$ and the fact that $\nabla\times\nabla\phi = 0$, the first term becomes:

$$\begin{equation*}\begin{split}
\nabla\cdot\left(\rho R^2\nabla u \times \nabla\phi\right) &= \nabla\phi\cdot\left[\nabla\times\left(\rho R^2\nabla u\right)\right] \\
&= \nabla\phi\cdot\left[\nabla(\rho R^2)\times\nabla u + \rho R^2\underbrace{\nabla\times\nabla u}_{=0}\right] \\
&= \nabla\phi\cdot\left[\nabla(\rho R^2)\times\nabla u\right] \\
&= \frac{1}{R}[\rho R^2, u]
\end{split}\end{equation*}$$

For the second term, using $\nabla\cdot\mathbf{B} = 0$, we have:

$$\begin{equation*}\begin{split}
-\nabla\cdot\left(\rho v_{||}\mathbf{B}\right) &= -\mathbf{B}\cdot\nabla(\rho v_{||}) - \rho v_{||}\underbrace{\nabla\cdot\mathbf{B}}_{=0} \\
&= -\mathbf{B}\cdot\nabla(\rho v_{||}) \\
&= -\frac{F}{R^2}\partial_3(\rho v_{||}) - \frac{1}{R}[\rho v_{||}, \psi]
\end{split}\end{equation*}$$

Combining these, the continuity equation in reduced MHD variables is:

$$\begin{equation}
\partial_t\rho = \frac{1}{R}[\rho R^2, u] - \frac{1}{R}[\rho v_{||}, \psi] - \frac{F}{R^2}\partial_3(\rho v_{||}) + \nabla\cdot(D\nabla\rho) + S_\rho
\end{equation}$$
