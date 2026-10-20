---
title: "The Poloidal Flux Equation Expression"
nav_order: 13
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the Poloidal Flux Equation

The starting point for the poloidal flux equation is the induction equation:

$$\begin{equation}
  \partial_t\mathbf{B} = -\nabla\times\mathbf{E}
\end{equation}$$

Expressing the magnetic field as $\mathbf{B}=\nabla\times\mathbf{A}$ and removing the curl leads to:

$$\begin{equation}\label{dAdt}
  \partial_t\mathbf{A} = -\mathbf{E} - \nabla\Phi + \alpha_L\nabla\phi
\end{equation}$$

where $\Phi$ is the electric potential and $\alpha_L$ represents the loop voltage. Multiplying this equation by $R^2\nabla\phi$, we get:

$$\begin{equation*}
  R^2\nabla\phi\cdot\partial_t\mathbf{A} = -R^2\nabla\phi\cdot\mathbf{E} - R^2\nabla\phi\cdot\nabla\Phi + R^2\nabla\phi\cdot\alpha_L\nabla\phi
\end{equation*}$$

Using the definition of the poloidal flux $\psi = R^2\mathbf{A}\cdot\nabla\phi = RA_\phi$ and the relation for the electric potential $F_0 u \equiv \Phi$, the terms are:

$$\begin{align*}
  R^2\nabla\phi\cdot\partial_t\mathbf{A} &= \partial_t\psi \\
  -R^2\nabla\phi\cdot\nabla\Phi &= -R^2\frac{1}{R^2}\partial_3\Phi = -\partial_3\Phi = -F_0\partial_\phi u \\
  R^2\nabla\phi\cdot\alpha_L\nabla\phi &= \alpha_L
\end{align*}$$

The electric field term $-R^2\nabla\phi\cdot\mathbf{E}$ is evaluated using Ohm's law $\mathbf{E} = -\mathbf{v}\times\mathbf{B} + \eta\mathbf{j}$:

$$\begin{equation*}\begin{split}
  -R^2\nabla\phi\cdot\mathbf{E} &= R^2\nabla\phi\cdot(\mathbf{v}\times\mathbf{B}) - R^2\nabla\phi\cdot\eta\mathbf{j} \\
  &= R^2\frac{1}{R}[\psi, u] - R^2\eta\frac{1}{R^2}(j - j_0) \\
  &= R[\psi, u] - \eta(j - j_0)
\end{split}\end{equation*}$$

where $j_0$ is the current source related to the loop voltage $\alpha_L$ via $\eta j_0 = \alpha_L$. 

Combining all the terms, we obtain the poloidal flux equation used in JOREK:

$$\begin{equation}
  \partial_t\psi = R[\psi, u] + \eta(j - j_0) - F_0\partial_\phi u
\end{equation}$$
