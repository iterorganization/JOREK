---
title: "The j x B Term Expression"
nav_order: 10
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the j x B Term

$$\begin{equation*}\begin{split}
\mathbf{j}\times\mathbf{B} 
=& \left(\nabla F\times\nabla\phi+\frac{1}{R^2}\nabla_\text{pol}\partial_3\psi-\Delta^*\psi\nabla\phi\right)\times(F\nabla\phi+\nabla\psi\times\nabla\phi) \\
=& -F\nabla\phi\times(\nabla F\times\nabla\phi) + (\nabla F\times\nabla\phi)\times(\nabla\psi\times\nabla\phi) + \frac{F}{R^2}\nabla_\text{pol}\partial_3\psi\times\nabla\phi \\
 & +\frac{1}{R^2}\nabla_\text{pol}\partial_3\psi\times(\nabla_\text{pol}\psi\times\nabla\phi)-F\Delta^*\psi\underbrace{\nabla\phi\times\nabla\phi}_{\equiv 0}-\Delta^*\psi\nabla\phi\times(\nabla\psi\times\nabla\phi) \\
=&-\frac{F}{R^2}\nabla_\text{pol}F+\frac{1}{R}\nabla\phi[F,\psi]-\frac{\Delta^*\psi}{R^2}\nabla_\text{pol}\psi+\frac{1}{R^2}\nabla_\text{pol}\partial_3\psi\times\mathbf{B}
\end{split}\end{equation*}$$
