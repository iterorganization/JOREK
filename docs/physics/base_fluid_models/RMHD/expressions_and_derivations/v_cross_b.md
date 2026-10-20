---
title: "The v x B Term Expression"
nav_order: 11
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the v x B Term

$$\begin{equation*}\begin{split}
  \mathbf{v}\times\mathbf{B}
&= (-R^2\nabla u\times\mathbf{a}^3 + v_{||}\mathbf{B})\times\mathbf{B} \\
&= (-R^2\nabla u\times\mathbf{a}^3)\times\mathbf{B} \\
&= (-R^2\nabla u\times\mathbf{a}^3)\times(\nabla\psi\times\mathbf{a}^3+F_0\mathbf{a}^3) \\
&= -R^2(\nabla u\times\mathbf{a}^3)\times(\nabla\psi\times\mathbf{a}^3)-R^2 F_0(\nabla u\times\mathbf{a}^3)\times\mathbf{a}^3 \\
&= -R^2[\nabla\psi\underbrace{(\nabla u\times\mathbf{a}^3)\cdot\mathbf{a}^3}_{=0}-\mathbf{a}^3((\nabla u\times\mathbf{a}^3)\cdot\nabla\psi)]+R^2 F_0[\nabla u\underbrace{(\mathbf{a}^3\cdot\mathbf{a}^3)}_{=1/R^2}-\mathbf{a}^3(\nabla u\cdot\mathbf{a}^3)] \\
&= R^2\mathbf{a}^3[\nabla\psi\cdot(\nabla u\times\mathbf{a}^3)]+F_0\nabla u-R^2 F_0\mathbf{a}^3(\partial_3 u\underbrace{\mathbf{a}^3\cdot\mathbf{a}^3}_{=1/R^2}) \\
&= R\mathbf{a}^3(\partial_2 u\partial_1\psi - \partial_1 u\partial_2\psi)+F_0(\nabla u-\mathbf{a}^3\partial_3 u) \\
&= R\mathbf{a}^3[\psi,u] + F_0\nabla_\text{pol} u
\end{split}\end{equation*}$$

Taking diamagnetic effects into account leads to:

$$\begin{equation*}\begin{split}
  \mathbf{v}\times\mathbf{B}
&= (-R^2\nabla u\times\mathbf{a}^3 + v_{||}\mathbf{B} {\color{grey}{- \frac{2 \tau_{IC}R^2}{\rho}\nabla p_i \times \mathbf{a}^3}})\times\mathbf{B} \\
&= R\mathbf{a}^3[\psi,u] + F_0\nabla_\text{pol} u {\color{grey}{+ \frac{2 \tau_{IC}R}{\rho}\mathbf{a}^3[\psi,p_i] + \frac{2 \tau_{IC}}{\rho}F_0\nabla_\text{pol}p_i}}
\end{split}\end{equation*}$$

Also:

$$\begin{equation*}\begin{split}
  \mathbf{a^3} \cdot \left( \mathbf{v}\times\mathbf{B} \right) =
  \mathbf{B} \cdot \left( \mathbf{a^3}\times\mathbf{v} \right) = - \mathbf{B} \cdot \nabla_{pol} u = - \frac{1}{R} \left[u,\psi \right] {\color{grey}{- \frac{2\tau_{IC}}{\rho R}F_0[p_i,\psi] }}
\end{split}\end{equation*}$$

Information from the [coordinates](../../../coordinates.md) as well as the [vector-identities](../../../vector-identities.md) pages was used.
