---
title: "The Diamagnetic Vorticity Expression"
nav_order: 14
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the Diamagnetic Vorticity

In JOREK models that include diamagnetic effects, the total vorticity $W$ is defined as the toroidal component of the curl of the ion velocity $\mathbf{v}_i$:

$$\begin{equation*}
  W \equiv (\nabla \times \mathbf{v}_i) \cdot \nabla \phi
\end{equation*}$$

The ion velocity is given by:

$$\begin{equation*}
  \mathbf{v}_i = \mathbf{v}_{E \times B} + \mathbf{v}_{dia,i} + v_{||} \mathbf{B}
\end{equation*}$$

Substituting the expressions for the $E \times B$ and diamagnetic velocities:

$$\begin{equation*}
  \mathbf{v}_i = -R^2 \nabla u \times \nabla \phi - \frac{\tau_{IC} R^2}{\rho} \nabla p_i \times \nabla \phi + v_{||} \mathbf{B}
\end{equation*}$$

The toroidal vorticity $W$ can be split into two main parts:

$$\begin{equation*}
  W = \omega + \omega^*
\end{equation*}$$

where $\omega$ is the standard vorticity related to the $E \times B$ flow:

$$\begin{equation*}
  \omega = \nabla \cdot \nabla_{pol} u
\end{equation*}$$

and $\omega^*$ is the part related to the diamagnetic flow and the parallel velocity term:

$$\begin{equation*}
  \omega^* = \nabla \cdot \left( \frac{\tau_{IC}}{\rho} \nabla_{pol} p_i - v_{||} \frac{1}{R^2} \nabla_{pol} \psi \right)
\end{equation*}$$

The full derivation of these terms, including the cancellation of terms with the gyro-viscous tensor, can be found in the [Vorticity Derivation](omega_deriv.md) and [RMHD Model](../rmhd_model.md) pages.
