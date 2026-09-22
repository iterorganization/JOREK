---
title: "Derivation of the Parallel Momentum Equation"
nav_order: 15
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Appendix: Derivation of the Parallel Momentum Equation

The parallel projection of the momentum equation is obtained by applying the operator $\nabla\phi \cdot (\dots)$ onto the momentum equation. 

The starting point is the momentum equation:

$$\begin{equation}
\rho \left( \frac{\partial \mathbf{v}}{\partial t} + \mathbf{v} \cdot \nabla \mathbf{v} \right) = \mathbf{j} \times \mathbf{B} - \nabla p + \mu \nabla^2 \mathbf{v}
\end{equation}$$

Projecting this in the parallel direction by multiplying with $F_0 \nabla\phi$, and assuming $F \approx F_0$ (constant), we analyze the terms.

### Time Derivative Term

The time derivative of the velocity $\mathbf{v} = \mathbf{v}_{pol} + v_{||} \mathbf{B}$ projected onto $F_0 \nabla\phi$ gives:

$$\begin{equation*}
F_0 \nabla\phi \cdot \rho \frac{\partial \mathbf{v}}{\partial t} \approx \rho \frac{\partial v_{||}}{\partial t}
\end{equation*}$$

### Advection Term

The advection term $\rho (\mathbf{v} \cdot \nabla \mathbf{v})$ projected onto $F_0 \nabla\phi$ is expanded as:

$$\begin{equation*}
F_0 \nabla\phi \cdot \rho (\mathbf{v} \cdot \nabla \mathbf{v}) = - R \rho \left[ u, v_{||} \right] + \frac{1}{R} \rho v_{||} \left[ \psi, v_{||} \right] - \frac{F_0}{R^2} \rho v_{||} \partial_3 v_{||}
\end{equation*}$$

### Force Terms

The $\mathbf{j} \times \mathbf{B}$ term and pressure gradient $\nabla p$ terms are projected as:

$$\begin{align*}
F_0 \nabla\phi \cdot (\mathbf{j} \times \mathbf{B}) &= - \frac{1}{2R^2} \frac{\partial}{\partial \phi} \left| \nabla \psi \right|^2 \\
- F_0 \nabla\phi \cdot \nabla p &= - \frac{\partial p}{\partial \phi}
\end{align*}$$

### Combined Equation

Combining these parts, we obtain the parallel momentum equation in the form used in JOREK:

$$\begin{equation*}\begin{split}
\rho \frac{\partial v_{\parallel}}{\partial t} = 
  & - R \rho \left[ u, v_{\parallel}\right] + \frac{1}{R} \rho v_{\parallel} \left[ \psi, v_{\parallel}\right] 
    - \frac{F_0}{R^2} \rho v_{\parallel} \partial_3 v_{\parallel}  \\
  & - \frac{1}{2R^2} \frac{1}{F_0} \frac{\partial}{\partial \phi} \left| \nabla \psi\right|^2 - \frac{1}{F_0} \frac{\partial p}{\partial \phi}
\end{split}\end{equation*}$$

This equation describes the evolution of the parallel velocity component in the reduced MHD limit.
