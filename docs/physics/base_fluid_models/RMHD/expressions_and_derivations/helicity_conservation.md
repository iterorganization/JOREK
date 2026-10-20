---
title: "Helicity Conservation"
nav_order: 17
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Helicity Conservation in JOREK

Magnetic helicity is defined as:

$$
\begin{equation*}
H=\mathbf{A}\cdot\mathbf{B}.
\end{equation*}
$$

It can be shown that in ideal MHD, the volumic integral of $H$ is a conserved quantity if the components of $\mathbf{B}$ and $\mathbf{v}$ normal to the boundary vanish (see for example, Section 2.2 of the book *Nonlinear Magnetohydrodynamics* by D. Biskamp).

In resistive MHD, $H$ is not exactly conserved but it is dissipated at a rate proportional to the plasma resistance (see Section 7.3.1 of Biskamp). During a fast MHD relaxation (like the thermal quench phase of a disruption, for example), it is expected that $H$ will be dissipated much more slowly than magnetic energy, which sets constraints to the global dynamics (see for example [this 2019 PPCF paper](https://iopscience.iop.org/article/10.1088/1361-6587/aaf293) by A. Boozer).

In the reduced MHD version of JOREK, an acceptable expression for the vector potential is:

$$
\begin{equation*}
\mathbf{A} = \psi \nabla \phi - F_0 \frac{Z}{R} \nabla R.
\end{equation*}
$$

Indeed,

$$
\begin{equation*}
\nabla \times \mathbf{A} = \nabla\psi\times \nabla\phi + F_0 \nabla\phi
\end{equation*}
$$

which corresponds to the expression for $\mathbf{B}$ (see [JOREK Reduced MHD Models](/JOREK/physics/base_fluid_models/RMHD/rmhd_model.html)), and

$$
\begin{equation*}
\partial_t \mathbf{A} = \partial_t \psi \nabla\phi
\end{equation*}
$$

which corresponds to the inductive part of the electric field.

The corresponding expression for the helicity is:

$$
\begin{equation*}
  \begin{split}
H &= \mathbf{A} \cdot \mathbf{B} \\
  &= (\psi \nabla \phi - F_0 \frac{Z}{R} \nabla R) \cdot (\nabla\psi\times \nabla\phi + F_0 \nabla\phi) \\
  &= F_0 \frac{\psi}{R^2} - F_0 \frac{Z}{R^2} \partial_Z \psi \\
  &= F_0 \frac{\psi}{R^2} - \frac{F_0}{R^2} [\partial_Z (Z \psi) - \psi] \\
  &= 2 F_0 \frac{\psi}{R^2} - \frac{F_0}{R^2} \partial_Z (Z \psi).
  \end{split}
\end{equation*}
$$

Let us now consider $\frac{d}{dt}\int H dV$, where the integral is a volumic integral over the computation domain (the latter being fixed in time):

$$
\begin{equation*}
  \begin{split}
\frac{d}{dt}\int H dV &= \int \partial_t H dV \\
  &= \int [2 F_0 \frac{\partial_t \psi}{R^2} - \frac{F_0}{R^2} \partial_Z (Z \partial_t \psi)] dV \\
  &= \int [2 F_0 \frac{\partial_t \psi}{R} - \frac{F_0}{R} \partial_Z (Z \partial_t \psi)] d\phi dS \\
  \end{split}
\end{equation*}
$$

where $dS$ is an infinitesimal surface element in the poloidal plane.

The second term can be readily integrated along $Z$ to give boundary terms which vanish if we have $\partial_t \psi = 0$ at the boundary, which is the case for a fixed boundary simulation.

To calculate the first term, we make use of the poloidal flux equation:

$$
\begin{equation*}
  \partial_t\psi = R[\psi,u] - F_0 \partial_\phi u + \eta j - R^2 \nabla \cdot \left( \frac{\eta_{num}}{R} \nabla_{\text{pol}} j \right)
\end{equation*}
$$

The term in $\partial_\phi u$ will vanish when integrated over $\phi$. As for the Poisson bracket:

$$
\begin{equation*}
  \begin{split}
  [\psi,u] &= \partial_R \psi \partial_Z u - \partial_Z \psi \partial_R u \\
           &= \partial_Z (u \partial_R \psi) - \partial_R (u \partial_Z \psi)
  \end{split}
\end{equation*}
$$

and both of these terms will vanish by integration over $R$ or $Z$ if we have $u=0$ at the boundary, which is the case for a fixed boundary simulation, i.e., the electric potential is constant on an ideally conducting wall corresponding to a vanishing normal velocity component into the wall.

The hyper-resistivity term becomes a boundary integral due to the divergence theorem:

$$
\begin{equation*}
 \int - \nabla \cdot \left( \frac{\eta_{num}}{R} \nabla_{\text{pol}} j \right)dV = -\oint \frac{\eta_{num}}{R} \mathbf{n}\cdot\nabla_{\text{pol}} j dS_{bnd}
\end{equation*}
$$

where $\mathbf{n}$ is a unit vector normal to the JOREK boundary (pointing outwards the domain). Finally the time derivative of the total helicity is

$$
\begin{equation*}
  \begin{split}
\frac{d}{dt}\int H dV &= 2 F_0 \int \frac{\eta j}{R} d\phi dS -2F_0\oint \frac{\eta_{num}}{R} \mathbf{n}\cdot\nabla_{\text{pol}} j dS_{bnd} \\
  \end{split}
\end{equation*}
$$

We therefore recover the classical result that helicity is dissipated at a resistive rate.
