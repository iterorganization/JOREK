---
title: "Stream Function and Electric Potential Expression"
nav_order: 9
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of Stream Function and Electric Potential

The connection between $u$ and $\Phi$ is derived in the following. The starting point is Faraday's law

$$\begin{equation}
  \partial_t\mathbf{B}=-\nabla\times\mathbf{E},
\end{equation}$$

which can also be expressed in the magnetic vector potential:

$$\begin{equation}\label{dAdt}
  \partial_t\mathbf{A}=-\mathbf{E}-\nabla\Phi+\alpha_L\nabla\phi.
\end{equation}$$

Here, the loop voltage $\alpha_L$ and the electric potential $\Phi$ are introduced as integration constants by removing the curl. From the requirement $\nabla\times(\alpha_L\nabla\phi)=0$ it follows, that $\nabla_\text{pol}\alpha_L$ must vanish. The electric field is given by

$$\begin{equation}\label{E}
  \mathbf{E}=-\mathbf{v}\times\mathbf{B}+\eta\mathbf{j}=-F\nabla_\text{pol}u - [\Psi,u]\mathbf{e}_\phi + \eta\mathbf{j},
\end{equation}$$

When plugging Eq. $\eqref{E}$ into Eq. $\eqref{dAdt}$, we obtain

$$\begin{equation}
  \partial_t\mathbf{A}=F\nabla_\text{pol}u + [\Psi,u]\mathbf{e}_\phi - \eta\mathbf{j} - \nabla\Phi + \alpha_L\nabla\phi.
\end{equation}$$

The loop voltage term is combined with the current in JOREK by introducing $\mathbf{j}_0$:

$$\begin{equation}\label{dAdt9}
  \partial_t\mathbf{A}=F\nabla_\text{pol}u + [\Psi,u]\mathbf{e}_\phi - \eta(\mathbf{j}-\mathbf{j}_0) - \nabla\Phi.
\end{equation}$$

To be fully consistent, $j_0$ needs to be $\frac{\alpha_L}{\eta}\nabla\phi$ with spatially constant $\alpha_L$.

## When taking $F\equiv F_0$ and neglecting $\eta(\mathbf{j}-\mathbf{j}_0)$:

We get $\partial_t\mathbf{A}=\partial_t\psi\nabla\phi$. Applying $\nabla\phi\times(\dots)$ to Eq. $\eqref{dAdt9}$, we obtain

$$\begin{equation}
  F_0\nabla\phi\times\nabla_\text{pol}u=\nabla\phi\times\nabla_\text{pol}\Phi
\end{equation}$$

which means that electric potential $\Phi$ and velocity stream function $u$ are connected by

$$\begin{equation}
  F_0 u \equiv \Phi.
\end{equation}$$

(In principle, we have an integration constant $c\phi$ in this equation which is chosen to be zero)

We can then calculate the local electric field from $E = -\nabla \Phi - \partial_t A$.

$$\begin{equation}
  \mathbf{E} = \left(- F_0 u_R, - F_0 u_Z, - F_0 u_\phi/R - \partial_t \psi/R\right),
\end{equation}$$

where subscripts denote differentiation and the vector components $(R,Z,\phi)$ are written down separately.
