---
title: "Vorticity Expression"
nav_order: 8
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Derivation of the Vorticity Expression

$$\begin{equation*}
\mathbf{W}=\nabla\times\mathbf{v}=\nabla\times\left(-R^2\nabla u\times\nabla\phi+v_{||}B\right)
\end{equation*}$$

## First Part (depending on $u$)

$$\begin{equation*}\begin{split}
\mathbf{W}_u & =-\nabla\times\left(R^2\nabla u\times\nabla\phi\right) \\
  &= -\nabla R^2\times(\nabla u\times\nabla\phi)-R^2\nabla\times(\nabla u\times\nabla\phi) \\
  &= -2\nabla R\times(\partial_2 u\mathbf{a}^1-\partial_1 u\mathbf{a}^2)+R^2\Delta^*u\nabla\phi-\nabla_\text{pol}\partial_3 u \\
  &= +2R\partial_1 u\nabla\phi+R^2\Delta^*u\nabla\phi-\nabla_\text{pol}\partial_3 u \\
  &= \left(2R\partial_1 u + R^2\Delta^* u\right)\nabla\phi - \nabla_\text{pol}\partial_3 u  \\
  &= R^2\Delta_\text{pol} u\nabla\phi - \nabla_\text{pol}\partial_3 u
\end{split}\end{equation*}$$

When multiplied by $\cdot\nabla\phi$, one gets

$$\begin{equation*}
  \omega\equiv \mathbf{W}_u\cdot\nabla\phi=\Delta_\text{pol}u
\end{equation*}$$

or:

$$\begin{equation*} \begin{split}
\mathbf{W}_u\cdot\nabla\phi&=\left(\nabla \times\left(- R^2 \nabla u \times \nabla \phi \right)\right) \cdot \nabla \phi =
\nabla \cdot \left(- R^2 \left(\nabla u \times \nabla \phi \right) \times \nabla \phi \right) \\ 
&=\nabla \cdot \left(+ R^2 \left(\nabla \phi \cdot \nabla \phi \right) \nabla u -R^2 \left( \nabla u \cdot \nabla \phi \right) \nabla \phi \right) =\nabla \cdot\nabla_{pol} u
\end{split}
\end{equation*}$$

with diamagnetic flow:

$$\begin{equation*} \begin{split}
\mathbf{W}_u\cdot\nabla\phi&=\left(\nabla \times\left(- R^2 \nabla u \times \nabla \phi -R^2 \frac{\tau_{IC}}{\rho} \nabla p \times \nabla \phi \right)\right) \cdot \nabla \phi \\
&=
\nabla \cdot \left(- R^2 \left(\nabla u \times \nabla \phi +\frac{\tau_{IC}}{\rho} \nabla p \times \nabla \phi \right) \times \nabla \phi \right) \\ 
&=\nabla \cdot \left(+ R^2 \left(\nabla \phi \cdot \nabla \phi \right) \left( \nabla u +\frac{\tau_{IC}}{\rho} \nabla p \right) -R^2  \left( \left(\nabla u +\frac{\tau_{IC}}{\rho}\nabla p \right) \cdot \nabla \phi \right) \nabla \phi  \right) \\
&=\nabla \cdot \left( \nabla_{pol} u +\frac{\tau_{IC}}{\rho} \nabla_\text{pol} p \right) 
\end{split}
\end{equation*}$$

## Second Part (depending on $v_{||}$)

$$\begin{equation*}\begin{split}
\mathbf{W}_{v_{||}} &= \nabla\times(v_{||}\mathbf{B}) \\
  &= \nabla v_{||}\times\mathbf{B}+v_{||}\underbrace{\nabla\times\mathbf{B}}_{\equiv\mathbf{j}} \\
  &= \left(\nabla_\text{pol}v_{||}+\partial_3v_{||}\nabla\phi\right)\times\left(F\nabla\phi+\nabla\psi\times\nabla\phi\right)+v_{||}\mathbf{j} \\
  &= F\nabla_\text{pol}v_{||}\times\nabla\phi+\nabla_\text{pol}v_{||}\times(\nabla\psi\times\nabla\phi)+\partial_3 v_{||}F \nabla\phi\times\nabla\phi+\partial_3 v_{||}\nabla\phi\times(\nabla\psi\times\nabla\phi)+v_{||}\mathbf{j}
\end{split}\end{equation*}$$

When multiplied by $\cdot\nabla\phi$, one gets (taking into account the definition of the current in JOREK: $j:=\Delta^*\psi$):

$$\begin{equation*}\begin{split}
\mathbf{W}_{v_{||}} \cdot \nabla \phi &= \left( \nabla \times\left( v_{||}\mathbf{B} \right) \right) \cdot \nabla \phi
 = \nabla \cdot \left(v_{||} \mathbf{B} \times \nabla \phi \right) =\nabla \cdot \left( v_{||} \left(\nabla \psi \times \nabla \phi \right) \times \nabla \phi \right) \\
&= - \nabla \cdot \left(v_{||} \frac{1}{R^2} \nabla_\text{pol} \psi \right) \\
&= -  \frac{1}{R^2}\nabla v_{||} \cdot \nabla_\text{pol} \psi-v_{||} \nabla \cdot\left( \frac{1}{R^2}\nabla_\text{pol} \psi\right) 
=  -  \frac{1}{R^2}\left( \nabla v_{||} \cdot \nabla_\text{pol} \psi + v_{||} j \right)
\end{split}\end{equation*}$$

## Third Part (depending on $v_{neo}$)

This derivation is exactly the same as for the $\vec{E}\times\vec{B}$ one above, just replacing $u$ with $T$:

$$\begin{equation*}\begin{split}
\mathbf{W}_{neo} & =-\tau_{ic}\nabla\times\left(R^2\nabla T\times\nabla\phi\right) \\
  &= \tau_{ic}R^2\Delta_\text{pol} T\nabla\phi - \nabla_\text{pol}\partial_3 T
\end{split}\end{equation*}$$

When multiplied by $\cdot\nabla\phi$, one gets

$$\begin{equation*}
  \omega_{neo} \equiv \mathbf{W}_{neo}\cdot\nabla\phi = \tau_{ic}\Delta_\text{pol}T
\end{equation*}$$

or:

$$\begin{equation*}
  \mathbf{W}_{neo}\cdot\nabla\phi = \tau_{ic}\nabla\cdot\nabla_{pol} T
\end{equation*}$$

## Combined

$$\begin{equation*}\begin{split}
\mathbf{W} \cdot \nabla \phi &= \nabla \cdot \left(\nabla_\text{pol}u +\frac{\tau_{IC}}{\rho} \nabla_\text{pol} p - v_{||} \frac{1}{R^2} \nabla_\text{pol} \psi + \tau_{ic}\nabla_{pol} T\right)
\end{split}\end{equation*}$$
