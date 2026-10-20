---
title: "RMHD Model"
nav_order: 1
parent: "Reduced MHD models"
grand_parent: "Base Fluid Models"
layout: default
render_with_liquid: false
---

# JOREK Reduced MHD Models

Refer also to the [models](../base_fluid_models.md), [coordinates](../../coordinates.md), [notation](../../notation.md), and [normalization](../../normalization.md) pages.

**Note:** The present JOREK model assumes that $F\equiv F_0$ is spatially and temporaly constant. The derivation shown here is valid for general $F(R,Z,t)$ which is, however, independent of $\phi$ to guarantee $\nabla\cdot\mathbf{B}=0$.

## Definitions

### Magnetic Field

The magnetic field is not normalized. A general form of $\mathbf{B}$ that satisfies $\nabla\cdot\mathbf{B}=0$ is given by:

$$\begin{align*}
  \mathbf{B} &= \nabla\psi\times \nabla\phi + F~\nabla\phi = \frac{1}{R}(-\partial_1\psi~\mathbf{a}^2 + \partial_2\psi~\mathbf{a}^1)+F\mathbf{a}^3 \\
\end{align*}$$

where $F$ and $\psi$ are scalars and $\partial_3 F=0$ to guarantee $\nabla\cdot\mathbf{B}=0$. The physical components (with respect to normalized basis vectors) are:

$$\begin{align*}
  B_R &= +\frac{1}{R}\frac{\partial\psi}{\partial Z} \\
  B_Z &= -\frac{1}{R}\frac{\partial\psi}{\partial R} \\
  B_\phi &= +\frac{F}{R}
\end{align*}$$

The poloidal magnetic flux is connected with the vector potential by:

$$\begin{equation*}
\psi = R^2\mathbf{A}\cdot\nabla\phi = R A_\phi
\end{equation*}$$

Some additional remarks about the [magnetic vector potential can be found here](expressions_and_derivations/mag_vec_pot.md) (not relevant for the derivation of the reduced MHD equations).

### Current

The [normalization](../../normalization.md) of the plasma current is given by $\mathbf{j}_\text{SI}=\mathbf{j}/\mu_0$. The toroidal current density $j$ (it is not really a current density due to the factor R!) is defined as

$$\begin{equation*}
  j \equiv - R~\mathbf{j}\cdot\mathbf{e}_\phi = - \mathbf{j}\cdot\mathbf{a}_3 = \Delta^*\psi
\end{equation*}$$

[The derivation is sketched here](expressions_and_derivations/j_psi.md). The physical component $j_\phi$ is therefore given by

$$\begin{equation*}
  j_\phi = -\frac{j}{R} = -\frac{\Delta^*\psi}{R}
\end{equation*}$$

Although not represented by JOREK variables, R- and Z-components of the current exist as well (see [here](expressions_and_derivations/j_psi.md)):

$$\begin{align*}
  j_R &= +\frac{1}{R}\partial_2 F+\frac{1}{R^2}\partial_{1,3}\psi \\
  j_Z &= -\frac{1}{R}\partial_1 F+\frac{1}{R^2}\partial_{2,3}\psi
\end{align*}$$

### Velocity

Normalization (see also [normalization table](../../normalization.md)): $\mathbf{v}_\text{SI}=\mathbf{v}/\sqrt{\mu_0~\rho_0}$.

##### Standard Velocity

$$\begin{equation*}
\begin{split}
\mathbf{v} &= -R^2\nabla u \times \mathbf{a}^3 + v_{||}\mathbf{B} \\
&= R(\partial_1 u~\mathbf{a}^2 - \partial_2 u~\mathbf{a}^1) + v_{||}\left[\frac{1}{R}(-\partial_1\psi~\mathbf{a}^2 + \partial_2\psi~\mathbf{a}^1)+F\mathbf{a}^3\right] \\
&= \left(-R\partial_2 u+\frac{v_{||}}{R}\partial_2\psi\right)\mathbf{a}^1 + \left(R\partial_1 u-\frac{v_{||}}{R}\partial_1\psi\right)\mathbf{a}^2+F~v_{||}\mathbf{a}^3
\end{split}
\end{equation*}$$

The physical components (with respect to normalized basis vectors) are:

$$\begin{align*}
  v_R    &= -R\partial_2 u+\frac{v_{||}}{R}\partial_2\psi \\
  v_Z    &= +R\partial_1 u-\frac{v_{||}}{R}\partial_1\psi \\
  v_\phi &= \frac{F~v_{||}}{R}
\end{align*}$$

##### Ion Velocity with Diamagnetic Component

$$\begin{equation*}
  \begin{split}
\mathbf{v} &= -R^2\nabla u \times \mathbf{a}^3 - \frac{\tau_{IC} R^2}{\rho}\nabla p_i\times\mathbf{a}^3 + v_{||}\mathbf{B} \\
&= R(\partial_1 u~\mathbf{a}^2 - \partial_2 u~\mathbf{a}^1) + v_{||}\left[\frac{1}{R}(-\partial_1\psi~\mathbf{a}^2 + \partial_2\psi~\mathbf{a}^1)+F\mathbf{a}^3\right] -\frac{\tau_{IC}R}{\rho} \left[ -\partial_1 p_i~\mathbf{a}^2 + \partial_2 p_i~\mathbf{a}^1 \right] \\
&= \left(-R\partial_2 u+\frac{v_{||}}{R}\partial_2\psi - \frac{\tau_{IC}R}{\rho}\partial_2 p_i\right)\mathbf{a}^1 + \left(R\partial_1 u-\frac{v_{||}}{R}\partial_1\psi + \frac{\tau_{IC}R}{\rho}\partial_1 p_i\right)\mathbf{a}^2+F~v_{||}\mathbf{a}^3
\end{split}
\end{equation*}$$

(to be added soon...)

### Vorticity

The toroidal component of the total vorticity $\mathbf{W}$ is: 

$$\begin{equation*}\begin{split}
\mathbf{W} \cdot \nabla \phi &= \nabla \cdot \left(\nabla_\text{pol}u +\frac{\tau_{IC}}{\rho} \nabla_\text{pol} p - v_{||} \frac{1}{R^2} \nabla_\text{pol} \psi \right)
\end{split}\end{equation*}$$

(see [here](expressions_and_derivations/omega.md) for the derivation).

The vorticity is normalized by $\mathbf{W}_\text{SI}=\mathbf{W}/\sqrt{\mu_0~\rho_0}$.

In the code, variable number 4, which we denote $\omega$, corresponds to the following part of the total vorticity:

$$\begin{equation*}
  \omega \equiv \nabla \cdot \nabla_\text{pol}u = \Delta_\text{pol} u
\end{equation*}$$

### Electric Potential

The velocity stream function $u$ is connected to the electric potential $\Phi$ by
$$\begin{equation*}
  u\equiv \Phi/F_0
\end{equation*}$$
when $F\equiv F_0$ and the resistive term is small. [The derivation is sketched here](expressions_and_derivations/u_phi.md).

## Useful Expressions

| Expression | Description |
|---|---|
| $\nabla\times\nabla\phi\equiv 0$ | |
| $\nabla\cdot\nabla\phi\equiv 0$ | |
| $\nabla a\cdot\nabla\phi=\frac{1}{R^2}\partial_3 a$ | For arbitrary scalars $a$, $b$ |
| $\Delta_\text{pol} u=\Delta^* u+\frac{2}{R}\partial_1 u$ | |
| $\nabla a\times\nabla\phi = \frac{1}{R}(\partial_2 a\mathbf{a}_1 - \partial_1 a\mathbf{a}_2)$ | |
| $\nabla\phi\times(\nabla a\times\nabla\phi)=\frac{1}{R^2}\nabla_\text{pol}a$ | |
| $(\nabla a\times\nabla\phi)\times(\nabla b\times\nabla\phi) = \frac{1}{R}[a,b]\nabla\phi$ | |
| $\nabla\times(\nabla a\times\nabla\phi)=-\Delta^* a\nabla\phi + \frac{1}{R^2}\nabla_\text{pol}\partial_3 a$ | |
| $\mathbf{B}\cdot\nabla a=\frac{F}{R^2}\partial_3 a+\frac{1}{R}[a,\psi]$ | |
| $\nabla a\cdot\nabla b = \nabla_\text{pol}a\cdot\nabla_\text{pol}b+\frac{1}{R^2}\partial_3 a\partial_3 b$ | |
| $\nabla\phi\cdot(\nabla C\times \nabla D)=\frac{1}{R}[C,D]$ | For arbitrary vectors $\mathbf{C}$, $\mathbf{D}$ |
| $\mathbf{B}\cdot(\mathbf{C}\times\mathbf{D})=\frac{F}{R}[C,D]-\frac{C}{R}[\nabla\psi,D]+\frac{D}{R}[\nabla\psi,C]$ | |
| $\mathbf{j} \equiv \nabla F\times\nabla\phi + \frac{1}{R^2}\partial_3(\nabla_\text{pol}\psi)-\Delta^*\psi~\nabla\phi$ | [see here](expressions_and_derivations/j_psi.md) |
| $j \equiv -\mathbf{j}\cdot\mathbf{a}_3=\Delta^*\psi$ | [see here](expressions_and_derivations/j_psi.md) |
| $\mathbf{j}\times\mathbf{B}= -\frac{F}{R^2}\nabla_\text{pol}F+\frac{1}{R}\nabla\phi[F,\psi]-\frac{\Delta^*\psi}{R^2}\nabla_\text{pol}\psi+\frac{1}{R^2}\nabla_\text{pol}\partial_3\psi\times\mathbf{B}$ | [see here](expressions_and_derivations/j_cross_b.md) |
| $\mathbf{v}\times\mathbf{B} = \frac{1}{R}\mathbf{a}_3 [\psi,u]+F\nabla_\text{pol} u$ | [see here](expressions_and_derivations/v_cross_b.md) |

## Equations

### Continuity Equation

$$\begin{equation}
  \partial_t\rho=-\nabla\cdot(\rho\mathbf{v})+\nabla\cdot(D\nabla\rho) + S_\rho
\end{equation}$$

[The derivation in reduced MHD form is given here](expressions_and_derivations/continuity_eq.md).

### Poloidal Flux Equation

$$\begin{equation}
  \partial_t\psi=R[\psi,u]+\eta(T)(j-j_0)-F_0\partial_\phi u
\end{equation}$$

where $j_0$ is implemented acting as a current source. [The derivation is sketched here](expressions_and_derivations/psi_eq.md).

### Perpendicular Momentum Equation

The perpendicular projection of the momentum equation is obtained by applying the operator $\nabla\phi\cdot\nabla\times(R^2...)$ onto the momentum equation.

##### Standard Equation

$$\begin{equation*}\begin{split}
\nabla\cdot\left[R^2\rho\nabla_{pol}\left(\frac{\partial U}{\partial t}\right)\right] = 
  &   \frac{1}{R} \left[R^4\rho\omega,U\right]
    - \frac{1}{2R} \left[R^2\rho,R^2\left|\nabla_{pol}U\right|^2\right] \\ 
  & - \frac{1}{R} \left[R^2,p\right] + \frac{1}{R} \left[\psi,j\right]
    - \frac{F_{o}}{R^2}\frac{\partial j}{\partial\phi} + \mu \nabla^2\omega
\end{split}\end{equation*}$$

##### Equation with Diamagnetic Components

This equation is obtained by replacing $\mathbf{v}_{pol} = R^2 \nabla\phi\times\nabla U$ with 
$$\mathbf{v}_{pol} = R^2\nabla\phi\times\nabla U + \frac{\tau_{IC}R^2}{\rho}\nabla\phi\times\nabla P_i$$
and taking out the terms that cancel out with the gyro-viscous tensor (see gyro-viscous cancellation, [Schnack, Phys. Plasmas 13, 058103, 2006]).

$$\begin{equation*}\begin{split}
\nabla\cdot\left[R^2\rho\nabla_{pol}\left(\frac{\partial U}{\partial t}\right)\right] = 
  &   \frac{1}{R} \left[R^4\rho W,U\right]
    - \frac{1}{2R} \left[R^2\rho,R^2\left|\nabla_{pol}U\right|^2\right] \\ 
  & - \frac{1}{R} \left[R^4\tau_{IC}\nabla^{2}_{pol}P_{i},U\right]
    + \frac{1}{R} \left[R^4\frac{\tau_{IC}}{\rho}\nabla_{pol}\rho\cdot\nabla_{pol}P_{i},U\right] \\
  & + \tau_{IC}R^3\left[W,P_i\right]
    + \tau_{IC}R^2 \nabla\cdot\left(\frac{\partial P_i}{\partial Z}\nabla_{pol}U\right) \\ 
  & - \tau_{IC}R^3 \left[\partial_{xy}U  \left(\partial_{xx}P_i-\partial_{yy}P_i\right)
                        -\partial_{xy}P_i\left(\partial_{xx}U  -\partial_{yy}U  \right)\right] \\
  & - \frac{1}{R} \left[R^{2},p\right] + \frac{1}{R} \left[\psi,j\right]
    - \frac{F_{o}}{R^2}\frac{\partial j}{\partial\phi}
    + \mu \nabla^{2}\left(\omega + \omega^*\right)
\end{split}\end{equation*}$$

where $\omega^*$ is the diamagnetic part of the vorticity ([see definition here](expressions_and_derivations/omega_star.md)).
Note that it is very important to include $\omega^*$ in the equation. In the code, we choose to define the total vorticity as

$$\begin{equation*}
  W = \omega + \omega^*
\end{equation*}$$

Note that introducing this new variable in the above equation means we have to extract $\omega^*$ from the term $\left[R^4\rho\omega,U\right]$.
Ideally, we would also need to extract $\omega^*$ from the term $\tau_{IC}R^4[\omega,P_i]$, but this would result in a term of order $\tau_{IC}^2$, which is considered small enough to be ignored. The resulting equation is (with the new definition of the total vorticity $W$):

$$\begin{equation*}\begin{split}
\nabla\cdot\left[R^2\rho\nabla_{pol}\left(\frac{\partial U}{\partial t}\right)\right] = 
  &   \frac{1}{R} \left[R^4\rho W,U\right]
    - \frac{1}{2R} \left[R^2\rho,R^2\left|\nabla_{pol}U\right|^2\right] \\ 
  & - \frac{1}{R} \left[R^4\tau_{IC}\nabla^{2}_{pol}P_{i},U\right]
    + \frac{1}{R} \left[R^4\frac{\tau_{IC}}{\rho}\nabla_{pol}\rho\cdot\nabla_{pol}P_{i},U\right] \\
  & + \tau_{IC}R^3\left[W,P_i\right]
    + \tau_{IC}R^2 \nabla\cdot\left(\frac{\partial P_i}{\partial Z}\nabla_{pol}U\right) \\ 
  & - \tau_{IC}R^3 \left[\partial_{xy}U  \left(\partial_{xx}P_i-\partial_{yy}P_i\right)
                        -\partial_{xy}P_i\left(\partial_{xx}U  -\partial_{yy}U  \right)\right] \\
  & - \frac{1}{R} \left[R^{2},p\right] + \frac{1}{R} \left[\psi,j\right]
    - \frac{F_{o}}{R^2}\frac{\partial j}{\partial\phi}
    + \mu \nabla^{2}W
\end{split}\end{equation*}$$

### Parallel Momentum Equation

The parallel projection of the momentum equation is obtained by applying the operator $\nabla\phi\cdot (...)$ onto the momentum equation.

[The derivation can be found here.](expressions_and_derivations/velocity_equation.md)

##### Standard Equation

$$\begin{equation*}\begin{split}
\rho \frac{\partial v_{\parallel}}{\partial t} = 
  & - R \rho \left[ u, v_{\parallel}\right] + \frac{1}{R} \rho v_{\parallel} \left[ \psi, v_{\parallel}\right] 
    - \frac{F_0}{R^2} \rho v_{\parallel} \partial_3 v_{\parallel}  \\
  & - \frac{1}{2R^2} \frac{1}{F_0} \frac{\partial}{\partial \phi} \left| \nabla \psi\right|^2 - \frac{1}{F_0} \frac{\partial p}{\partial \phi}
\end{split}\end{equation*}$$

### Energy Equation
*(To be filled soon...)*

### Current Definition Equation
*(To be filled soon...)*

### Vorticity Definition Equation
*(To be filled soon...)*
