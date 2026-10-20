---
title: "Run with Sheath Heat-Flux BC"
nav_order: 2
parent: "Boundary Conditions"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# Sheath Heat-Flux Boundary Conditions

This page covers:

- [choosing where to apply the sheath boundary conditions](#choosing-the-boundaries);
- [the one-temperature sheath condition and its implementation](#one-temperature-models);
- [the electron and ion sheath conditions for two-temperature models](#two-temperature-models);
- [the treatment of grazing-incidence field lines](#treatment-at-grazing-angles).

The related [$v_\parallel$ boundary-condition smoothing](./vpar_smooth.md) can be used to smooth the Mach-one parallel-velocity condition near grazing angles.

## Choosing the boundaries

The sheath boundary conditions can be applied to selected boundary types in any grid by using the `bcs` structure. See [Choose Boundary Conditions](./choose_boundary_conditions.md) for the boundary-type indices and the complete selection rules.

For example, to apply the one-temperature sheath heat-flux condition on boundary type `i`, use:

```fortran
bc_natural_open      = .true. ! Enable boundary-integral calls
bcs(i)%dirichlet%T   = .false.
bcs(i)%natural%T     = .true.
gamma_stangeby       = 8.d0
```

Replace `i` with each boundary-type index on which the condition should be imposed. The settings must be written separately for every selected boundary type; `bcs(:)` cannot be used. For separate ion and electron temperatures, select `bcs(i)%natural%Ti` and `bcs(i)%natural%Te` instead, as described [below](#two-temperature-models).

## One-temperature models

This section applies to models 303, 333 and 500, and to model 600 when it is run with a single temperature.

In most boundary types, the temperature is fixed at the boundary by default (Dirichlet B.C.). Instead, you can impose the normal heat flux (to the boundary) to be

$$
\begin{align*}
  \mathbf{Q}\cdot\mathbf{n} &= \gamma_{sh} n_e k_B T_e \mathbf{v}_\parallel\cdot\mathbf{n} \\
\end{align*}
$$

which is given by P. Stangeby, *The Plasma Boundary of Magnetic Fusion Devices* (2000), Chapter 2, Equation 2.94. Typical values of $\gamma_{sh}$ are 7–8 for $T_e=T_i$.

### How the condition is implemented

The user-facing `gamma_stangeby` is converted internally into the coefficient $\gamma_{sh}^{JOREK}$ used in the natural boundary condition. The conversion follows from matching the heat flux represented by the JOREK energy equation to the Stangeby heat flux.

In JOREK, the implemented natural B.C. is 

$$
\begin{align*}
  -\kappa\nabla T \cdot\mathbf{n} &= (\gamma_{sh}^{JOREK}-1) n_e k_B T \mathbf{v}_\parallel\cdot\mathbf{n} \\
\end{align*}
$$

where $T=T_i+T_e$. The total kinetic and thermal MHD heat flux (see Ideal MHD or Goedbloed books) is

$$
\begin{align*}
  \mathbf{Q}\cdot\mathbf{n} &= -\frac{\kappa}{\gamma-1}\nabla T \cdot\mathbf{n} + \frac{\gamma}{\gamma-1}n_e k_B T \mathbf{v}_\parallel\cdot\mathbf{n} + \frac{\rho}{2}v^2 \mathbf{v}_\parallel\cdot\mathbf{n} \\
\end{align*}
$$

For $T_i=T_e$, one has $T=2T_e$. At the sheath entrance, $v^2=c_s^2=\gamma k_BT/m_i$, so that

$$
\begin{align*}
  \frac{\rho}{2}v^2 &= \frac{\gamma}{2}n_e k_BT .
\end{align*}
$$

Substituting the coded boundary condition into the MHD heat flux therefore gives

$$
\begin{align*}
  \mathbf{Q}\cdot\mathbf{n}
  &= \left[
       \frac{\gamma_{sh}^{JOREK}-1+\gamma}{\gamma-1}
       + \frac{\gamma}{2}
     \right]
     n_e k_BT\,\mathbf{v}_\parallel\cdot\mathbf{n} .
\end{align*}
$$

The Stangeby expression, written using the JOREK temperature, is

$$
\begin{align*}
  \mathbf{Q}\cdot\mathbf{n}
  &= \frac{\gamma_{sh}}{2}
     n_e k_BT\,\mathbf{v}_\parallel\cdot\mathbf{n} .
\end{align*}
$$

Equating the two expressions yields the conversion implemented in the code:

$$
\begin{align*}
  \gamma_{sh}^{JOREK} &= (\gamma-1) \left(\frac{\gamma_{sh}}{2}-1 -\frac{\gamma}{2} \right) \\
\end{align*}
$$

For example, $\gamma=5/3$ and `gamma_stangeby = 8.d0` give $\gamma_{sh}^{JOREK}=1.44$. The internal coefficient can also be supplied directly through `gamma_sheath`, but `gamma_stangeby` expresses the physical sheath-transmission coefficient and is normally the clearer input.

## Two-temperature models

This section applies to model 600 when separate ion and electron temperatures are evolved. 

Select the natural boundary condition independently for each temperature on every desired boundary type, following [Choose Boundary Conditions](./choose_boundary_conditions.md). For example:

```fortran
bc_natural_open      = .true.
bcs(i)%dirichlet%Te  = .false.
bcs(i)%natural%Te    = .true.
bcs(i)%dirichlet%Ti  = .false.
bcs(i)%natural%Ti    = .true.
gamma_e_stangeby     = 5.d0  ! Stangeby book end of section 2.8
gamma_i_stangeby     = 3.d0  ! Stangeby book end of section 2.8
```

The physical electron and ion sheath heat fluxes are written separately as

$$
\begin{align*}
  \mathbf{Q}_e\cdot\mathbf{n}
    &= \gamma_{sh,e} n_e k_BT_e\,\mathbf{v}_\parallel\cdot\mathbf{n}, \\
  \mathbf{Q}_i\cdot\mathbf{n}
    &= \gamma_{sh,i} n_i k_BT_i\,\mathbf{v}_\parallel\cdot\mathbf{n}.
\end{align*}
$$

JOREK implements the conductive parts through

$$
\begin{align*}
  -\kappa_e\nabla T_e\cdot\mathbf{n}
    &= (\gamma_{sh,e}^{JOREK}-1)n_e k_BT_e\,\mathbf{v}_\parallel\cdot\mathbf{n}, \\
  -\kappa_i\nabla T_i\cdot\mathbf{n}
    &= (\gamma_{sh,i}^{JOREK}-1)n_i k_BT_i\,\mathbf{v}_\parallel\cdot\mathbf{n}.
\end{align*}
$$

### Electron coefficient

The electron energy flux contains conductive and convective contributions,

$$
\begin{align*}
  \mathbf{Q}_e\cdot\mathbf{n}
  &= -\frac{\kappa_e}{\gamma-1}\nabla T_e\cdot\mathbf{n}
     + \frac{\gamma}{\gamma-1}n_e k_BT_e\,
       \mathbf{v}_\parallel\cdot\mathbf{n} \\
  &= \frac{\gamma_{sh,e}^{JOREK}-1+\gamma}{\gamma-1}
     n_e k_BT_e\,\mathbf{v}_\parallel\cdot\mathbf{n}.
\end{align*}
$$

Matching this expression to $\gamma_{sh,e}n_e k_BT_e\mathbf{v}_\parallel\cdot\mathbf{n}$ gives


$$
\begin{align*}
  \gamma_{sh,e}^{JOREK} &= (\gamma-1) \left(\gamma_{sh,e}-1 \right) \\
\end{align*}
$$

This is the conversion from `gamma_e_stangeby` to the internal `gamma_sheath_e`.

### Ion coefficient

The ion energy flux also carries the bulk kinetic energy:

$$
\begin{align*}
  \mathbf{Q}_i\cdot\mathbf{n}
  &= -\frac{\kappa_i}{\gamma-1}\nabla T_i\cdot\mathbf{n}
     + \frac{\gamma}{\gamma-1}n_i k_BT_i\,
       \mathbf{v}_\parallel\cdot\mathbf{n}
     + \frac{\rho}{2}v^2\mathbf{v}_\parallel\cdot\mathbf{n}.
\end{align*}
$$

Under the same $T_i=T_e$ sheath assumption, $v^2=c_s^2=\gamma k_B(T_i+T_e)/m_i=2\gamma k_BT_i/m_i$. Hence

$$
\begin{align*}
  \frac{\rho}{2}v^2 &= \gamma n_i k_BT_i,
\end{align*}
$$

and substitution of the coded boundary condition gives

$$
\begin{align*}
  \mathbf{Q}_i\cdot\mathbf{n}
  &= \left[
       \frac{\gamma_{sh,i}^{JOREK}-1+\gamma}{\gamma-1}
       + \gamma
     \right]
     n_i k_BT_i\,\mathbf{v}_\parallel\cdot\mathbf{n}.
\end{align*}
$$

Matching this expression to the physical ion sheath heat flux gives

$$
\begin{align*}
  \gamma_{sh,i}^{JOREK} &= (\gamma-1) \left(\gamma_{sh,i}-1 - \gamma\right) \\
\end{align*}
$$

This is the conversion from `gamma_i_stangeby` to the internal `gamma_sheath_i`. The extra $\gamma$ term compared with the electron expression comes from the ion kinetic-energy flux.

## Treatment at grazing angles

As shown in the Stangeby equation at the top of the page, the heat flux vanishes in situations where $v_\parallel=0$. These situations may occur when magnetic field lines are totally parallel to the wall (grazing angles). In those cases, artificial energy accumulation may occur, which is unphysical. For that, you can use a minimum heat flux, which is purely diffusive and is chosen as

$$
\begin{align*}
  \left.\mathbf{Q}\cdot\mathbf{n}\right|_{min} &= \gamma_{sh} n_e k_B T_e c_s \sin (\theta_{min}) = -\frac{\kappa}{\gamma-1}\nabla T \cdot\mathbf{n} \\
\end{align*}
$$

where the sound speed $c_s$ is calculated using the temperature directly (not $v_\parallel$). Note that the total heat flux now has an extra contribution due to this minimum heat flux:

$$
\begin{align*}
  \mathbf{Q}\cdot\mathbf{n} &= \gamma_{sh} n_e k_B T_e \mathbf{v}_\parallel\cdot\mathbf{n} + \left.\mathbf{Q}\cdot\mathbf{n}\right|_{min} \\
\end{align*}
$$

You can select the value of $\theta_{min}$ from the JOREK input file by setting

```fortran
min_sheath_angle = 1 ! Example of 1 degree (not in radians!)
```

**Important:** In order to calculate these fluxes, your grid resolution or $\kappa_\perp$ must be sufficiently large to resolve the numerical gradients imposed by the boundary condition. Note also that, by default, this treatment is applied to the particle flux:

$$
\begin{align*}
  \mathbf{\Gamma}\cdot\mathbf{n} &= n_i \mathbf{v}_\parallel\cdot\mathbf{n} + \left.\mathbf{\Gamma}\cdot\mathbf{n}\right|_{min} \\
\end{align*}
$$

where

$$
\begin{align*}
  \left.\mathbf{\Gamma}\cdot\mathbf{n}\right|_{min} &= n_i c_s \sin (\theta_{min}) = -D \nabla n_i\cdot\mathbf{n} \\
\end{align*}
$$

This approach was used in the paper below:

1. [F. J. Artola et al., 2021, *Plasma Physics and Controlled Fusion* **63**, 064004](https://iopscience.iop.org/article/10.1088/1361-6587/abf620)
