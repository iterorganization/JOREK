---
title: "Check Energy Conservation"
nav_order: 2
parent: "Diagnostics and post-processing"
layout: default
render_with_liquid: false
---

# Checking Energy Conservation

Live-data diagnostics can be used to monitor how well JOREK conserves energy during a simulation.

> **Important:** The energy-conservation diagnostic must be adapted to each model. Adding physical terms can introduce new energy channels, boundary fluxes, sources, or sinks. A balance implemented for one model must not be assumed complete for another.

## Total Energy Conservation

For a basic visco-resistive MHD model, the local total-energy balance can be written as

$$
\partial_t w+\nabla\cdot\boldsymbol{\Gamma}=\tau_{\mathrm{nc}},
\tag{1}
$$

where, omitting factors of $\mu_0$, the energy density and energy flux are

$$
\begin{aligned}
w
&\equiv
\frac{\rho}{2}\lvert\boldsymbol{v}\rvert^2
+\frac{p}{\gamma-1}
+\frac{\lvert\boldsymbol{B}\rvert^2}{2}, \\
\boldsymbol{\Gamma}
&\equiv
\left[
\frac{\rho}{2}\lvert\boldsymbol{v}\rvert^2
+\frac{\gamma}{\gamma-1}p
\right]\boldsymbol{v}
+\frac{\boldsymbol{q}}{\gamma-1}
+\boldsymbol{E}\times\boldsymbol{B}.
\end{aligned}
\tag{2}
$$

The right-hand side of Eq. (1), $\tau_{\mathrm{nc}}$, contains non-conservative terms such as heating sources and radiation sinks. Dissipative resistive and viscous terms should not appear in $\tau_{\mathrm{nc}}$ when the dissipated magnetic or kinetic energy is converted consistently into thermal energy.

For example, if $\eta\neq\eta_{\mathrm{ohmic}}$, part of the magnetic-energy loss is not restored as thermal energy by Ohmic heating. The mismatch therefore contributes to the non-conservative balance.

> **Reference for the derivation:** The standard total-energy equation follows by combining the MHD mass, momentum, internal-energy, and induction equations. A general derivation and discussion of the underlying MHD conservation laws can be found in J. P. Goedbloed and S. Poedts, [*Principles of Magnetohydrodynamics: With Applications to Laboratory and Astrophysical Plasmas*](https://www.cambridge.org/core/books/principles-of-magnetohydrodynamics/847AF12C94451B41D8F71C1F11EB308A), Cambridge University Press (2004). This page focuses on the JOREK diagnostic and its model-dependent energy accounting rather than repeating the general derivation.

Integrating Eq. (1) over the JOREK domain gives

$$
\partial_t E_{\mathrm{tot}}
=\text{boundary fluxes}
+\text{non-conservative/source terms},
\tag{3}
$$

with

$$
\begin{aligned}
E_{\mathrm{tot}}
&\equiv \int_V w\,dV, \\
\text{boundary fluxes}
&\equiv
-\oint_S\boldsymbol{\Gamma}\cdot\boldsymbol{n}\,dS, \\
\text{non-conservative/source terms}
&\equiv
\int_V\tau_{\mathrm{nc}}\,dV.
\end{aligned}
\tag{4}
$$

### Plotting the Total-Energy Balance

Plot the left- and right-hand sides of Eq. (3) with

```bash
jorek_folder/util/plot_live_data.sh -q energy_conservation
```

For good energy conservation, the two curves should overlap within the expected numerical accuracy. A persistent difference represents an unaccounted energy channel or a numerical conservation error.

### Integrated Energies and Time Derivatives

Plot the different energy components and their time derivatives with

```bash
jorek_folder/util/plot_live_data.sh -q integrated_energies
jorek_folder/util/plot_live_data.sh -q dEdt
```

The plotted labels correspond to

$$
\begin{aligned}
\text{Total energy}
&\equiv E_{\mathrm{tot}}, \\
\text{Magnetic}
&\equiv
\int_V\frac{\lvert\boldsymbol{B}\rvert^2}{2}\,dV, \\
\text{Thermal energy}
&\equiv
\int_V\frac{p}{\gamma-1}\,dV, \\
\text{Kinetic parallel}
&\equiv
\int_V\frac{\rho}{2}
\lvert\boldsymbol{v}_{\parallel}\rvert^2\,dV, \\
\text{Kinetic perpendicular}
&\equiv
\int_V\frac{\rho}{2}
\lvert\boldsymbol{v}_u\rvert^2\,dV,
\end{aligned}
\tag{5}
$$

where

$$
\boldsymbol{v}_u=-R^2\nabla u\times\nabla\phi.
$$

### Boundary Fluxes

Plot the individual boundary-flux contributions with

```bash
jorek_folder/util/plot_live_data.sh -q bnd_fluxes
```

The labels correspond to

$$
\begin{aligned}
\text{p vn}
&\equiv
\int_S\frac{\gamma}{\gamma-1}
p\,\boldsymbol{v}\cdot\boldsymbol{n}\,dS, \\
\text{kinpar-flux}
&\equiv
\int_S\frac{\rho}{2}\lvert\boldsymbol{v}\rvert^2
\boldsymbol{v}\cdot\boldsymbol{n}\,dS, \\
\text{qn-par}
&\equiv
\int_S
\frac{\boldsymbol{q}_{\parallel}\cdot\boldsymbol{n}}
{\gamma-1}\,dS, \\
\text{qn-perp}
&\equiv
\int_S
\frac{\boldsymbol{q}_{\perp}\cdot\boldsymbol{n}}
{\gamma-1}\,dS.
\end{aligned}
\tag{6}
$$

If sheath boundary conditions are imposed on the entire domain boundary, not only on the divertor, the sum of these four fluxes corresponds to the energy flowing through the sheath according to the Stangeby expression. See [Running with Sheath Boundary Conditions for the Heat Flux](../howto/physics_options/boundary_conditions/sheath_heatflux_bc.html).

The Poynting flux is typically zero in fixed-boundary simulations. Its contribution can be inspected with the magnetic-energy balance described below.

### Non-Conservative and Source Terms

For the RMHD diagnostic described here, the non-conservative and source terms are

$$
\text{non-conservative/source terms}
=\int_V
\left[
\frac{S_T}{\gamma-1}
+S_{\mathrm{mag}}
+(\eta_{\mathrm{ohmic}}-\eta)
\lvert\boldsymbol{J}\rvert^2
-\boldsymbol{\nu}_{\parallel}
\cdot\boldsymbol{v}_{\parallel}
\right]\,dV.
\tag{7}
$$

Here:

- $S_T=S_{T_i}+S_{T_e}$ is the combined ion and electron heating source used in the temperature equations.
- $S_{\mathrm{mag}}$ is the magnetic-energy source associated with the current source.
- $(\eta_{\mathrm{ohmic}}-\eta)\lvert\boldsymbol{J}\rvert^2$ accounts for a mismatch between resistive magnetic-energy dissipation and Ohmic heating.
- The final term represents energy dissipated through parallel viscosity.

Plot these contributions with

```bash
jorek_folder/util/plot_live_data.sh -q dissipative_terms
```

All these quantities can also be recovered through [`jorek2_postproc`](../howto/introduction_to_jorek_diagnostics.html#jorek2_postproc) using the `expressions_int` command. Within `jorek2_postproc`:

```text
help expressions_int
expressions_int
```

The first command displays usage information. The second lists the available integrated expressions. You can get all of them and more with the `zeroD_quantities` command.

## Magnetic Energy Conservation

Let $W\equiv E_{\mathrm{mag}}$ denote the magnetic component of the total energy already defined in Eq. (5). In this section, $\mu_0=1$ is used for compactness.

> **Model dependence:** The derivation below gives the reduced-MHD balance implemented by this diagnostic. Its terms must be checked against the equations and physical effects enabled in the selected model.

### General Balance

The time evolution of the magnetic energy is

$$
\partial_t W
=
\int
\mathbf{B}\cdot\partial_t\mathbf{B}
\,dV
=
\int
\mathbf{B}\cdot\nabla\times\partial_t\mathbf{A}
\,dV.
\tag{8}
$$

Using the identity below with $\boldsymbol{C}=\partial_t\boldsymbol{A}$,

$$
\mathbf{B}\cdot\nabla\times\mathbf{C}
=
\nabla\cdot(\mathbf{C}\times\mathbf{B})
+
\mathbf{C}\cdot\nabla\times\mathbf{B},
$$

together with Ampère's law, we obtain

$$
\partial_t W
=
\int
\partial_t\mathbf{A}\cdot\mathbf{J}
\,dV
+
\int
\nabla\cdot
\left(
\partial_t\mathbf{A}\times\mathbf{B}
\right)
\,dV.
\tag{9}
$$

Using Gauss' divergence theorem,

$$
\partial_t W
=
\int
\partial_t\mathbf{A}\cdot\mathbf{J}
\,dV
+
\oint
\partial_t\mathbf{A}\times\mathbf{B}
\cdot d\mathbf{S}.
\tag{10}
$$

The surface term must be retained for the finite JOREK domain. It vanishes only when imposed boundary conditions make it zero, or when the boundary is taken sufficiently far away.

Writing the latter form in terms of the electric field,

$$
\mathbf{E}
=
-\partial_t\mathbf{A}
-
\nabla\Phi,
$$

gives

$$
\partial_t W
=
-
\int
\mathbf{E}\cdot\mathbf{J}
\,dV
-
\oint
\mathbf{E}\times\mathbf{B}
\cdot d\mathbf{S}
-
\int
\nabla\Phi\cdot\mathbf{J}
\,dV
-
\oint
\nabla\Phi\times\mathbf{B}
\cdot d\mathbf{S}.
\tag{11}
$$

The last two terms cancel. This follows by integrating the volume term by parts and using

$$
\nabla\Phi\cdot\nabla\times\mathbf{B}
=
-\nabla\cdot
\left(
\nabla\Phi\times\mathbf{B}
\right).
$$

We therefore obtain the standard form

$$
\boxed{
\partial_t W
=
-
\int
\mathbf{E}\cdot\mathbf{J}
\,dV
-
\oint
\mathbf{E}\times\mathbf{B}
\cdot d\mathbf{S}
}
\tag{12}
$$

The first term in Eq. (12) describes magnetic-energy exchange with the plasma through $\boldsymbol{E}\cdot\boldsymbol{J}$. The second is the outward Poynting flux through the boundary.

### Reduced-MHD Specialization

In reduced MHD,

$$
\partial_t\mathbf{A}
=
\partial_t A_\phi\,\hat{\mathbf e}_\phi
=
\partial_t\left(\frac{\psi}{R}\right)
\hat{\mathbf e}_\phi.
$$

Also using the JOREK variable

$$
j = -J_\phi R,
$$

the volume term of Eq. (10) becomes

$$
\mathcal{V}_{\mathrm{mag}}
\equiv
\int
\partial_t\mathbf{A}\cdot\mathbf{J}
\,dV
=
-
\int
\frac{1}{R^2}
j\,\partial_t\psi
\,dV.
\tag{13}
$$

The integrand can be written in terms of the evolved variables using the reduced-MHD poloidal-flux equation,

$$
j\,\partial_t\psi
=
j
\left(
R[\psi,u]
+
\eta(j-j_S)
-
F\frac{\partial u}{\partial\phi}
\right).
\tag{14}
$$

Here, $j_S$ is the current source and $F=RB_\phi$. See the [derivation of the poloidal-flux equation](../physics/base_fluid_models/RMHD/expressions_and_derivations/psi_eq.html) for the sign and variable conventions.

This gives

$$
\mathcal{V}_{\mathrm{mag}}
=
\int
\partial_t\mathbf{A}\cdot\mathbf{J}
\,dV
=
-
\int
[\psi,u]
\frac{j}{R}
\,dV
-
\int
\eta
\left(
\frac{j}{R}
\right)^2
\,dV
+
\int
\eta
\frac{jj_S}{R^2}
\,dV
+
\int
\frac{Fj}{R^2}
\frac{\partial u}{\partial\phi}
\,dV.
\tag{15}
$$

The boundary flux is

$$
\partial_t\mathbf{A}\times\mathbf{B}
=
\partial_t\psi
\nabla\phi
\times
\left(
\nabla\psi\times\nabla\phi
\right)
=
\frac{\partial_t\psi}{R^2}
\nabla_{\mathrm{pol}}\psi.
\tag{16}
$$

Using the flux equation, the boundary flux becomes

$$
\partial_t\mathbf{A}\times\mathbf{B}
=
\frac{1}{R}
[\psi,u]
\nabla_{\mathrm{pol}}\psi
+
\frac{1}{R^2}
\left(
\eta(j-j_S)
-
F\frac{\partial u}{\partial\phi}
\right)
\nabla_{\mathrm{pol}}\psi.
\tag{17}
$$

The first term can be shown to be

$$
-\mathbf{v}_{\perp}^{\mathrm{pol}}
B_{\mathrm{pol}}^2.
$$

Combining the volume contribution $\mathcal{V}_{\mathrm{mag}}$ from Eq. (15) with the boundary contribution from Eq. (16) gives the magnetic-energy balance without repeating the expanded volume terms:

$$
\boxed{
\partial_t W
=
\mathcal{V}_{\mathrm{mag}}
+
\oint
\frac{\partial_t\psi}{R^2}
\nabla_{\mathrm{pol}}\psi
\cdot\mathbf{n}
\,dS
}
\tag{18}
$$

### Diagnostic Labels

By running

```bash
jorek_folder/util/plot_live_data.sh -q mag_energy_balance
```

you can compare the different terms of Eq. (18) as functions of time.

The diagnostic names map onto the derivation as follows:

| Live-data label | Term in the magnetic-energy balance |
|:---|:---|
| `magSource` | The third term in Eq. (15), representing magnetic-energy input from the current source $j_S$. |
| `Ohmic` | The second term in Eq. (15), representing resistive loss of magnetic energy. |
| `JxB.v` | The sum of the first and fourth terms in Eq. (15), representing magnetic-to-kinetic energy transfer. |
| `Poynting` | The surface integral in Eq. (18), representing magnetic energy crossing the domain boundary. |
| `dWmagdt` | The left-hand side of Eq. (18), $\partial_t W$. |
| `sum all...` | The sum `magSource + Ohmic + JxB.v + Poynting`, i.e. the complete right-hand side of Eq. (18). |

For a closed balance, `dWmagdt` and `sum all...` should agree within the temporal, spatial, and diagnostic-integration accuracy. A systematic difference indicates a missing physical term, an inconsistent coefficient, an unaccounted boundary contribution, or insufficient numerical accuracy.
