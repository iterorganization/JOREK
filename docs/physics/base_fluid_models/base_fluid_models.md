---
title: "Base Fluid Models"
nav_order: 5
parent: "Physics Models"
has_children: true
nav_fold: true
layout: default
render_with_liquid: false
---

Andres working on this now

# JOREK Physics Models

Note that on top of the models described below, **multiple switches are available to activate/deactive different terms and boundary conditions in the equations**. For this, go to [howto/physics_options](/howto/physics_options.md).

## [Reduced MHD](reduced_mhd.md): 

__model 600__ is the workhorse for reduced MHD as it allows users to combine the different physics ($v_\parallel$, separate $T_e$ and $T_i$, fluid neutrals, fluid impurities, ...) they wish to keep in their simulations at [hard-coded parameters](/JOREK/compiling/compilation_and_libraries/hard-coded_parameters.html) time. The variables that are always present in __model 600__ are:

| Variable | Symbol | variable number | 
| --- | --- | --- |
| Poloidal magnetic flux: | $\psi$ | 1 |
| Velocity stream function: | $u$ | 2 |
| Toroidal current density: | $zj$ | 3 |
| Vorticity: | $\omega$ | 4 |
| Mass density: | $\rho$ | 5 |
| Temperature: | $T$ or $T_i$ | 6 |

When $v_\parallel$ is used, it is always `var_vpar=7`.

For __model 600__, the variable numbers are assigned at compile time.

In addition to __model 600__, there are older, legacy, reduced MHD models shown below. __However__, new developments are not necessarily included in these models and, as such, users are suggested to use __model 600__.

| Model | Description | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| model199 | No parallel flows | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T$ | | | | |
| model303/333 | Parallel flows | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T$ | $v_{\vert\vert}$ | | | |
| model305 | ECCD 1 current model | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T$ | $v_{\vert\vert}$ | $J_{ECCD}$ | | |
| model306 | ECCD 2 current model | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T$ | $v_{\vert\vert}$ | $J_{EC1}$ | $J_{EC2}$ | |
| model400 | Parallel flows; separate $T_e$ and $T_i$ | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T_i$ | $v_{\vert\vert}$ | $T_e$ | | |
| model500/501/555 | Parallel flows and neutral/impurity density | $\psi$ | $u$ | $j$ | $\omega$ | $\rho$ | $T$ | $v_{\vert\vert}$ | $\rho_{n/imp}$ | | |


## [Full MHD](full_mhd.md)

| Model | Variables: 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| model710 | $A_3$ | $A_R$ | $A_Z$ | $u_R$ | $u_Z$ | $u_\phi$ | $\rho$ | $T$ | | |
| model711 | $A_3$ | $A_R$ | $A_Z$ | $u_R$ | $u_Z$ | $u_\phi$ | $\rho$ | $Ti$ | $Te$ | |
| model712 | $A_3$ | $A_R$ | $A_Z$ | $u_R$ | $u_Z$ | $u_\phi$ | $\rho$ | $Ti$ | $Te$ | $\rho_n$ |


Refer also to the page about [Notation Conventions](../notation.md).

## Derivations by Emmanuel Franck et al. 

  * [Part I: Document on two-fluid and single-fluid MHD](assets/hierarchymhd.pdf)
  * [Part II: Document on reduced MHD](assets/reduced_mhd.pdf)
