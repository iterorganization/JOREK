---
title: "Variable MultiScale (VMS) Stabilization"
nav_order: 6
parent: "Numerics and Stabilization"
grand_parent: "Howto"
layout: default
render_with_liquid: false
---

# Variable MultiScale (VMS) Stabilization

VMS stabilization is currently implemented for full MHD model 750.

## Usage

Set `use_vms = .true.` in the JOREK input file and assign suitable values to the coefficients for the evolved variables:

```fortran
use_vms = .true.

vms_coeff_AR     = 1.d-3
vms_coeff_AZ     = 1.d-3
vms_coeff_A3     = 1.d-3
vms_coeff_UR     = 1.d-3
vms_coeff_UZ     = 1.d-3
vms_coeff_Up     = 1.d-3
vms_coeff_rho    = 1.d-3
vms_coeff_T      = 1.d-3
vms_coeff_Ti     = 1.d-3
vms_coeff_Te     = 1.d-3
vms_coeff_rhon   = 1.d-3
vms_coeff_rhoimp = 1.d-3
```

The preset value is `use_vms = .false.`, so VMS stabilization is disabled by default. For `n_order = 3`, the recommended starting value is `0.001` for each active variable.

## Derivation and Implementation

Let $\mathcal{R}(\boldsymbol{w})$ be the residual

$$
\begin{aligned}
\mathcal{R}(\boldsymbol{w})
&:=\mathbb{M}(\boldsymbol{w})\partial_t\boldsymbol{w}
  +\mathcal{L}(\boldsymbol{w}) \\
&=\mathbb{M}(\boldsymbol{w})\partial_t\boldsymbol{w}
  +\mathbb{C}\nabla\boldsymbol{w}
  -\nabla\cdot(\mathbb{D}\nabla\boldsymbol{w})
  -\mathbb{S}(\boldsymbol{w})
=0.
\end{aligned}
$$

Here, $\boldsymbol{w}\subset\boldsymbol{W}$ forms a well-defined system of MHD equations, with

$$
\boldsymbol{W}
=\left(
A_R,\,A_Z,\,\psi,\,u_R,\,u_Z,\,u_\phi,\,
\rho,\,T,\,T_i,\,T_e,\,\rho_n,\,\rho_f
\right)^\mathrm{T}.
$$

The Galerkin finite-element method (FEM) for this system of partial differential equations is

$$
\int \boldsymbol{w}^*\cdot\mathcal{R}(\boldsymbol{w})\,d\Omega=0.
$$

The VMS-stabilized FEM formulation is

$$
\int \boldsymbol{w}^*\cdot\mathcal{R}(\boldsymbol{w})\,d\Omega
-\int
\mathcal{L}^\mathrm{T}\boldsymbol{w}^*
\cdot\mathbb{T}\mathcal{R}(\boldsymbol{w})\,d\Omega'
=0.
$$

The first term is the finite-element weak form and the second is the VMS stabilization.

> **Implementation remark:** In the implicit method, $\mathcal{L}^\mathrm{T}\boldsymbol{w}^*$ is assumed to be evaluated at the previous time step. Its Jacobian is therefore not taken into account.

## Residual and Practical Difficulties

The residual appearing in the VMS stabilization terms is approximated as

$$
\mathcal{R}
\approx\widetilde{\mathcal{R}}
=\mathbb{C}\nabla\boldsymbol{w}
-\nabla\cdot(\mathbb{D}\nabla\boldsymbol{w})
-\mathbb{S}(\boldsymbol{w}).
$$

The following simplifications are made:

- Transient (time-derivative) terms are ignored.
- The $\boldsymbol{J}\times\boldsymbol{B}$ term is ignored for simplicity. Including it would require extending the Fourier method.
- Diffusivities parallel to the magnetic field are ignored for simplicity.

## Stabilization Operator

Consider the two-temperature full MHD model with impurities. Its partial differential-equation operator is

$$
\mathcal{L}=
\begin{bmatrix}
-(\boldsymbol{v}\times\nabla\times)+\nabla^2
  & \boldsymbol{0} & 0 & 0 & 0 & 0 \\
\boldsymbol{0}
  & \rho\boldsymbol{v}\cdot\nabla-\nabla^2
  & (T_i+T_e)\nabla & \rho\nabla & \rho\nabla & 0 \\
\boldsymbol{0}
  & \rho\nabla\cdot
  & \boldsymbol{v}\cdot\nabla-\nabla^2 & 0 & 0 & 0 \\
\boldsymbol{0}
  & \gamma\rho T_i\nabla\cdot
  & T_i\boldsymbol{v}\cdot\nabla
  & \rho_{\mathrm{eff}}^i\boldsymbol{v}\cdot\nabla
    -\nabla\cdot(\boldsymbol{\kappa}_i\nabla)
  & 0
  & T_{\mathrm{eff}}^i\boldsymbol{v}\cdot\nabla \\
\boldsymbol{0}
  & \gamma\rho T_e\nabla\cdot
  & 0
  & T_e\boldsymbol{v}\cdot\nabla
  & \rho_{\mathrm{eff}}^i\boldsymbol{v}\cdot\nabla
    -\nabla\cdot(\boldsymbol{\kappa}_e\nabla)
  & T_{\mathrm{eff}}^e\boldsymbol{v}\cdot\nabla \\
\boldsymbol{0}
  & \rho_f\nabla\cdot
  & 0 & 0 & 0
  & \boldsymbol{v}\cdot\nabla-\nabla\cdot(\mathbb{D}_f\nabla)
\end{bmatrix}.
$$

## Stabilization Matrix

The simplest choice is a diagonal stabilization matrix:

$$
\mathbb{T}
=\frac{h_e}{\max(\lambda_e)}
\begin{bmatrix}
C_{A_R} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & C_{A_Z} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & C_{A_3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & C_{v_R} & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & C_{v_Z} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & C_{v_\phi} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & C_{\rho} & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & C_{T_i} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & C_{T_e} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & C_{\rho_f}
\end{bmatrix}.
$$

Here, $h_e$ is the characteristic element size and $\lambda_e$ is the Alfvén-wave speed evaluated element-wise. The real parameters $C_{A_R}$ and so on correspond to input parameters such as `vms_coeff_AR`.

In principle, the stabilization matrix should be

$$
\mathbb{T}
=\mathcal{L}^{-1}
=\left(
\mathbb{M}\partial_t
+\mathbb{A}\nabla
+\nabla\cdot(\mathbb{D}\nabla)
\right)^{-1}.
$$

However, evaluating this expression is harder than solving the original problem.

## Additional Remarks

- Different formulations and approximations of the residual produce a hierarchy of VMS stabilization methods.
- Different choices of stabilization operator also produce a hierarchy of methods. The simplest choice, for example, is to use the diagonal part of $\mathcal{L}$ without diffusive terms.
- For details and demonstrations, see [Computers & Mathematics with Applications (2023)](https://doi.org/10.1016/j.camwa.2023.04.034).
