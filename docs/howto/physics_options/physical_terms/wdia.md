---
title: "Include Diamagnetic Drift in Viscosity"
nav_order: 2
parent: "Activation of Physical Terms"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# How to include diamagnetic drift in the viscosity term

- **When is this useful?** When using diamagnetic terms with low resistivity.
- **Implemented for:** Models 303, 307, 333, 500, and 600.
- **How is it activated?** Set `Wdia = .true.` in the namelist input file. The contribution also requires nonzero `tauIC` and `visco` values.


# Viscous diamagnetic term derivation

In the current equations, we have the viscous term as

$$
-\int \mu \nabla v \cdot \nabla W \, dV,
$$

so that the actual physics term, before integration by parts, is

$$
+\int v \, \nabla \cdot \left(\mu \nabla W\right)\, dV.
$$

That is, the viscosity is inside the divergence.

For the viscous diamagnetic term, in weak form, we start with

$$
\int v \, \nabla \cdot \left(\mu \nabla W_{\mathrm{dia}}\right)\, dV.
$$

Integration by parts gives

$$
-\int \mu \nabla v \cdot \nabla W_{\mathrm{dia}}\, dV
+
\int \nabla \cdot
\left(
\mu v \nabla W_{\mathrm{dia}}
\right)\, dV.
$$

Assuming the boundary term vanishes,

$$
-\int \mu \nabla v \cdot \nabla W_{\mathrm{dia}}\, dV.
$$

A second integration by parts gives

$$
+\int
\left[
\nabla \cdot \left(\mu \nabla v\right)
\right]
W_{\mathrm{dia}}\, dV
-
\int
\nabla \cdot
\left(
\mu \nabla v \, W_{\mathrm{dia}}
\right)\, dV.
$$

Again assuming the boundary term vanishes, the final form is

$$
+\int
W_{\mathrm{dia}}
\left[
\nabla \cdot \left(\mu \nabla v\right)
\right]
\, dV.
$$

Expanding the divergence,

$$
+\int
W_{\mathrm{dia}}
\nabla\mu \cdot \nabla v
\, dV
+
\int
\mu W_{\mathrm{dia}}
\nabla^2 v
\, dV.
$$

If $\mu$ depends on temperature,

$$
\nabla \mu
=
\frac{\partial \mu}{\partial T}\nabla T,
$$

and therefore

$$
+\int
W_{\mathrm{dia}}
\frac{\partial\mu}{\partial T}
\nabla T \cdot \nabla v
\, dV
+
\int
\mu W_{\mathrm{dia}}
\nabla^2 v
\, dV.
$$

## JOREK implementation

In JOREK, with

$$
dV = R\,\mathrm{xjac},
$$

the corresponding terms are

$$
+
R
\frac{\partial\mu}{\partial T}
W_{\mathrm{dia}}
\left(
v_x T_x + v_y T_y
\right)
\mathrm{xjac},
$$

and

$$
+
R \mu W_{\mathrm{dia}}
\left(
v_{xx}
+
\frac{v_x}{R}
+
v_{yy}
\right)
\mathrm{xjac}.
$$

Here,

$$
W_{\mathrm{dia}}
=
\frac{\tau}{\rho}
\left(
p_{xx}
+
\frac{p_x}{R}
+
p_{yy}
\right)
-
\frac{\tau}{\rho^2}
\left(
\rho_x p_x
+
\rho_y p_y
\right).
$$