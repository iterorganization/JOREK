---
title: "Parallel Electric Field in RMHD"
nav_order: 16
parent: "Expressions and Derivations"
grand_parent: "Reduced MHD models"
layout: default
render_with_liquid: false
---

# Electric field in JOREK

In an actual discharge, the full shape and size of the plasma does not appear at once.

The toroidal field is first switched on, there is gas, and then the "transformer voltage" is applied. This causes a small "spot" of plasma near the inboard side (inboard because of the minimum toroidal path there, causing the maximum volts/meter). This "spot" is then grown quickly to its full size together with the plasma current.

It is done this way because during this growth, since the plasma is still not too hot, there is no problem of E-field/B-field diffusion into the plasma. However, once we reach the hot-plasma state with full $I_p$, the magnetic diffusion timescale is too long for any external electromagnetic changes to cause quick changes in the plasma core.

Therefore, on the timescales of interest, once we reach the hot-plasma state with full $I_p$, the $\partial\mathbf{A}/\partial t$ of the transformer does not enter into the plasma core (although the edge region may see it).

What we have in the full-$I_p$, hot-plasma state is therefore a plasma current $J_0$ with an associated curl-free toroidal electric field $\mathbf{E}_0$ (curl-free because of steady state). One can easily show that the curl-free condition leads to an equilibrium toroidal electric field of the form

$$
E_0 = \frac{C_0}{R},
$$

where $C_0$ is a constant.

The main point of the above is that the equilibrium toroidal electric field $\mathbf{E}_0$ in the plasma cannot be associated or equated with a local $\partial\mathbf{A}/\partial t$ within the plasma space.

Also, in JOREK we do not track/include the time derivative of the boundary flux caused by transformer action. Then, on what premise does the "current source" term enter Faraday's law? This can be seen rather straightforwardly if we uncurl Faraday's law to accommodate the initial no-flow 2D equilibrium state of the plasma.

Uncurling Faraday's law, in the strict interior of the plasma, we can write

$$
\frac{\partial \mathbf{A}}{\partial t}
=
-\mathbf{E}
-
F_0\nabla u
-
\nabla v_l,
$$

without any loss of generality.

The extra term $\nabla v_l$ is necessary because JOREK's usage and interpretation of the electric potential $F_0u$ is incomplete in the sense that it does not include/accommodate an initial electric-potential gradient in the toroidal direction.

The above equation should satisfy the initial condition in our simulations, which is a steady 2D equilibrium with no flows, i.e.

$$
\frac{\partial}{\partial t}(\ldots)=0,
\qquad
\nabla u = 0.
$$

This gives

$$
-\nabla v_l
=
\mathbf{E}_{t=0}
=
E_0\hat{\mathbf e}_{\phi}
=
\eta_0\left(J_0-J_{b0}\right)\hat{\mathbf e}_{\phi}.
$$

The last equality arises by equating the electric force and collisional drag on the electrons, because of negligible electron inertia.

Here,

$$
\eta_0 = \eta(T_0,Z_{\mathrm{eff},0})
$$

is the resistivity computed at the equilibrium/initial state (typically $Z_{\mathrm{eff},0}=1$ in our simulations), $J_0$ is the total equilibrium plasma current density, and $J_{b0}$ is the equilibrium bootstrap current density.

Since the bootstrap current is "transport induced", we cannot subject it to resistivity as used in JOREK. Otherwise, we should also include the collisional force between the passing and trapped particles that causes the bootstrap current in the first place. We therefore do not include both the collisional cause and collisional drag, but instead use an ad-hoc model to simply estimate the bootstrap current.

Substituting the above expression for $\nabla v_l$ into the uncurled Faraday's law gives

$$
\frac{\partial\mathbf{A}}{\partial t}
=
-\mathbf{E}
-
F_0\nabla u
+
\eta_0\left(J_0-J_{b0}\right)\hat{\mathbf e}_{\phi}.
\tag{18.1}
$$

Reaching the above equation in this way avoids the need for introducing the loop-voltage term as an out-of-the-blue term with ad-hoc reasons such as "artificial term to maintain a constant current in our simulations because that is the reality", etc.

Therefore, the last term in the above equation should instead be interpreted simply as the initial ($t=0$) equilibrium electric field in the plasma domain in our simulations.

In the general case, the electric field can then be trivially computed via the necessary projection(s) of

$$
\boxed{
\mathbf{E}
=
-\frac{\partial\mathbf{A}}{\partial t}
-
F_0\nabla u
+
\eta_0\left(J_0-J_{b0}\right)\hat{\mathbf e}_{\phi}
}
\tag{18.2}
$$

Thus, the background loop-voltage term is always included when we compute the electric field $\mathbf{E}$.


## Expression for the parallel electric field

The equation, in normalized form, for the poloidal flux $\psi$ implemented in JOREK is

$$
\frac{1}{R^2}\frac{\partial\psi}{\partial t}
=
\frac{\eta}{R^2}
\left[
j-j_0-(j_b-j_{b0})
\right]
-
\frac{1}{R}[u,\psi]
-
\frac{F_0}{R^2}\frac{\partial u}{\partial\phi}
+
\frac{F_0^2}{R^2B^2}
\frac{\tau_{\mathrm{IC}}}{\rho}
\left[
\frac{1}{R}[p,\psi]
+
\frac{F_0}{R^2}\frac{\partial p}{\partial\phi}
\right].
\tag{18.3}
$$

Here the hyperresistive term, which is not important for the present discussion, has been excluded.

However, as seen from the discussion above, the physically consistent form that should be implemented is

$$
\frac{1}{R^2}\frac{\partial\psi}{\partial t}
=
\frac{\eta}{R^2}
\left(j-j_b\right)
-
\frac{\eta_0}{R^2}
\left(j_0-j_{b0}\right)
-
\frac{1}{R}[u,\psi]
-
\frac{F_0}{R^2}\frac{\partial u}{\partial\phi}
+
\frac{F_0^2}{R^2B^2}
\frac{\tau_{\mathrm{IC}}}{\rho}
\left[
\frac{1}{R}[p,\psi]
+
\frac{F_0}{R^2}\frac{\partial p}{\partial\phi}
\right].
\tag{18.4}
$$

Uncurling Faraday's law gave us the following expression for the electric field:

$$
\mathbf{E}
=
-\frac{1}{R}
\frac{\partial\psi}{\partial t}
\hat{\mathbf e}_{\phi}
-
F_0\nabla u
+
\eta_0
\left(J_0-J_{b0}\right)
\hat{\mathbf e}_{\phi}.
$$

A parallel projection, $\hat{\mathbf b}\cdot(\ldots)$, then easily gives, in normalized form,

$$
E_{\parallel}
=
-\frac{F_0}{B}
\left[
\frac{1}{R^2}\frac{\partial\psi}{\partial t}
+
\frac{1}{R}[u,\psi]
+
\frac{F_0}{R^2}\frac{\partial u}{\partial\phi}
+
\frac{\eta_0}{R^2}(j_0-j_{b0})
\right].
\tag{18.5}
$$

The terms within the brackets on the right-hand side of Eq. (18.5) can be rewritten in an alternative form using Eq. (18.4). This leads to the following expression for the parallel electric field, in normalized form:

$$
E_{\parallel}
=
-\frac{F_0}{B}
\left[
\frac{\eta}{R^2}(j-j_b)
+
\frac{F_0^2}{R^2B^2}
\frac{\tau_{\mathrm{IC}}}{\rho}
\left(
\frac{1}{R}[p,\psi]
+
\frac{F_0}{R^2}\frac{\partial p}{\partial\phi}
\right)
\right].
\tag{18.6}
$$

Alternatively, one can simply use the generalized Ohm's law to obtain an expression for the parallel electric field:

$$
\mathbf{E}
=
-\mathbf{v}\times\mathbf{B}
+
\eta(\mathbf{J}-\mathbf{J}_b)
+
\frac{\mathbf{J}\times\mathbf{B}-\nabla p_e}{ne}.
$$

Therefore,

$$
E_{\parallel}
=
\frac{\mathbf{E}\cdot\mathbf{B}}{B}
=
\frac{\eta(\mathbf{J}-\mathbf{J}_b)\cdot\mathbf{B}}{B}
-
\frac{\mathbf{B}\cdot\nabla p_e}{Bne}.
$$

Now,

$$
\mathbf{B}\cdot\nabla p_e
=
\frac{1}{R}[p_e,\psi]
+
\frac{F_0}{R^2}\frac{\partial p_e}{\partial\phi},
$$

and it is a safe approximation to take

$$
(\mathbf{J}-\mathbf{J}_b)\cdot\mathbf{B}
\approx
(J_\phi-J_{\phi,b})B_\phi
\approx
(J_\phi-J_{\phi,b})\frac{F_0}{R}.
$$

This leads to

$$
E_{\parallel}
=
\eta(J_\phi-J_{\phi,b})\frac{F_0}{BR}
-
\frac{1}{enB}
\left(
\frac{1}{R}[p_e,\psi]
+
\frac{F_0}{R^2}\frac{\partial p_e}{\partial\phi}
\right).
$$

Normalizing the above gives

$$
E_{\parallel}
=
-\frac{F_0}{B}
\left[
\frac{\eta}{R^2}(j-j_b)
+
\frac{\tau_{\mathrm{IC}}}{\rho}
\left(
\frac{1}{R}[p,\psi]
+
\frac{F_0}{R^2}\frac{\partial p}{\partial\phi}
\right)
\right].
\tag{18.7}
$$

Therefore, either Eq. (18.5), Eq. (18.6), or Eq. (18.7) can be used to evaluate $E_{\parallel}$ in JOREK.

The only difference between Eqs. (18.6) and (18.7) is the rather not-so-significant factor

$$
\frac{F_0^2}{R^2B^2}.
$$