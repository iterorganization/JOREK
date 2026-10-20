---
title: "Run with Ohmic Heating"
nav_order: 4
parent: "Activation of Physical Terms"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

## Models 199/303/333/400/**600**

By default, the Ohmic heating term is switched off in the simulations. **But note that this breaks energy conservation**. When switched on, the following term is added to the energy equation for the electrons
$$
\eta_{ohm}  \, j_{ohm}^2
$$
where $j_{ohm}=j-j_{RE}$ is the ohmic current density and $\eta_{ohm}$ is the plasma resistivity. Note that a factor ($\gamma-1$) must be added in the electron pressure form of the equation.   

The Ohmic heating term can be included by adding the following parameters to the input file:

```ini
# Setting normal resistivity in electric field equation
eta           = 1.d-7
T_max_eta     = 1.d99

# Setting resistivity for the ohmic heating term
eta_ohmic     = 1.d-7
T_max_eta_ohm = 1.d99
```

`eta_ohmic` is the resistivity value, in JOREK units, at the equilibrium plasma core that is used in the Ohmic heating term. 
Any temperature dependence for the resistivity `eta` is automatically replicated for `eta_ohmic`. The parameters used to setup the ohmic heating resistivity are analogous to that of the [$\eta$ term in the equation for the electric field](/JOREK/howto/physics_options/plasma_functions_and_coefficients/spitzer_resistivity.html) (just the suffix "_ohmic" is added). Whenever `eta` differs from `eta_ohmic`, a warning message will appear, since in this case the magnetic energy is not properly transformed into thermal energy.


## Models 710/711/712

Similar to models 199/303/333/400/600, the Ohmic heating term can be included in full MHD models using the `eta_ohmic` parameter in the input file.

The implementation of Ohmic heating terms in full MHD models is different from that in reduced MHD models. For this, the second derivatives of the magnetic vector potential are needed, which requires extending the Fourier spectral method in the $\phi$ direction to include the second derivatives of the basis functions.

The components of the current vector $\boldsymbol{J} = \nabla \times \boldsymbol{B}$ are:

$$
\begin{align*}
J_R ={}& \frac{\partial B_\phi}{\partial Z} - \frac{1}{R}\frac{\partial B_Z}{\partial \phi} \\
J_Z ={}& \frac{1}{R} \left(\frac{\partial B_R}{\partial \phi} - \frac{\partial (R B_\phi)}{\partial R} \right) \\
J_\phi ={}& \frac{\partial B_Z}{\partial R} - \frac{\partial B_R}{\partial Z}
\end{align*}
$$

These can be further rewritten as $\boldsymbol{J} = \nabla \times \nabla \times \boldsymbol{A}$:

$$
\begin{align*}
J_R ={}&
\left(
\frac{\partial^2 A_Z}{\partial Z \partial R}
- \frac{\partial^2 A_R}{\partial Z^2}
+ \frac{1}{R} \frac{\partial F}{\partial \psi} \frac{\partial \psi}{\partial Z}
\right)
- \frac{1}{R^2}
\left(
\frac{\partial^2 A_R}{\partial \phi^2}
- \frac{\partial^2 \psi}{\partial \phi \partial R}
\right)
\\[0.5em]
J_Z ={}&
\frac{1}{R^2}
\left(
\frac{\partial^2 \psi}{\partial \phi \partial Z}
- \frac{\partial^2 A_Z}{\partial \phi^2}
\right)
-
\left(
\frac{\partial^2 A_Z}{\partial R^2}
- \frac{\partial^2 A_R}{\partial R \partial Z}
+ \frac{1}{R} \frac{\partial F}{\partial \psi} \frac{\partial \psi}{\partial R}
- \frac{F}{R^2}
\right)
\\
&-
\frac{1}{R}
\left(
\frac{\partial A_Z}{\partial R}
- \frac{\partial A_R}{\partial Z}
+ \frac{F(\psi)}{R}
\right)
\\[0.5em]
J_\phi ={}&
\frac{1}{R}
\left(
\frac{\partial^2 A_R}{\partial R \partial \phi}
- \frac{\partial^2 \psi}{\partial R^2}
\right)
-
\frac{1}{R^2}
\left(
\frac{\partial A_R}{\partial \phi}
- \frac{\partial \psi}{\partial R}
\right)
-
\frac{1}{R}
\left(
\frac{\partial^2 \psi}{\partial Z^2}
- \frac{\partial^2 A_Z}{\partial Z \partial \phi}
\right)
\end{align*}
$$

For the details of the implementation, please see: [ohmic_heating_in_fmhd.pdf](assets/ohmic_heating/ohmic_heating_in_fmhd.pdf)
