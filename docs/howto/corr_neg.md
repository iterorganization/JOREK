---
title: "Correct Negative Densities / Temperatures"
nav_order: 2
parent: "Numerics and Stabilization"
grand_parent: "Howto"
layout: default
render_with_liquid: false
---

# Correcting Negative Densities and Temperatures

Negative densities or temperatures can cause floating-point exceptions when they are used in expressions such as $T^{-3/2}$. JOREK therefore provides smooth correction functions that can supply a positive value to calculations which cannot accept a negative input.

> **Warning:** This correction is a numerical workaround, not a physical solution. Negative density or temperature usually indicates a problem such as insufficient resolution, an excessive time step, poorly resolved sources or boundary layers, or a nonlinear solver that has not converged adequately. Always investigate that problem first. The correction prevents some invalid evaluations; it does not make the underlying solution physical.

The functions are implemented in `models/corr_neg.f90` and `models/corr_neg_include.f90`:

- `corr_neg_temp` corrects a temperature supplied to it.
- `corr_neg_dens` corrects a density supplied to it.

For example, code which requires a positive temperature can use

```fortran
corr_neg_temp(T0)**(-1.5d0)
```

The correction returns a separate value. It does **not** clip or overwrite the evolved temperature or density field. Consequently, a run can still contain negative values even though selected coefficients are evaluated using corrected ones.

## Shape of the smooth correction

For an input value $x$, the correction used by both functions is

$$
f(x) =
\begin{cases}
x, & x \geq L_1 + L_2, \\
L_1 + L_2 \exp\left(\dfrac{x-(L_1+L_2)}{L_2}\right), & x < L_1 + L_2,
\end{cases}
$$

where

$$
L_1=c_1x_{\min}, \qquad L_2=c_2x_{\min}.
$$

Here $x_{\min}$ is the temperature or density correction scale, and $c_1$ and $c_2$ are the corresponding correction coefficients. This description replaces the old plot:

- Above $L_1+L_2$, the value is unchanged.
- Below that point, the corrected value bends smoothly away from the original value.
- For increasingly negative input, the corrected value approaches $L_1$.
- The value and its first derivative are continuous at the transition.

The default coefficients are $c_1=c_2=0.5$. With these defaults, the correction starts at $x=x_{\min}$ and approaches $x_{\min}/2$ for very negative input. The first coefficient therefore controls the asymptotic lower value, while the second controls the width of the transition. In normal use, keep the second coefficient positive; changing these coefficients should be done only after checking the resulting correction shape.

The input parameters are named `corr_neg_temp_coef` and `corr_neg_dens_coef`:

```fortran
corr_neg_temp_coef = 0.5d0, 0.5d0
corr_neg_dens_coef = 0.5d0, 0.5d0
```

Within the source code, a call can override the input coefficients:

```fortran
corr_neg_temp(T0, (/ 0.5d0, 0.5d0 /))
```

An overload with a third argument can also supply a local correction scale instead of the global one:

```fortran
corr_neg_temp(T0, (/ 0.5d0, 0.5d0 /), local_T_min)
```

## Choosing the correction scales

Set positive correction scales explicitly in the input file using `T_min_neg` and `rho_min_neg`. For example:

```fortran
T_min_neg   = 4.02d-4 ! 2.01d-5 * central_density * Tmin_eV
rho_min_neg = 1.d-3
```

For `central_density = 1`, the temperature value in this example corresponds to `Tmin_eV = 20`. Both parameters are in JOREK units.

If `T_min_neg` or `rho_min_neg` remains negative, JOREK instead uses the lower equilibrium-profile value, `T_1` or `rho_1`, as the respective correction scale and prints a warning. Explicit positive values are preferable because they express a deliberate, physics-based choice rather than depending on the initialization profiles.

`T_min_neg` and `rho_min_neg` are independent of `T_min` and `rho_min`. The `_min_neg` parameters define the scales of the smooth correction described above. The latter parameters are separate lower limits used in selected calculations, including temperature-dependent coefficients such as parallel heat conductivity and viscosity.

## Enhanced diffusion below a threshold

Another way to remove a localized negative excursion is to increase the heat or particle diffusion in the affected region. These controls are different from `corr_neg_temp` and `corr_neg_dens`: they change transport coefficients in the evolution equations rather than only supplying corrected values to coefficient calculations.

Unlike `implicit_heat_source`, these transport terms do not create energy or particles locally:

- Heat conduction redistributes thermal energy and is globally energy-conserving, apart from heat transported through the domain boundary and the usual discretization error.
- Particle diffusion redistributes particles and conserves the corresponding particle inventory, apart from boundary fluxes. Because density is coupled to the other equations, its effect on the complete energy balance should nevertheless be monitored.

Enhanced transport can therefore be preferable to adding artificial heat. It can still alter the physical solution substantially, especially by increasing losses through an open boundary.

> **Convergence warning:** The activation is a sharp conditional switch, not a smoothly varying coefficient. The threshold crossing itself is not treated implicitly in the Jacobian. During time iterations, a point can repeatedly move across the threshold and alternate between the normal and replacement coefficients. Large coefficient jumps can therefore make solver convergence extremely slow. Increase these coefficients cautiously and check time-step and resolution sensitivity.

The configured value **replaces** the normal local coefficient below its threshold; it is not added to it. It only enhances diffusion if it is larger than the coefficient that would otherwise be used at that location.

### Heat diffusion for a one-temperature model

| Parameter | Default | Function |
|:---|---:|:---|
| `ZK_prof_neg_thresh` | `0.d0` | Activate the perpendicular replacement when $T$ is below this value. |
| `ZK_prof_neg` | `1.d-5` | Replacement perpendicular heat-diffusion coefficient. |
| `ZK_par_neg_thresh` | `0.d0` | Activate the parallel replacement when $T$ is below this value. |
| `ZK_par_neg` | `1.d-3` | Replacement parallel heat-diffusion coefficient. |

For example:

```fortran
ZK_prof_neg_thresh = 0.d0
ZK_prof_neg        = 1.d-5
ZK_par_neg_thresh  = 0.d0
ZK_par_neg          = 1.d-3
```

With zero thresholds, the replacement is selected only where the raw temperature is negative. A small positive threshold activates it before the temperature reaches zero, but also modifies a larger part of the physical low-temperature region.

In model 600, `ZK_par_neg` overrides the normal constant or temperature-dependent parallel conductivity below `ZK_par_neg_thresh`. When `use_zkperp_times_density = .true.`, the perpendicular replacement follows the same density scaling as the ordinary perpendicular coefficient.

### Heat diffusion for a two-temperature model

The two-temperature formulation provides independent controls for ions and electrons:

| Species and direction | Activation threshold | Replacement coefficient | Defaults |
|:---|:---|:---|:---|
| Ion, perpendicular | `ZK_i_prof_neg_thresh` | `ZK_i_prof_neg` | `0.d0`, `1.d-5` |
| Ion, parallel | `ZK_i_par_neg_thresh` | `ZK_i_par_neg` | `0.d0`, `1.d-3` |
| Electron, perpendicular | `ZK_e_prof_neg_thresh` | `ZK_e_prof_neg` | `0.d0`, `1.d-5` |
| Electron, parallel | `ZK_e_par_neg_thresh` | `ZK_e_par_neg` | `0.d0`, `1.d-3` |

Each threshold is compared with its own raw species temperature. For example, a negative electron temperature activates only the electron replacements; it does not automatically change the ion coefficients.

### Particle diffusion for low or negative density

`D_prof_neg` is the analogous replacement for density diffusion. Despite appearing alongside the negative-temperature controls, it is activated by density thresholds, not by temperature.

| Parameter | Default | Model-600 behavior |
|:---|---:|:---|
| `D_prof_neg` | `1.d-5` | Replacement particle-diffusion coefficient. For an activated species, it replaces both the local perpendicular and parallel coefficients. |
| `D_prof_neg_thresh` | `0.d0` | Activate for the background-species density when $\rho-\rho_{\mathrm{imp}}$ is below this threshold. Without impurities, this reduces to the evolved density. |
| `D_prof_imp_neg_thresh` | `-1.d3` | Activate `D_prof_neg` for the impurity density below this threshold. The very negative default disables this safeguard in ordinary runs because it can cause convergence problems near zero impurity density. |
| `D_prof_tot_neg_thresh` | `0.d0` | Activate `D_prof_neg` for both background and impurity transport when the total evolved density is below this threshold, provided the background-species condition has not already activated. |
| `D_imp_extra_neg` | `1.d-6` | Replacement diffusion coefficient for the additional impurity/neutral-density transport channels available in the relevant extended models. |
| `D_imp_extra_neg_thresh` | `-1.d3` | Activation threshold for `D_imp_extra_neg`; its very negative default normally disables it. |

The exact density represented by the impurity variable differs between extended models, so check the selected model before using the impurity thresholds. In particular, activating a correction as soon as an impurity density crosses zero can cause the coefficient to switch repeatedly while an initially absent impurity species oscillates around zero.

### Recommended use

Start with the default zero temperature or background-density thresholds, and choose replacement coefficients by comparing them with the normal local transport coefficients. Then:

1. Check whether the negative region shrinks through redistribution rather than simply moving to a boundary.
2. Monitor boundary heat and particle fluxes; conservative diffusion can still increase energy or particle loss from the computational domain.
3. Monitor nonlinear iteration counts. If convergence deteriorates, reduce the jump between the ordinary and replacement coefficients or reduce the time step.
4. Repeat with smaller replacement coefficients and thresholds closer to zero to establish that the physical conclusions do not depend on this safeguard.

## Implicit heat source in model 600

Model 600 also provides `implicit_heat_source`. This is different from merely evaluating a coefficient with `corr_neg_temp`: it adds a positive term directly to the temperature equation wherever the temperature is below the selected threshold.

> **Severe warning — this source can hide real numerical or physical problems and can completely break energy conservation.** It injects thermal energy without a compensating loss elsewhere. A run may appear stable only because this artificial source masks negative-temperature excursions. Do not enable it as a substitute for resolving time-step, mesh, transport, source, boundary-condition, or solver-convergence problems. Results involving energy balances are not trustworthy unless the injected energy has been explicitly measured and shown to be negligible.

The default setting disables the source:

```fortran
implicit_heat_source = 0.d0
```

Setting it to `1.d0` fully enables the implemented stabilization:

```fortran
implicit_heat_source = 1.d0 ! Emergency numerical stabilization only
```

For the one-temperature formulation, define $T_*=\texttt{T_min_neg}$. Apart from the weak-form and normalization factors, the added source has the shape

$$
S_{\mathrm{impl}}(T) = a(\gamma-1)\left[
\frac{T_*}{2}\left(1+
\exp\left(\frac{\min(T,T_*)-T_*}{T_*/2}\right)\right)
-\min(T,T_*)
\right],
$$

where $a=\texttt{implicit_heat_source}$. The term is zero for $T>T_* $. It is included consistently in the nonlinear residual and its Jacobian, which is why it is called *implicit*; that does not make it a physical or energy-conserving source.

In the two-temperature formulation, the same expression is applied separately to $T_i$ and $T_e$, using $T_*=\texttt{T_min_neg}/2$ for each species.

Before considering this source, first try to remove the negative-temperature excursions by reducing the time step, checking nonlinear convergence, improving spatial resolution, and reviewing physical sources, boundary conditions, and transport or stabilization coefficients. If the source is temporarily required for diagnosis or recovery:

1. Record the value used and why it was needed.
2. Inspect the model-600 RHS terms `T_Eq__impl_heating`, `Ti_Eq__implicit_heating`, and `Te_Eq__implicit_heating`; see [Plotting Separate Equation Terms in VTK](../diagnostics/plot_rhs_terms.html).
3. Quantify the artificial energy input and report it with any energy balance.
4. Repeat the run while reducing `implicit_heat_source` toward zero, and do not rely on conclusions that disappear when the source is removed.
