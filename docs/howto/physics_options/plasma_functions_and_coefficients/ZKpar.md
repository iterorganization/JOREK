---
title: "Set Up Parallel Heat Conduction (ZKpar)"
nav_order: 3
parent: "Plasma Functions and Coefficients"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# Setting Up Parallel Heat Conduction

The `ZKpar` parameters control the parallel heat-conduction coefficient in JOREK units. You can use a constant coefficient or the Spitzer–Härm-like temperature dependence $T^{5/2}$. In models with separate ion and electron temperatures, the two coefficients are configured independently.

This page describes the common `conductivity_parallel` implementation used by models 600 and 750.

## One-temperature setup

For a model that evolves a single MHD temperature $T=T_i+T_e$, use `ZK_par`:

```fortran
ZK_par             = 1.d0  ! Central parallel heat-conduction coefficient
ZKpar_T_dependent  = .true.
ZK_par_max         = 1.d20
```

Here, `ZK_par` is the coefficient at the equilibrium central temperature `T_0`. When `ZKpar_T_dependent = .true.`, JOREK evaluates

$$
\kappa_{\parallel}(T)
  = \kappa_{\parallel,0}
    \left(\frac{T_{\mathrm{corr}}}{T_0}\right)^{5/2}.
$$

## Separate ion and electron temperatures

When $T_i$ and $T_e$ are evolved separately, use `ZK_i_par` and `ZK_e_par`:

```fortran
ZK_i_par            = 1.d0 ! Central ion coefficient
ZK_e_par            = 1.d0 ! Central electron coefficient
ZKpar_T_dependent   = .true.
ZK_par_max          = 1.d20
```

The common routine is applied independently to each species:

$$
\begin{aligned}
\kappa_{i,\parallel}(T_i)
  &= \kappa_{i,\parallel,0}
     \left(\frac{T_{i,\mathrm{corr}}}{T_{i0}}\right)^{5/2}, \\
\kappa_{e,\parallel}(T_e)
  &= \kappa_{e,\parallel,0}
     \left(\frac{T_{e,\mathrm{corr}}}{T_{e0}}\right)^{5/2}.
\end{aligned}
$$

`ZK_par_max` is shared by the ion and electron coefficients.

## Choosing the central coefficient

JOREK calculates nominal Spitzer–Härm coefficients from the central plasma parameters and prints them in the log for reference. For separate-temperature models, look for

```text
ZK_e_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)
ZK_i_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)
```

For a one-temperature model, look for

```text
ZK_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)
```

These are diagnostic reference values, not input parameters. JOREK does **not** automatically assign them to `ZK_par`, `ZK_i_par`, or `ZK_e_par`. To use the nominal Spitzer–Härm value, read it from a short initialization run and copy it to the corresponding input parameter.

## Available options

### Constant parallel conduction

Disable the temperature scaling with

```fortran
ZKpar_T_dependent = .false.
```

The coefficient is then equal to `ZK_par`, or to `ZK_i_par` and `ZK_e_par` in a separate-temperature model. The negative-temperature override described below is evaluated independently and can still replace this value.

### Maximum coefficient

When the $T^{5/2}$ scaling produces a value larger than `ZK_par_max`, the coefficient is clipped:

$$
\kappa_{\parallel}=\kappa_{\parallel,max}.
$$

Configure the limit with, for example,

```fortran
ZK_par_max = 1.d20 ! Default
```

This limit is mainly a numerical safeguard against excessively large parallel conduction in hot regions. It is applied to the single-temperature, ion, and electron coefficients.

### Minimum temperature used in the scaling

The $T^{5/2}$ expression uses the positive corrected temperature. See [Correct Negative Densities / Temperatures](/howto/corr_neg.html) for the temperature correction.

When `ZKpar_T_dependent = .true.`, the decision to freeze the coefficient is based on the raw, uncorrected temperature. If $T_{raw}<T_{min,ZK\parallel}$, JOREK uses

$$
\kappa_{\parallel}
  = \kappa_{\parallel,0}
    \left(\frac{T_{min,ZK\parallel}}{T_0}\right)^{5/2}.
$$

The corresponding inputs are

```fortran
T_min_ZKpar  = 1.d-6 ! Single-temperature model
Ti_min_ZKpar = 1.d-6 ! Ion temperature
Te_min_ZKpar = 1.d-6 ! Electron temperature
```

All values are in JOREK temperature units. Their preset sentinel value is `-1.d12`; during initialization, any value below `-1.d10` is replaced by the general `T_min` setting. Therefore, if these parameters are not specified, their effective value is `T_min`.


## Input parameters at a glance

| Parameter | Purpose | Preset value |
|---|---|---:|
| `ZK_par` | Central coefficient for a single temperature | `1.d0` |
| `ZK_i_par` | Central ion coefficient | `1.d0` |
| `ZK_e_par` | Central electron coefficient | `1.d0` |
| `ZKpar_T_dependent` | Apply the $T^{5/2}$ dependence | `.true.` |
| `ZK_par_max` | Maximum temperature-dependent coefficient | `1.d20` |
| `T_min_ZKpar` | Minimum single temperature used in the scaling | `-1.d12`, then replaced by `T_min` |
| `Ti_min_ZKpar` | Minimum ion temperature used in the scaling | `-1.d12`, then replaced by `T_min` |
| `Te_min_ZKpar` | Minimum electron temperature used in the scaling | `-1.d12`, then replaced by `T_min` |
