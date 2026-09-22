---
title: "Run with Taylor-Galerkin Stabilization"
nav_order: 1
parent: "Numerics and Stabilization"
grand_parent: "Howto"
layout: default
render_with_liquid: false
---

# Taylor–Galerkin stabilization

Taylor–Galerkin stabilization is available in `model600`. It adds stabilization to selected advective terms, with a separate coefficient for each supported equation. All coefficients default to zero, which disables the corresponding stabilization.

The current implementation does not include the diamagnetic velocity in the Taylor–Galerkin terms.

## Input parameters for `model600`

Unlike older models, `model600` does **not** use the equation-indexed `tgnum` array. Because `model600` is a model family whose evolved variables depend on compile-time extensions, it uses named parameters instead:

| Parameter | Equation | When it is relevant |
|---|---|---|
| `tgnum_u` | Perpendicular momentum | Always available |
| `tgnum_rho` | Main-plasma density | Always available |
| `tgnum_T` | Total temperature | Single-temperature model: `with_TiTe=.false.` |
| `tgnum_Ti` | Ion temperature | Two-temperature model: `with_TiTe=.true.` |
| `tgnum_Te` | Electron temperature | Two-temperature model: `with_TiTe=.true.` |
| `tgnum_vpar` | Parallel momentum | `with_vpar=.true.` |
| `tgnum_rhoimp` | Impurity density | `with_impurities=.true.` |

The general namelist also exposes `tgnum_psi`, `tgnum_zj`, `tgnum_w`, `tgnum_rhon`, `tgnum_nre`, `tgnum_AR`, `tgnum_AZ`, and `tgnum_A3`. The current `model600` element-matrix implementation does not contain Taylor–Galerkin terms controlled by these parameters, so setting them has no effect in `model600`.

## Basic example

For the base, single-temperature `model600` configuration, a possible starting point is:

```fortran
tgnum_u   = 0.5d0
tgnum_rho = 0.5d0
tgnum_T   = 0.5d0
```

For a two-temperature configuration with parallel velocity, use the coefficients associated with the compiled equations:

```fortran
tgnum_u    = 0.5d0
tgnum_rho  = 0.5d0
tgnum_vpar = 0.5d0
tgnum_Ti   = 0.5d0
tgnum_Te   = 0.5d0
```

If the impurity-density extension is compiled, its stabilization can be enabled separately:

```fortran
tgnum_rhoimp = 0.5d0
```

Do not set both `tgnum_T` and `tgnum_Ti`/`tgnum_Te` for the same executable: only the coefficient or coefficients corresponding to the compiled temperature model are used.

## Choosing the coefficients

A value of `0.5` is a reasonable initial value for testing, but it is not universal. The appropriate value depends on the case, time step, and spatial resolution.

Perform a parameter scan that includes the unstabilized case and verify that the selected coefficients do not significantly change resolved physical results, particularly linear growth rates. Set only the coefficients needed for equations showing numerical instability; the coefficients do not need to have identical values.

The active values are printed in the JOREK logfile at startup as `tgnum_u`, `tgnum_rho`, and the other named parameters. A printed nonzero value for a compile-time-disabled or unsupported equation does not make that stabilization term active.
