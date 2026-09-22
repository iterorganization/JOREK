---
title: "Set Up Resistivity"
nav_order: 2
parent: "Plasma Functions and Coefficients"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---
# Setting Up the Resistivity

The input parameter `eta` sets the resistivity at the equilibrium plasma centre ($\psi_n=0$), in JOREK units. The resistivity used elsewhere in the plasma can be constant or can follow a Spitzer-like temperature dependence. Variations of the Coulomb logarithm and, in impurity models, $Z_{eff}$ can also be included.

## Simple Spitzer-like setup

A typical input setup is

```fortran
eta                  = 1.d-8 ! Central resistivity in JOREK units
eta_T_dependent      = .true.
eta_coul_log_dep     = .true.
T_max_eta            = 1.d99
```

Here, `eta` is the value at the equilibrium central electron temperature $T_{e0}$. When `eta_T_dependent = .true.`, JOREK evaluates

$$
\eta_T = \eta\left(\frac{T_{e,\mathrm{corr}}}{T_{e0}}\right)^{-3/2}.
$$

For a model with separate temperatures, the evolved $T_e$ is used. For a one-temperature model, JOREK uses $T_e=T/2$; the factor of two cancels in the ratio above.

If `eta_coul_log_dep = .true.`, the result is additionally multiplied by the local Coulomb logarithm normalized to its central value:

$$
\eta_T \longrightarrow \eta_T\frac{\ln\Lambda}{\ln\Lambda_0}.
$$

Thus, `eta` remains the central value even when the Coulomb-logarithm dependence is enabled. This option is currently available in models that call the common resistivity function, including models 600 and 750.

## Choosing `eta`

The physical Spitzer resistivity used for reference is

$$
\eta_{\mathrm{Spitzer}}^{SI}
  = 1.65\times10^{-9}\,\ln\Lambda_0
    \left(T_{e0}[\mathrm{keV}]\right)^{-3/2}\;\Omega\,\mathrm{m}.
$$

you will need to convert it to JOREK units using the [normalization](/JOREK/physics/normalization.md). For convenience, in the logfile you will find

```text
eta_Spitzer (not input parameter; printed for reference in JOREK units)
```

`eta_Spitzer` is a diagnostic reference, not an input parameter: JOREK does **not** automatically assign it to `eta`. To use the nominal Spitzer value, read the printed value from a short initialization run and set `eta` to that value in the input file.

You can deliberately choose a different `eta` to change the resistive timescale, provided that the consequences for the Lundquist number and Ohmic heating are understood.

## Available options

### Constant resistivity

For a resistivity that is constant in both temperature and density, disable both dependencies:

```fortran
eta                  = 1.d-6
eta_T_dependent      = .false.
eta_coul_log_dep     = .false.
```

Disabling only `eta_T_dependent` is not sufficient for a strictly constant profile when `eta_coul_log_dep` remains enabled, because local variations of $\ln\Lambda$ still modify `eta`.

### Temperature-dependent resistivity

Enable the Spitzer-like $T_e^{-3/2}$ scaling with

```fortran
eta_T_dependent = .true.
```

The function uses the positive corrected temperature $T_{e,\mathrm{corr}}$. The correction prevents a non-integer power from being evaluated with a negative temperature; see [Correct Negative Densities / Temperatures](/howto/corr_neg.html).

When the raw temperature falls below `T_min`, the resistivity is frozen at

$$
\eta_T = \eta\left(\frac{T_{min}}{T_{e0}}\right)^{-3/2}.
$$

This avoids an unbounded resistivity in very cold or numerically negative-temperature regions.

### Upper temperature cutoff

`T_max_eta` limits the temperature used in the Spitzer scaling:

```fortran
T_max_eta = 1.d99 ! Effectively no upper cutoff (default)
```

For $T_{e,\mathrm{corr}}>T_{max,\eta}$, JOREK uses

$$
\eta_T = \eta\left(\frac{T_{max,\eta}}{T_{e0}}\right)^{-3/2}.
$$

This imposes a minimum resistivity in the hottest regions. `T_max_eta` is expressed in the same normalized temperature units as the evolved temperature and should be changed only for a specific numerical reason.

### Coulomb-logarithm dependence

Enable local density and temperature effects in $\ln\Lambda$ with

```fortran
eta_coul_log_dep = .true.
```

The common plasma-function implementation calculates $n_e$ in $\mathrm{cm}^{-3}$ and the corrected $T_e$ in eV. It uses

$$
\ln\Lambda =
\begin{cases}
23.0 - \frac{1}{2}\ln n_e + \frac{3}{2}\ln T_e,
  & T_e < 10\;\mathrm{eV}, \\
24.1513 - \frac{1}{2}\ln n_e + \ln T_e,
  & T_e \ge 10\;\mathrm{eV}.
\end{cases}
$$

The density passed to this expression is bounded from below at $10^{10}\;\mathrm{cm}^{-3}$ to prevent nonphysical values of the logarithm.

### Impurity correction

When the selected model evolves impurities, the temperature-dependent resistivity is automatically multiplied by a $Z_{eff}$ correction. The factor implemented in the common plasma function is

$$
f(Z_{eff}) =
\frac{
  Z_{eff}\left(1+1.198Z_{eff}+0.222Z_{eff}^2\right)
  }{
  1+2.966Z_{eff}+0.753Z_{eff}^2
  }
\left[
\frac{1+1.198+0.222}{1+2.966+0.753}
\right]^{-1},
$$

which is normalized so that $f(1)=1$. This factor is applied by the current implementation when `eta_T_dependent = .true.`.

## Input parameters at a glance

| Parameter | Purpose | Default |
|---|---|---:|
| `eta` | Central resistivity in JOREK units | `1.d-5` |
| `eta_T_dependent` | Apply the $T_e^{-3/2}$ dependence | `.true.` |
| `eta_coul_log_dep` | Multiply by $\ln\Lambda/\ln\Lambda_0$ where supported | `.true.` |
| `T_max_eta` | Upper temperature used to evaluate `eta` | `1.d99` |

## Other considerations

* In RMHD there is formally no poloidal resistivity
* Neoclassical effects may increase the resisitivity via trapped particles

Trapped electrons do not carry a current, thus neoclassical reduction leads to an additional radial dependency:
$$
\eta_{||,neo}=\eta_{||}\;\left(1-\sqrt{r/R_0}\right)^{-2}
$$
For $r/R=0.3$ the correction factor is $5$, but this also depends on the collisionality and de-trapping time scale.

