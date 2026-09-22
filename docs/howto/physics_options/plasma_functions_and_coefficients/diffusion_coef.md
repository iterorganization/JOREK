---
title: "Get D_perp and ZK_perp for Stationary Profiles"
nav_order: 1
parent: "Plasma Functions and Coefficients"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# Calculator for Diffusivity and Conductivity Coefficients

Many JOREK users struggle to maintain stationary H-mode ne and T profiles over long simulations (n=0). Achieving this requires selecting transport coefficients consistent with:

- the imposed particle source and heating sources, and
- the desired target or equilibrium profiles (assumed flux functions).

We explain below how to use some tools for doing this.

## Step 1: Get Metric, Source, and Target Profiles from the Equilibrium

For that, we will simply use [`jorek2_postproc`](../../introduction_to_jorek_diagnostics.md#jorek2_postproc), with the command `pol_contour_integrals`. This is what you will have to do inside postproc:

```text
namelist input
set surfaces 300
for step 0 do
  pol_contour_integrals
done
```

This will export a file with all necessary profiles (metric, target, sources) to obtain your diffusion coefficients:

```text
postproc/pol_contour_integrals_s000000.dat
```

## Step 2: Get the Coefficients with Python

Run the Python script in `jorek_repo/util/obtain_diff_coefficients.py`:

```bash
python jorek_repo/util/obtain_diff_coefficients.py postproc/pol_contour_integrals_s000000.dat --skiprows 2
```

This will generate your ready-to-use coefficients in the files:

```text
D_profile.dat
kappa_i.dat
kappa_e.dat
kappa_t.dat   # Total temperature
```

Add those to your input file in `d_perp_file`, `zk_i_perp_file`, `zk_e_perp_file`, and `zk_perp_file`. Then you are ready to go.

## Solved Equation

You can find the employed method here. It is basically a 1D diffusion + sources equation in equilibrium, which has an analytical solution.

**J. Artola**  
*October 29, 2021*

Let's take the equilibrium continuity equation

$$
\nabla \cdot \left(D\nabla \rho\right) = -S,
\tag{1}
$$

where $S$ includes all sinks/sources, including flows.

By transforming the equation into flux-aligned coordinates and doing some algebra, we find an analytical solution for $D$:

$$
D(\psi_N)
=
-
\frac{(\psi_{\mathrm{bnd}}-\psi_{\mathrm{axis}})^2}
{\alpha(\psi_N)\,
\displaystyle\frac{\partial \rho}{\partial \psi_N}}
\int_0^{\psi_N}
\gamma(\psi_N')
\left\langle S \right\rangle_{\psi_N'}
\,d\psi_N'.
\tag{2}
$$

Here, $\left\langle S \right\rangle_{\psi_N'}$ is the flux-surface-averaged net source. One can also derive this equation from the divergence theorem of calculus.

We have defined the following geometrical coefficients:

$$
\alpha(\psi_N)
=
\oint_{\psi_N}
R |\nabla\psi|\,dl_{\mathrm{pol}},
\tag{3}
$$

and

$$
\gamma(\psi_N)
=
\oint_{\psi_N}
\frac{R}{|\nabla\psi|}
\,dl_{\mathrm{pol}},
\tag{4}
$$

where $dl_{\mathrm{pol}}$ is the differential line element of the poloidal contour of the flux surface.

## Remarks

- Equation (2) does not work in the SOL. In that region, the diffusion has to be manually set in order to obtain a reasonable $\rho$ at the separatrix.

- When the SOL is included, Eq. (2) does not completely constrain $\rho$, since Eq. (1) is free up to a global constant in $\rho$. SOL sources and sinks, or boundary conditions, must be tuned for that purpose.

- Equation (2) diverges at $\psi_N = 0$ or when

  $$
  \rho' = 0.
  $$

  The profile must be smoothed in those situations. Also, a source must be given; otherwise a unique solution does not exist.

- Flow terms can be included in $\left\langle S \right\rangle$ by averaging them with `jorek2_postproc`.

- The same formula can be applied to the temperature equation.

## Caveats and Hardcoded Parameters in the Python Script

The Python script is intentionally configured via hardcoded constants near the top of the file (rather than many command-line flags). Before using the tool on a new case, you should review these parameters.

The most important settings are the MIN/MAX bounds for numerical stability.

**SOLUTIONS of the equation do not exist if sources are active but gradients are zero**. For example, flat profiles may exist only in the absence of sources. **Non-monotonic target profiles have no solutions in a diffusion(only)+source problem.**

### 1. Numerical Tolerances (`EPS_*`)

These control how the script detects regions where the analytic formula is ill-conditioned.

- `EPS_INT`: Threshold for considering the cumulative integral of the source term "approximately zero".
- `EPS_DER_RHO`, `EPS_DER_Ti`, `EPS_DER_Te`, `EPS_DER_T`: Thresholds for considering derivatives "approximately zero".

If the `EPS_DER_*` thresholds are too small, numerical noise can produce extremely large coefficients.
If they are too large, you may incorrectly classify broad regions as "flat/unconstrained".

### 2. MIN/MAX Bounds (Critical for Numerical Stability)

After computing and cleaning the profiles, the script clips the coefficients to hardcoded ranges:

- Density diffusivity bounds: `D_MIN_RHO`, `D_MAX_RHO`.
- Heat conductivity bounds: `KAPPA_I_MIN`, `KAPPA_I_MAX`, `KAPPA_E_MIN`, `KAPPA_E_MAX`, `KAPPA_T_MIN`, `KAPPA_T_MAX`.

Use these limits to avoid numerical problems (too small or too large diffusion).

### 3. Core Handling (`EXTR_CORE_PSIN`)

Near the magnetic axis, profiles/metrics can be poorly defined or noisy. The script enforces a clean core behavior on the FINAL computed coefficient profiles:

- It fits a straight line using the first two points with `psi_N >= EXTR_CORE_PSIN`.
- It overwrites all coefficient values for `psi_N < EXTR_CORE_PSIN` using that line.

Adjust `EXTR_CORE_PSIN` if needed:

- Increase it if the axis region is unreliable/noisy.
- Decrease it if you trust the core values.

### 4. SOL Extension (`EXTEND_TO_PSIN` and Ramp Parameters)

The script extends each coefficient profile up to `EXTEND_TO_PSIN` by:

- Adding `N_SOL_POINTS` points beyond the last input point.
- Using a constant tail equal to the last computed value (default).
- Optionally applying a smooth arctan ramp in the SOL.

Controls:

- `EXTEND_TO_PSIN`, `N_SOL_POINTS`.
- `RAMP_AMP` (0 means no ramp), `RAMP_CENTER`, `RAMP_WIDTH`.

This is useful if you want enhanced transport in the SOL while leaving the core/edge (psi_N <= 1) unchanged.

## Practical Troubleshooting

If you see spikes in the coefficients:

- Tighten the relevant MAX values.
- Consider increasing `EPS_DER_*` (avoid dividing by very small derivatives).

If profiles still evolve in time:

- Check that sources and target gradients are physically compatible.
- Consider using `--hard-fail-impossible` to identify incompatible regions early.
