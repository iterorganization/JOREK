---
title: "Run with Neoclassical Effects"
nav_order: 3
parent: "Activation of Physical Terms"
grand_parent: "Physics Options"
layout: default
render_with_liquid: false
---

# How to run JOREK with neoclassical coefficient profiles

## 1. Calculate neoclassical coefficient profiles

### 1.1. Obtain JOREK data for calculating the coefficients

Running the JOREK equilibrium with `neo=.false.` produces the following files:

- `qprofile.dat`, calculated in the `qprofile` subroutine, contains the normalized poloidal flux $\psi_{N,\mathrm{pol}}$ and the safety factor $q(\psi_N)$.
- `T_rho_profiles.dat`, calculated in the `equilibrium` subroutine, contains $T(\psi_N)$ and $\rho(\psi_N)$ in JOREK units.

The neoclassical-coefficient calculation requires:

- the normalized poloidal flux $\psi_N$;
- the safety factor $q(\psi_N)$;
- the electron temperature $T_e(\psi_N)$;
- the electron density $N_e(\psi_N)$.

In JOREK, the poloidal flux and corresponding safety-factor profile are calculated by the `qprofile` subroutine, which is called from `equilibrium` and writes `qprofile.dat`. The $\psi$-dependent temperature and density profiles must have the same number of points as $q(\psi_N)$. The `equilibrium` subroutine writes these profiles to the ASCII file `T_rho_profiles.dat`.

For cleaner profiles, remove all points outside the separatrix from both `qprofile.dat` and `T_rho_profiles.dat`.

### 1.2. Run `neoclassical_program`

The program and its dedicated Makefile are located in `util/neoclass_program/`. This utility is not built by JOREK's top-level Makefile.

#### Prepare the input file

The namelist input for `neoclassical_program` must specify:

- the central density in units of $10^{20}\,\mathrm{m}^{-3}$;
- the magnetic field $B_t$ in tesla;
- the major radius $R_0$ and minor radius $a_\mathrm{min}$;
- the files containing the safety-factor and temperature-density profiles.

For example, for ITER-like parameters:

```fortran
&in1
  central_density  = 1.0                 ! Units: 10^20 m^-3
  Bt               = 5.3d0               ! Magnetic field at the axis
  R0               = 6.19476d0           ! Major radius at the axis
  amin             = 2.d0                ! Minor radius of the equilibrium
  qprofile_file    = 'qprofile_ITER2.dat'
  T_Ne_profile_file = 'Te_ne_profiles_ITER.dat'
/
```

The utility directory contains the example files `qprofile_ITER2.dat` and `Te_ne_profiles_ITER.dat`. For a new JOREK case, replace these names with the paths to the corresponding `qprofile.dat` and `T_rho_profiles.dat` files, or copy those files into the utility's working directory. Both profile files must contain the same number of lines.

#### Compile `neoclassical_program`

The program requires BLAS. Its Makefile currently defaults to the Intel Fortran compiler and Intel MKL, using `MKLROOT`; adapt `FC`, `BLAS_HOME`, and `BLASLIB` if your system uses a different compiler or BLAS implementation. Compile from the utility directory:

```bash
cd util/neoclass_program
make neoclassical_program
```

#### Run `neoclassical_program`

```bash
./neoclassical_program < input_neo_ITER
```

The program generates `neoclass_coef.dat`, which contains the normalized poloidal flux $\psi_N$, the neoclassical friction coefficient $\mu_{i,\mathrm{neo}}$, and the neoclassical heat-flux coefficient $k_i$ in JOREK units.

### 1.3. Use the coefficient profiles in JOREK

Add the following parameters to the JOREK input file:

```fortran
neo      = .t.
neo_file = 'neoclass_coef.dat'
```

## 2. Alternative: constant neoclassical coefficients

Instead of using profiles, constant values can be used for $\mu_{i,\mathrm{neo}}$ and $k_i$. In this case, set typical coefficients in the JOREK input file:

```fortran
neo           = .t.
aki_neo_const = -1.
amu_neo_const =  1.e-5
```
