---
title: "Run"
nav_order: 3
parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Running JOREK

## Run the Code

- Use the namelist input file as standard input for the binary:

```bash
./jorek_model303 < ./input
```

- For OpenMP, e.g. with 16 threads, set the following before running the code:

```bash
export OMP_NUM_THREADS=16
```

- Run with MPI parallelization, e.g. with 4 processes:
  - The number of MPI processes must be a multiple of `(n_tor+1)/2` [see hard-coded parameters](compilation_and_libraries/hard-coded_parameters.md)

```bash
mpirun -n 4 ./jorek_model303 < input
```

## Namelist Input File

- JOREK simulations are controlled by input parameters that are set in a namelist input file of the following form:

```fortran
&in1
  restart=.false.
  nstep_n=0
  tstep_n=1.e3
  ...
/
```

- **Some input parameters are model-specific.** Check `models/modelXXX/initialise_parameters.f90` for details.
- **Default values are assigned to the input parameters** by `models/preset_parameters.f90`.
- The input parameters are read by MPI thread 0 and sent to all other MPI threads via `trunk/communication/broadcast_phys.f90`
- **Example input files** can be found in the repository in the directories `trunk/namelist/modelXXX/`.
- All input parameters are written out automatically at the start of a simulation. This allows checking that the correct parameters are used, especially in case default values are used.

### General Parameters

| Input parameter | Explanation | Comment |
|----------|--------|--------|
| `restart` | Restart an interrupted simulation | Allows to restart a JOREK simulation from the `jorek_restart.rst` file. Default: `.false.` |
| `rst_format` | Use this format number in the restart file. Valid options: 0, 1, >1 | See `export_restart.f90` |
| `regrid` | Realign the grid on flux surfaces | Never used |
| `nout` | Output a restart file every `nout` steps | |

### Physics Parameters

| Input parameter | Explanation | Comment |
|----------|--------|--------|
| `tauIC` | Normalised inverse cyclotron frequency (extended MHD) | |
| `eta` | Resistivity | |
| `visco(_par)` | Viscosity of the perpendicular(parallel) velocity | |
| `D_par` | Parallel particle diffusion | = zero |
| `D_perp(i)` | Perpendicular particle diffusion $D_{base} = D_1 \left(( 1 - D_2) + \frac{D_2}{2} \left(1 - \tanh \left(\frac{\psi - D_5}{D_4} \right) \right)\right)$ if model < 300 then $D = D_{base}$ else $D = D_{base} + \frac{D_6 D_2}{2} \left(1 - \tanh \left( \frac{D_5 + D_3 -\psi}{D_4}  \right) \right)$ | in the centre, (1-D(2)) in the pedestal, width of the pedestal region, width of transition, start of pedestal, diffusion outside of pedestal respectively (1-6) |
| `ZK_par` | Parallel heat diffusion | |
| `ZK_perp(i)` | Perpendicular heat diffusion | See particle diffusion |
| `F0` | R*B_phi | |
| `FF_0` | FF' in the plasma centre | |
| `FF_1` | FF' at the plasma edge | |
| `FF_coef(i)` | `FF'(psi) = [(FF_0 - FF_1)*(1 + coef(1) * psi + coef(2)*psi**2 + coef(3)*psi**3) + coef(6)/cosh((psi - coef(7))/coef(8))**2 / (2.d0 * coef(8)) / (psi_1-psi_0) ]*0.5*(1- tanh((psi-coef(5))/coef(4)))  +  FF_1` | |
| `heatsource` | Amplitude of energy source | heating, in normalized units |
| `particlesource` | Amplitude of particle source | |
| `T_0, T_1` | Initial temperature profile, on axis (0) and at boundary | |
| `T_coef(i)` | `T(psi) = (T_0 - T_1)*(1 + coef(1) * psi + coef(2)*psi**2 + coef(3)*psi**3)  *0.5*(1- tanh((psi-coef(5))/coef(4)))  +  T_1` | |
| `rho_(0,1,coef(i))` | Initial density profile (see temperature) | |
| `V_(0,1,coef(i))` | Initial toroidal velocity (rotation frequency if `normalized_velocity_profile = .t.`) profile (see temperature) | |
| `central_density` | Central number density | Used in normalization of other quantities |
| `central_mass` | Central average ion mass | Deuterium mass if not defined |

### Numerics Parameters

| Input parameter | Explanation | Comment |
|----------|--------|--------|
| `tstep_n` | Size of the time step $\Delta t$ | Up to seven values possible (corresponding to several `nstep_n` values) |
| `nstep_n` | Number of time steps | Up to seven values possible (corresponding to several `tstep_n` values) |
| `tstep` | Timestep when not using `nstep_n` | |
| `nstep` | Number of `tstep` steps to make | |
| `eta_num` | Numerical hyper-resistivity | It is recommended to use small values according to `eta_num \lesssim eta^2` and perform a scan of `eta_num` to ensure it does not influence results. |
| `visco(_par)_num` | Numerical hyper-viscosity | It is recommended to use small values according to `visco_num \lesssim visco^2` and perform a scan of `visco_num` to ensure it does not influence results. |
| `D_perp_num` | Hyper-diffusion (particles, numerical stabilization) | |
| `ZK_perp_num` | Hyper-diffusion (heat, numerical stabilization) | |
| `tgnum(i)` | Factor for Taylor-Galerkin stabilization in equation i | See information about using the stabilization. |
| `iter_precon` | redo factorisation if gmres iterations > `iter_precon` | GMRES parameters |
| `gmres_m` | gmres dimension | |
| `gmres_4` | scaling of error | see source code |
| `gmres_max_iter` | maximum number of gmres iterations | |
| `gmres_tol` | convergence when error < `gmres_tol` | |

### Geometry Parameters

| Input parameter | Explanation | Comment |
|----------|--------|--------|
| `mf` | Number of harmonics in `fbnd` | if `mf=0` use ellipticity etc |
| `fbnd(i)` | Radius of plasma boundary (Fourier series) in polar angle | has `mf` terms |
| `ellip` | Ellipticity of initial polar grid | |
| `tria_u/l` | Upper/Lower triangularity | |
| `quad_u/l` | Upper/Lower quadrangularity | |
| `xampl` | Analytic definition of `psi(theta)` on the boundary | amplitude |
| `xsig` | | |
| `xwidth` | width of perturbation `psi(theta)` in `theta` | |
| `xtheta` | position of x-point in angle | |
| `xleft` | Shifts the radial position of the plasma | |
| `xshift` | Shifts the vertical plasma position | |
| `xpoint` | Plasma has x-point | |
| `R_begin` | Square Bezier grid minimum R | |
| `R_end` | Square Bezier grid maximum R | |
| `Z_begin` | Square Bezier grid minimum Z | |
| `Z_end` | Square Bezier grid maximum Z | |
| `R_geo` | R-coordinate of the centre of the polar grid | |
| `Z_geo` | Z-coordinate of the centre of the polar grid | |
| `amin` | The minor radius of the polar grid | |
| `n_R` | The number of horizontal nodes in the square Bezier grid | Leads to `n_R-1` elements in that direction |
| `n_Z` | The number of vertical nodes in the square Bezier grid | Leads to `n_Z-1` elements |
| `n_radial` | The number of radial nodes in the polar Bezier grid (grid 1) | |
| `n_pol` | The number of poloidal nodes in the polar Bezier grid (grid 1) | |
| `n_flux` | The number of radial nodes on the flux aligned grid inside the plasma boundary (grid 2) | |
| `n_tht` | The number of poloidal nodes (grid 2) | |
| `n_open` | The number of radial nodes, from magnetic axis to separatrix (grid 2) | |
| `n_leg` | The number of nodes along the separatrix legs, from x-point to divertor wall (grid 2) | |
| `n_private` | The number of radial nodes in the private region, region below the separatrix legs (grid 2) | |
