---
title: "Jorek2_connection_flux_aligned"
nav_order: 23
render_with_liquid: false
parent: "General"
---

## Jorek2_connection_flux_aligned

#### Overview
Several implementations of jorek2_connection exist in JOREK already, producing output for field lines in VTK format. This diagnostic is a  modification that produces flux surface averaged connection lengths - defined as the distance a field line travels before striking a given boundary, normally the wall. In such a way, it can be used to create radially averaged $\Psi_N-t$ plots to assess ergodicity and parallel transport.

#### Running the code
The code is run as

  ```
  mpirun -np X jorek2_connection_flux_aligned < ./input_file
  ```
 
#### Output

The output produces:
- connections.txt: Individual field line data for the connection length.
- strikes.txt: Final location of individual field lines in (R, Z, phi). The initial temperature and normalised poloidal flux are saved as a measure of the start point. The connection length is recorded as well as a boolean, in_domain, that defines whether the field line has struck the defined boundary.
- poinc_R-Z.dat/poinc_rho-theta.dat: Poincares based on field line tracing. By experience, it is better to use jorek2_poincare for poincares because the number of field lines necessary for statistically significant connection length data is too large for a good poincare plot.

#### Program settings
A connect.nml file in the run directory is used to set the parameters for the routine:
  &connect_params
   n_turns     = 150                          ! number of toroidal turns to follow a fieldline
   n_phi       = 1000                         ! number of steps per toroidal turn
   tol         = 1.d-6                        ! tolerance when stepping from element to element
   phi_start   = 0.0                          ! starting toroidal angle
   ntheta      = 200                          ! number of poloidal starting points on each flux surface
   n_rcoord       = 21                           ! number of radial starting points on each flux surface
   rcoord_range_min = 0.001                     ! minimum psi normalised
   rcoord_range_max = 0.999                     ! maximum psi normalised (< 0.999)
   rcoord_strike_bnd = 999.0                    ! Boundary to define strike points on (if psin>psin_max then the simulation boundary is used)
   /
