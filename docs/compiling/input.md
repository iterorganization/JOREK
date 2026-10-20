---
title: "List of input parameters"
nav_order: 4
parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

<!--
  AUTO-GENERATED FILE - DO NOT EDIT BY HAND.
  Regenerate by running ./util/parameter-overview.py from the repo root.
-->

<style>
  /* Target the tables directly rather than a `.params-table` class: kramdown's
     `{: .params-table}` IAL does not reliably attach to very large tables.
     All column widths use `em` (no percentages): mixing `%` text columns with
     `em` numeric columns under `table-layout: fixed` + `min-width: 100%` makes
     Firefox resolve the percentages against the wrong basis and collapse the
     text columns. With uniform `em` widths both engines distribute leftover
     space proportionally, keeping the description column widest. Scoped to
     .main-content so it only affects this page's parameter tables. */
  .main-content table { table-layout: fixed; width: 100%; font-size: 0.9em; }
  .main-content table th, .main-content table td {
      word-wrap: break-word; overflow-wrap: break-word; vertical-align: top; padding: 4px 6px;
  }
  .main-content table th:nth-child(1), .main-content table td:nth-child(1) { width: 10em; }
  .main-content table th:nth-child(2), .main-content table td:nth-child(2) { width: 8em; }
  .main-content table th:nth-child(3), .main-content table td:nth-child(3) { width: 30em; }
  .main-content table th:nth-child(n+4), .main-content table td:nth-child(n+4) {
      width: 2.2em; text-align: center;
  }
</style>

# List of input parameters

> **⚠️ Auto-generated file — do not edit manually.**
> This page is generated from the Fortran sources. Any manual changes will be
> overwritten the next time the generator runs. To update it, run
> `./util/parameter-overview.py` from the repository root and commit the result.

## phys_module

| parameter | default | description | 180 | 183 | 199 | 600 | 710 | 711 | 712 | 750 |
|---|---|---|---|---|---|---|---|---|---|---|
| **eta** | 1.d-5 | Resistivity at plasma cener (normalized) | x | x | x | x | x | x | x | x |
| **eta_ohmic** | 0. | Resistivity at core for the Ohmic heating term | x | x | x | x | x | x | x | x |
| **eta_T_dependent** | .true. | Resistivity dependent on temperature? Otherwise constant | x | x | x | x | x | x | x | x |
| **eta_coul_log_dep** | .true. | Resistivity dependent on variations of the Coulomb logarithm? |  |  |  | x |  |  |  | x |
| **T_max_eta** | 1.d99 | Temperature above which the resistivity is truncated (use with care; only for numerical reasons) |  |  |  | x |  |  |  |  |
| **T_max_eta_ohm** | 1.d99 | Temperature above which the resistivity used in the Ohmic heating term is truncated (use with care; only for numerical reasons) |  |  |  | x |  |  |  |  |
| **T_max_visco** | 1.d99 | Temperature above which the viscosity is truncated; It is aimed for keeping the Prandtl number constant when T_max_eta is activated. |  |  |  | x |  |  |  |  |
| **visco** | 1.d-5 | Viscosity at plasma center (normalized) | x | x | x | x | x | x | x | x |
| **visco_heating** | 0. | Viscosity used in the perpendicular viscous heating term |  |  |  | x |  |  |  |  |
| **visco_T_dependent** | .true. | Viscosity dependent on temperature? Otherwise constant. | x | x | x | x | x | x | x | x |
| **visco_old_setup** | .false. | If true, the old perp. viscosity treatment is used for compatibility (old visco depends on R^2) |  |  |  | x |  |  |  |  |
| **visco_par** | 1.d-5 | Cross B-field viscosity acting on parallel flow (normalized) | x | x | x | x | x | x | x | x |
| **visco_par_par** | 0. | B-field Parallel viscosity acting on parallel flow (normalized) |  | x |  | x |  |  |  |  |
| **visco_par_heating** | 0. | Parallel viscosity used in the parallel viscous heating term (normalized) |  |  |  | x |  |  |  |  |
| **TiTe_ratio** | 0.5 | ratio to set ion and electron temperature from T (in model 180): Ti=TiTe_ratio*T; Te=(1.0-TiTe_ratio)*T | x |  |  |  |  |  |  |  |
| **F0** | 10. | Determines fixed toroidal magnetic field: $ B_\phi = F_0/R $ | x | x | x | x | x | x | x | x |
| **central_density** | 1. | particle density at the magnetic axis (in units of $10^{20} m^{-3}$) | x | x | x | x | x | x | x | x |
| **central_mass** | 2.01410177811 | average ion mass in atomic mass units (constant in time and space, including electron mass) | x | x | x | x | x | x | x | x |
| **gamma** | 5. / 3. | ratio of specific heat (typically 5/3) |  |  |  |  | x | x | x | x |
| **tauIC** | 0. | Scaling factor for diamagnetic terms (see [[diamag\|diamagnetic]]) | x | x | x | x | x | x | x | x |
| **Wdia** | .false. | Include diamagnetic flows in viscosity terms? (see [[wdia\|here]]) |  |  |  | x |  |  |  |  |
| **gamma_sheath** | 4.5 | sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead! | x | x | x | x | x | x | x | x |
| **gamma_stangeby** | -1.d99 | Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices) |  |  |  | x | x | x | x | x |
| **gamma_sheath_e** | 3.00 | sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead! |  |  |  | x |  | x | x | x |
| **gamma_e_stangeby** | -1.d99 | Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices) |  |  |  | x |  |  |  |  |
| **gamma_sheath_i** | -1.11d-1 | sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead! |  |  |  | x |  | x | x | x |
| **gamma_i_stangeby** | -1.d99 | Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices) |  |  |  | x |  |  |  |  |
| **density_reflection** | 0. | density reflection coeefficient on open fieldlines | x | x | x | x |  |  |  |  |
| **neutral_reflection** | 0. | reflection coefficient of ions into neutrals (model500) |  |  |  | x |  |  | x | x |
| **loop_voltage** | 0. | Apply a loop voltage at the boundary of the computational domain (in V; works only for fixed boundary) |  |  |  | x |  |  |  |  |
| **old_deuterium_atomic** | .false. | use old fit to calculate atomic coefficients for D (ionization, recombination, radiation), otherwise a better fit is used |  |  |  | x |  |  |  |  |
| **deuterium_adas** | .false. | use OPEN ADAS to calculate ionization, recombination and radiation coeffients for deuterium |  |  |  | x |  |  |  |  |
| **deuterium_adas_1e20** | .false. | use OPEN ADAS with fixed density=1e20 to calculate ionization, recombination and radiation coeffients for deuterium |  |  |  | x |  |  |  |  |
| **mach_one_bnd_integral** | .false. |  |  |  |  | x |  |  |  |  |
| **vpar_smoothing** | .false. | apply a smoothing function to smooth jumps in Vpar at B.n=0 |  |  |  | x |  |  |  |  |
| **vpar_smoothing_coef** | 0.01, 0., 0. | coefficients for the smoothing profile of the parallel velocity |  |  |  | x |  |  |  |  |
| **min_sheath_angle** | 1. | For sheath boundary conditions: Minimum incident angle for heat and particle fluxes (in degrees) |  |  |  | x |  |  |  |  |
| **mode** |  | Toroidal mode number corresponding to the JOREK modes, e.g., for n_period=8 and n_tor=3, mode(:)=0,8,8 | x | x | x | x | x | x | x | x |
| **nout** | 9999999 | Output a restart file every nout timesteps | x | x | x | x | x | x | x | x |
| **nout_projection** | -1 | Output particle projection every nout_projection timesteps (only for diagnostics) |  |  |  | x | x | x | x | x |
| **xcase** | LOWER_XPOINT | 1->LowerXpoint. 2->UpperXpoint. 3->doubleNull | x | x | x | x | x | x | x | x |
| **forceSDN** | .false. | Force a symmetric double null, within the accuracy of SDN_threshold |  |  | x | x | x | x | x | x |
| **SDN_threshold** | 1.d-4 | threshold, in absolute psi, for a symmetric-double-null grid construction | x | x | x | x | x | x | x | x |
| **rst_format** | 0 | 0 == old format, 1 == new format for restart file | x | x | x | x |  |  |  |  |
| **restart** | .false. | Restart a code run from the restart file jorek_restart.h5? | x | x | x | x | x | x | x | x |
| **regrid** | .false. | Re-generate the flux-aligned grid (does not work currently)? | x | x | x | x | x | x | x | x |
| **regrid_from_rz** | .false. | Re-generate the flux-aligned grid from an rz equilibrium | x | x | x | x | x | x | x | x |
| **xpoint** | .false. | X-point plasma or not? see also xcase | x | x | x | x | x | x | x | x |
| **Z_xpoint_limit** | -0.4 0.4 | Search the lower X-point in the region Z < Z_xpoint_limit(1) and the upper X-point in the region Z > Z_xpoint_limit(2) |  |  |  | x | x | x | x | x |
| **xpoint_search_tries** | 500 | The number of candidate elements to check for being the element containing the upper or lower X-point. |  |  | x | x | x | x | x | x |
| **bootstrap** | .false. | Evolve the Bootstrap current consistently with time? | x | x | x | x |  |  |  |  |
| **bootstrap_psin_cutoff** | 0.9995 |  | x | x | x | x |  |  |  |  |
| **refinement** | .false. | Use mesh refinement? (not presently available) | x | x | x | x | x | x | x | x |
| **force_central_node** | .true. | Force all nodes in the center to have the same values in flux aligned grids or independent values? | x | x | x | x | x |  |  | x |
| **fix_axis_nodes** | .false. | Fix t-derivative and cross st-derivative on axis to avoid noise | x | x | x | x | x | x | x | x |
| **treat_axis** | .false. | > Flag for chosing grid axis treatment (see grids/mod_axis_treatment.f90) | x | x | x | x | x | x | x | x |
| **bc_natural_flux** | .false. | boundary conditions for flux surface boundaries (2 and 3) |  |  |  |  | x | x | x | x |
| **bc_natural_open** | .false. | use natural boundary conditions on the open fieldlines | x | x | x | x | x | x | x | x |
| **produce_live_data** | .true. | Write data 'macroscopic_vars.dat' during the code run allowing to use plot_live_data.sh? | x | x | x | x | x | x | x | x |
| **grid_to_wall** | .false. | extend the grid to a physical wall | x | x | x | x | x | x | x | x |
| **RZ_grid_inside_wall** | .false. | build the rectangular grid inside first wall |  |  |  | x | x | x | x | x |
| **RZ_grid_jump_thres** | 0.85 | threshold to change R-resolution as RZ-grid gets sqeezed by limiter contour |  |  |  | x | x | x | x | x |
| **manipulate_psi_map** | 0. 99. 99. 0.1 0.1 | Option to manipulate Psi_boundary for the initial grid | x | x | x | x | x | x | x | x |
| **adaptive_time** | .false. | (presently not useful) | x | x | x | x | x | x | x | x |
| **equil** | .true. | compute equilibrium | x | x | x | x | x | x | x | x |
| **no_mach1_bc** | .false. | Never apply Mach-1 BCs |  |  |  | x | x | x | x | x |
| **Mach1_openBC** | .true. | Full-MHD: Apply Mach-1 BCs inside mod_boundary_matrix_open.f90 (or mod_boundary_conditions.f90) |  |  |  |  | x | x | x | x |
| **Mach1_fix_B** | .true. | Full-MHD: Use the initial magnetic field for Mach1 BCs on targets, ie. without AR and AZ variations |  |  |  |  | x |  |  | x |
| **export_polar_boundary** | .false. | Option to export boundary.txt even in the case of a polar boundary. |  |  |  | x |  |  |  |  |
| **eta_ARAZ_const** | 0. | Use uniform resistivity for AR and AZ equations, used only if eta_ARAZ_on=.false. |  |  |  |  | x | x | x | x |
| **eta_ARAZ_on** | .true. | Full-MHD: to switch on/off resistive terms for AR and AZ equations |  |  |  |  | x | x | x | x |
| **eta_ARAZ_simple** | .false. | Full-MHD: remove the Fprof dependence of Bphi in the resistive terms for AR and AZ (which should be compensated by current source anyway) |  |  |  |  | x | x | x | x |
| **tauIC_ARAZ_on** | .true. | Full-MHD: to switch on/off diamagnetic terms for AR and AZ equations |  |  |  |  | x | x | x | x |
| **bench_without_plot** | .false. | if .true., do not produce certain output plots (e.g., for benchmarking) | x | x | x | x | x | x | x | x |
| **gmres** |  | Use iterative GMRES solver | x | x | x | x | x | x | x | x |
| **gmres_max_iter** | 200 | Maximum number of GMRES iterations | x | x | x | x | x | x | x | x |
| **keep_n0_const** | .false. | Perform a linear run where the equilibrium quantities (i_tor=1) do not change with time? |  |  | x | x | x | x | x | x |
| **linear_run** | .false. | Same as keep_n0_const, to be replaced soon by true linear run where modes are independent |  |  | x | x | x | x | x | x |
| **export_for_nemec** | .false. | Export equilibrium information for the NEMEC code? | x | x | x | x | x | x | x | x |
| **export_aux_node_list** | .true. | Include the aux_node_list for particle projections in the restart files |  |  | x | x | x | x | x | x |
| **use_murge** | .false. | (Deprecated, Cannot be used any more) | x | x | x | x | x | x | x | x |
| **use_murge_element** | .false. | (Deprecated, Cannot be used any more) | x | x | x | x | x | x | x | x |
| **use_BLR_compression** | .false. | Use Block-Low-Rank (BLR) compression in MUMPS / PaStiX 6 solvers | x | x | x | x | x | x | x | x |
| **epsilon_BLR** | 0. | Accuracy of BLR compression | x | x | x | x | x | x | x | x |
| **just_in_time_BLR** | .true. | Use Just-in-time strategy for BLR compression (speed optimized) | x | x | x | x | x | x | x | x |
| **write_ps** | .true. | Write postscript file at the end of the run | x | x | x | x | x | x | x | x |
| **use_mumps** | .false. | Use Mumps solver | x | x | x | x | x | x | x | x |
| **use_pastix** | .false. | Use Pastix solver | x | x | x | x | x | x | x | x |
| **use_strumpack** | .true. | Use Strumpack solver | x | x | x | x | x | x | x | x |
| **use_mumps_eq** | .false. | Use Mumps equilibrium solver | x | x | x | x | x | x | x | x |
| **use_pastix_eq** | .false. | Use Pastix equilibrium solver | x | x | x | x | x | x | x | x |
| **use_strumpack_eq** | .false. | Use Strumpack equilibrium solver | x | x | x | x | x | x | x | x |
| **use_mumps_prj** | .true. | Use Mumps projection solver | x | x | x | x | x | x | x | x |
| **use_pastix_prj** | .false. | Use Pastix projection solver | x | x | x | x | x | x | x | x |
| **use_strumpack_prj** | .false. | Use Strumpack projection solver | x | x | x | x | x | x | x | x |
| **use_wsmp** | .false. | Use WSMP solver | x | x | x | x |  |  |  |  |
| **centralize_harm_mat** | .true. | Centralize harmonic matrices on toridal master ranks; switch for STRUMPACK solver | x | x | x | x | x | x | x | x |
| **mumps_ordering** | 7 | MUMPS ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7 | x | x | x | x | x | x | x | x |
| **pastix_maxthrd** | 1024 | maximum number of threads used by pastix solver (could be beneficial to use the reduced number) | x | x | x | x | x | x | x | x |
| **pastix_pivot** |  | Pastix epsilon for magnitude control (pivot threshold) | x | x | x | x |  |  |  |  |
| **use_newton** | .false. | Use inexact Newton method |  |  | x | x | x | x | x |  |
| **maxNewton** | 20 | maximum number of Newton iterations |  |  | x | x | x | x | x |  |
| **gamma_Newton** | 0.5 | Newton gamma-parameter: gmres_tol = gamma_Newton*(normRHScurrent/normRHSprevious)**alpha_Newton |  |  | x | x | x | x | x |  |
| **alpha_Newton** | 2. | Newton alpha-parameter: gmres_tol = gamma_Newton*(normRHScurrent/normRHSprevious)**alpha_Newton |  |  | x | x | x | x | x |  |
| **strumpack_matching** | .false. | Perform maximum-diagonal-product reordering algorithm in STRUMPACK solver (improves direct solver, but use matrix centralization) |  |  | x | x | x | x | x | x |
| **bcs** | see [its documentation page](/JOREK/howto/physics_options/boundary_conditions/choose_boundary_conditions.html) |  |  |  |  |  |  |  |  |  |
| **n_limiter** | 0 | Number of limiter points | x | x | x | x | x | x | x | x |
| **R_limiter** | 0. | R-positions of the limiter points | x | x | x | x | x | x | x | x |
| **Z_limiter** | 0. | Z-positions of the limiter points | x | x | x | x | x | x | x | x |
| **first_target_point** |  |  | x | x | x | x | x | x | x | x |
| **last_target_point** |  |  | x | x | x | x | x | x | x | x |
| **gvec_grid_import** | .false. | Generate grid fourier representation with GVEC | x | x | x |  |  |  |  |  |
| **extended_boundary** | .false. | Choose if extended boundary conditions (Biot-Savart version) should be used, default (false) is grad_chi with Dommaschk potentials | x |  |  |  |  |  |  |  |
| **j_cutoff_rcoord** | 99.0 | Radial location from which the current is set to zero as it approaches the boundary - rcoord corresponds to the normalised toroidal flux | x |  |  |  |  |  |  |  |
| **j_cutoff_sig** | 0.025 | Radial width over which the current is ramped down to zero towards the boundary | x |  |  |  |  |  |  |  |
| **eqdsk_psi_fact** | 1. | multiply eqdsk psi by factor for grid_inside_wall |  |  |  | x | x | x | x | x |
| **extend_existing_grid** | .false. | Add patches to existing grid from restart file |  |  |  | x | x | x | x | x |
| **n_wall_blocks** | 0 | Number of blocks |  |  |  | x | x | x | x | x |
| **corner_block** | 0 | =1 for a corner block ("left" side will also be wall-aligned) |  |  |  | x | x | x | x | x |
| **n_ext_block** | 0 | Number of 'radial' grid points from the outermost flux surface to wall) |  |  |  | x | x | x | x | x |
| **n_ext_equidistant** | .false. | if true, radial spacing of grid points will be equidistant (not adapted) |  |  |  | x | x | x | x | x |
| **n_block_points_left** | 0 | Number of points on left side of block |  |  |  | x | x | x | x | x |
| **R_block_points_left** | 0. | R-positions of points on left side of block |  |  |  | x | x | x | x | x |
| **Z_block_points_left** | 0. | Z-positions of points on left side of block |  |  |  | x | x | x | x | x |
| **n_block_points_right** | 0 | Number of points on left side of block |  |  |  | x | x | x | x | x |
| **R_block_points_right** | 0. | R-positions of points on left side of block |  |  |  | x | x | x | x | x |
| **Z_block_points_right** | 0. | Z-positions of points on left side of block |  |  |  | x | x | x | x | x |
| **use_simple_bnd_types** | .false. | convert Stan's bnd_types to Guido's bnd_types |  |  |  | x | x | x | x | x |
| **xampl** | 0. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **xwidth** | 0. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **xsig** | 1. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **xtheta** | 0. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **xshift** | 0. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **xleft** | 0. | Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition) | x | x | x | x | x | x | x | x |
| **particlesource** | 1.e-5 | Particle source amplitude | x | x | x | x | x | x | x | x |
| **particlesource_psin** | 1.0 | Position around which the source is ramped down | x | x | x | x | x | x | x | x |
| **particlesource_sig** | 0.1 | Width over which the source is ramped down | x | x | x | x | x | x | x | x |
| **particlesource_gauss** | 0. | Additional Gaussian particle source amplitude | x | x | x | x | x | x | x | x |
| **particlesource_gauss_psin** | 0.9 | Position around which Gaussian source is set | x | x | x | x | x | x | x | x |
| **particlesource_gauss_sig** | 0.1 | Width over which Gaussian source is set | x | x | x | x | x | x | x | x |
| **edgeparticlesource** | 0. | Edge particle source amplitude | x | x | x | x | x | x | x | x |
| **edgeparticlesource_psin** | 0.98 | Position around which the edge particle source is located | x | x | x | x | x | x | x | x |
| **edgeparticlesource_sig** | 0.01 | Width over which edge particle source extends | x | x | x | x | x | x | x | x |
| **neutral_line_source** | 0. | neutral inflow source |  |  |  | x |  |  | x | x |
| **neutral_line_R_start** | 1.d20 | neutral inflow source (starting point of line source) |  |  |  | x |  |  | x | x |
| **neutral_line_Z_start** | 1.d20 | neutral inflow source |  |  |  | x |  |  | x | x |
| **neutral_line_R_end** | 2.d20 | neutral inflow source (end point of line source) |  |  |  | x |  |  | x | x |
| **neutral_line_Z_end** | 2.d20 | neutral inflow source |  |  |  | x |  |  | x | x |
| **heatsource** | 1.e-7 | Heat source amplitude | x | x | x | x | x | x | x | x |
| **heatsource_e** | 0.5e-7 | Electron heat source amplitude | x | x |  | x |  | x | x | x |
| **heatsource_i** | 0.5e-7 | Ion heat source amplitude | x | x |  | x |  | x | x | x |
| **heatsource_psin** | 1.0 | Position around which the source is ramped down | x | x | x | x | x | x | x | x |
| **heatsource_sig** | 0.1 | Width over which the source is ramped down | x | x | x | x | x | x | x | x |
| **heatsource_e_psin** | 1.0 | Position around which the electron source is ramped down | x |  |  | x |  |  |  |  |
| **heatsource_e_sig** | 0.1 | Width over which the electron source is ramped down | x |  |  | x |  |  |  |  |
| **heatsource_i_psin** | 1.0 | Position around which the ion source is ramped down | x |  |  | x |  |  |  |  |
| **heatsource_i_sig** | 0.1 | Width over which the ion source is ramped down | x |  |  | x |  |  |  |  |
| **heatsource_gauss** | 0. | Additional Gaussian heat source amplitude | x | x | x | x | x | x | x | x |
| **heatsource_gauss_psin** | 0.9 | Position around which Gaussian source is located | x | x | x | x | x | x | x | x |
| **heatsource_gauss_sig** | 0.1 | Width over which Gaussian source extends | x | x | x | x | x | x | x | x |
| **heatsource_gauss_e** | 0. | Gaussian heat source for electrons | x | x |  | x |  | x | x | x |
| **heatsource_gauss_i** | 0. | Gaussian heat source for ions | x | x |  | x |  | x | x | x |
| **heatsource_gauss_e_psin** | 0.9 | Position around which electrons Gaussian source is located | x | x |  | x |  |  |  |  |
| **heatsource_gauss_e_sig** | 0.1 | Width over which electrons Gaussian source extends | x | x |  | x |  |  |  |  |
| **heatsource_gauss_i_psin** | 0.9 | Position around which ions Gaussian source is located | x | x |  | x |  |  |  |  |
| **heatsource_gauss_i_sig** | 0.1 | Width over which ions Gaussian source extends | x | x |  | x |  |  |  |  |
| **constant_imp_source** | 0. | Adds a constant impurity source |  |  |  | x |  |  |  |  |
| **eta_num** | 0. |  | x | x | x | x | x | x | x | x |
| **visco_num** | 0. |  | x | x | x | x | x | x | x | x |
| **visco_par_num** | 0. |  | x | x | x | x | x | x | x | x |
| **Dn_perp_num** | 0. |  |  |  |  | x |  |  | x | x |
| **maintain_profiles** | .false. | Add artificial sources to maintain initial rho and T profiles |  | x |  |  |  |  |  |  |
| **use_sc** | .false. | Use shock-capturing stabilization |  |  |  | x |  |  |  | x |
| **D_perp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **D_par_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **Dn_pol_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **Dn_p_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **D_perp_imp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **D_par_imp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_perp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_par_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_i_perp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_i_par_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_e_perp_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **ZK_e_par_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **visco_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **visco_par_sc_num** | 0. |  |  |  |  | x |  |  |  | x |
| **eta_num_T_dependent** | .false. | Hyper-resistivity dependent on temperature? Otherwise constant. |  |  |  | x |  |  |  |  |
| **eta_num_psin_dependent** | .false. | Give profile for Hyper-resistivity as function of \psi_N? Useful for 2D current flattening |  |  |  | x |  |  |  |  |
| **eta_num_prof** | 0. 0.8 0.03 | Coefficients to specify \psi_N profile for hyper-resistivity |  |  |  | x |  |  |  |  |
| **visco_num_T_dependent** | .false. |  |  |  |  | x |  |  |  |  |
| **add_sources_in_sc** | .false. | Whether to add effect of sources in shock-capturing stabilization or not |  |  |  | x |  |  |  | x |
| **use_vms** | .false. | Use VMS stabilization in model 750 only |  |  |  |  |  |  |  | x |
| **vms_coeff_AR** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_AZ** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_A3** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_UR** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_UZ** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_Up** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_T** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_Te** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_Ti** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_rho** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_rhon** | 0. |  |  |  |  |  |  |  |  | x |
| **vms_coeff_rhoimp** | 0. |  |  |  |  |  |  |  |  | x |
| **tstep** | 1. | Size of the timesteps ($ \Delta t $) | x | x | x | x | x | x | x | x |
| **tstep_n** | 1. | Alternative to tstep: Up to ten values may be given | x | x | x | x | x | x | x | x |
| **nstep** | 0 | Number of timesteps to perform | x | x | x | x | x | x | x | x |
| **nstep_n** | 0 | Alternative to nstep: Up to ten values may be given | x | x | x | x | x | x | x | x |
| **time_evol_scheme** | 'Crank-Nicholson' | Time evolution scheme to use (see [[time-integration\|time_integration]]) | x | x | x | x | x | x | x | x |
| **time_evol_theta** |  | Time evolution parameter theta (see [[time-integration\|time_integration]]) |  |  |  |  | x | x | x | x |
| **time_evol_zeta** |  | Time evolution parameter zeta (see [[time-integration\|time_integration]]) |  |  |  |  | x | x | x | x |
| **rst_hdf5** | 1 | Write hdf5 restart files if set to 1 | x | x | x | x | x | x | x | x |
| **rst_hdf5_version** |  | Write which version of hdf5 files? | x | x | x | x | x | x | x | x |
| **tokamak_device** | 'none' | Name of the tokamak device we are simulating | x | x | x | x | x | x | x | x |
| **amin** | 1. | Minor radius for polar grid construction, set to 1 if boundary is specified with R,Z points | x | x | x | x | x | x | x | x |
| **ellip** | 1. | Ellipticity of polar grid (see analytical definition in phys_module.f90) | x | x | x | x | x | x | x | x |
| **tria_u** | 0. | Upper triangularity of polar grid (see analytical definition in phys_module.f90) | x | x | x | x | x | x | x | x |
| **tria_l** | 0. | Lower triangularity of polar grid (see analytical definition in phys_module.f90) | x | x | x | x | x | x | x | x |
| **quad_u** | 0. | Upper quadrangularity of polar grid (see analytical definition in phys_module.f90) | x | x | x | x | x | x | x | x |
| **quad_l** | 0. | Lower quadrangularity of polar grid (see analytical definition in phys_module.f90) | x | x | x | x | x | x | x | x |
| **mf** | 2 | Number of entries in fbnd and fpsi | x | x | x | x | x | x | x | x |
| **fbnd** | 0.;   fbnd(1)  = 2. | Fourier expansion of boundary | x | x | x | x | x | x | x | x |
| **fpsi** |  | Fourier expansion of the poloidal flux at the boundary | x | x | x | x | x | x | x | x |
| **n_boundary** | 0 | Number of points in R_boundary, Z_boundary, psi_boundary. | x | x | x | x | x | x | x | x |
| **R_boundary** | 0. | Numerical R values defining the boundary | x | x | x | x | x | x | x | x |
| **Z_boundary** | 0. | Numerical Z values defining the boundary | x | x | x | x | x | x | x | x |
| **psi_boundary** | 0. | Numerical values giving the poloidal flux at the boundary | x | x | x | x | x | x | x | x |
| **n_pfc** | 0 | Number of coils, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **Rmin_pfc** | 0. | Minimum R of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **Rmax_pfc** | 0. | Maximum R of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **Zmin_pfc** | 0. | Minimum Z of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **Zmax_pfc** | 0. | Maximum Z of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **current_pfc** | 0. | Current density in the coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall\|JOREK-STARWALL]] | x | x | x | x |  |  |  |  |
| **n_jropes** | 0 | Number of ropes, |  |  |  | x |  |  |  |  |
| **R_jropes** | 0. | R centre of rope |  |  |  | x |  |  |  |  |
| **Z_jropes** | 0. | Z centre of rope |  |  |  | x |  |  |  |  |
| **w_jropes** | 0. | width of rope |  |  |  | x |  |  |  |  |
| **current_jropes** | 0. | Current inside the rope |  |  |  | x |  |  |  |  |
| **pellet_amplitude** | 0. | amplitude of density source (when pellet modelled as density source) | x | x | x | x | x | x | x | x |
| **pellet_R** | 3.8 | major radius position pellet | x | x | x | x | x | x | x | x |
| **pellet_Z** | 0.0 | Z position pellet | x | x | x | x | x | x | x | x |
| **pellet_phi** | 1.57 | width of the pellet cloud (density source) in toroidal angle | x | x | x | x | x | x | x | x |
| **pellet_ellipse** | 5. | the ellipticity of the pellet source | x | x | x | x |  |  |  |  |
| **pellet_radius** | 0.08 | radius of the simulation pellet | x | x | x | x | x | x | x | x |
| **pellet_sig** | 0.02 | width of smoothing of density source (arctan( (r-pellet_radius)/pellet_sig) ) | x | x | x | x | x | x | x | x |
| **pellet_length** | 0.785 | width of smoothing of density source in toroidal angle | x | x | x | x | x | x | x | x |
| **pellet_theta** | 0. | orientation of the pellet ellipse | x | x | x | x |  |  |  |  |
| **pellet_psi** | 1.0 | pellet_width in poloidal flux | x | x | x | x | x | x | x | x |
| **pellet_delta_psi** | 999. | width of smoothing in poloidal flux | x | x | x | x | x | x | x | x |
| **pellet_velocity_R** | 0. | pellet velocity component radial direction | x | x | x | x |  |  |  |  |
| **pellet_velocity_Z** | 0. | pellet velocity component Z direction | x | x | x | x |  |  |  |  |
| **pellet_density** | 5.985d8 | pellet atom number density (in units $10^{20} m^{-3}$) | x | x | x | x |  |  |  | x |
| **pellet_density_bg** | 5.958d8 | background species pellet atom number density (in units $10^{20} m^{-3}$) |  |  |  | x |  |  |  | x |
| **pellet_particles** | 0. | the number of particles in the pellet (in units of $10^{20}$) | x | x | x | x |  |  | x | x |
| **use_pellet** | .false. |  | x | x | x | x |  |  | x | x |
| **t_ns** | 2.d3 | MGI onset time (JOREK units) |  |  |  | x |  |  | x | x |
| **ns_amplitude** | 0. | Amplitude of gas source |  |  |  | x |  |  | x | x |
| **ns_R** | 3.2 | R position of gas source |  |  |  | x |  |  | x | x |
| **ns_Z** | 1.5 | Z position of gas source |  |  |  | x |  |  | x | x |
| **ns_phi** | 1.57 | Phi position of gas source |  |  |  | x |  |  | x | x |
| **ns_radius** | 0.08 | Poloidal radius of gas source |  |  |  | x |  |  | x | x |
| **ns_deltaphi** | 0.5 | Toroidal extension of gas source |  |  |  | x |  |  | x | x |
| **ns_delta_minor_rad** | 0. | Extension of gas source in the minor radial direction (if greater than 0.) |  |  |  | x |  |  | x | x |
| **drift_distance** | 0. | Shift the R position of the neutral deposition outward by drift_distance (in meters) for plasmoid drift |  |  |  | x |  |  | x | x |
| **energy_teleported** | 0. | Energy (in eV) teleported per atom to consider plasmoid drift effects |  |  |  | x |  |  | x | x |
| **imp_type** | ' ' | Type of injected material or background impurity species: Argon, neon, ... |  |  |  | x |  |  |  | x |
| **use_imp_adas** | .true. | Use open adas to calculate ionization, recombination and radiation coeffients for impurities |  |  |  | x |  |  |  |  |
| **JET_MGI** | .false. | Switch to use a JET-like MGI |  |  |  | x |  |  | x | x |
| **ASDEX_MGI** | .false. | Switch to use an ASDEX-like MGI |  |  |  | x |  |  | x | x |
| **V_Dmv** | 9.75d-4 | Volume of the DMV reservoir |  |  |  | x |  |  | x | x |
| **P_Dmv** |  | Pressure in the DMV reservoir (bar) |  |  |  | x |  |  | x | x |
| **A_Dmv** | 1.77d-2 | Cross sectional area of DMV (Disruption mitigation valve) pipe |  |  |  | x |  |  | x | x |
| **K_Dmv** | 4.d-2 | Correction parameter describing the gas expansion near the pipe orifice |  |  |  | x |  |  | x | x |
| **L_tube** | 0. | Pipe length |  |  |  | x |  |  | x | x |
| **ksi_ion** | 1.84d-24 | Energy cost of each ionization, ksi_ion / mu_0 / (gamma-1) / e = 13.7 eV |  |  |  | x |  |  | x | x |
| **delta_n_convection** | 0 | Switch to activate the convection term for neutrals (at the plasma velocity) |  |  |  | x |  |  | x | x |
| **nimp_bg** | 0. | Density of background impurities (in $m^{-3}$) |  |  |  | x |  |  | x | x |
| **index_main_imp** | 0 | Index of the main impurity species (in imp_type and nimp_bg) solved with continuity equation |  |  |  | x |  |  |  | x |
| **using_spi** | .false. | This determines whether to use SPI or traditional MGI; see [[spi_tutorial\|SPI Tutorial]] |  |  |  | x |  |  | x | x |
| **spi_Vel_Rref** | 0.0 | Reference velocity of pellet center along R upon injection |  |  |  | x |  |  | x | x |
| **spi_Vel_Zref** | 0.0 | Reference velocity of pellet center along Z upon injection |  |  |  | x |  |  | x | x |
| **spi_Vel_RxZref** | 0.0 | Reference velocity of pellet center along RxZ direction upon injection |  |  |  | x |  |  | x | x |
| **spi_quantity** | 0.0 | Total injected atom number for impurity SPI |  |  |  | x |  |  | x | x |
| **spi_quantity_bg** | 0.0 | Total injected atom number for background species SPI |  |  |  | x |  |  |  | x |
| **ns_radius_ratio** | 1.4 | We are assuming a constant ratio between the radius of NG clouds |  |  |  | x |  |  | x | x |
| **spi_Vel_diff** | 0.0 | The velocity difference from the reference velocity |  |  |  | x |  |  | x | x |
| **spi_angle** | 0.0 | The vertex angle of spi spreading in terms of rad |  |  |  | x |  |  | x | x |
| **spi_L_inj** | 0.25 | Distance between SPI nozzle and ns_R, ns_Z, ns_phi |  |  |  | x |  |  | x | x |
| **spi_L_inj_diff** | 0.0 | The position difference with respect to the point (ns_R, ns_Z, ns_phi) |  |  |  | x |  |  | x | x |
| **tor_frequency** | 0.0 | The rigid body rotation frequency |  |  |  | x |  |  | x | x |
| **ns_radius_min** | 8.d-2 | This defines the minimum radius of neutral cloud for numerical reasons (in m) |  |  |  | x |  |  | x | x |
| **spi_abl_history_old** | .false. | If this is .t., convert the old spi_abl_history format to the new one upon restart. |  |  |  | x |  |  |  |  |
| **n_spi** | 0 1 | Number of shattered fragment injected for each injection |  |  |  | x |  |  | x | x |
| **n_inj** | 1 | Number of injections |  |  |  | x |  |  | x | x |
| **spi_abl_model** | -1 | Determine which type of ablation model is used. |  |  |  | x |  |  | x | x |
| **spi_rnd_seed** | 0 | Random seed array used for the generation of the SPI velocity spread |  |  |  | x |  |  | x | x |
| **spi_shard_file** | 'none' | The name of the shard size file |  |  |  | x |  |  |  | x |
| **spi_plume_file** | 'none' | The name of the shard information datafile (array) |  |  |  | x |  |  |  | x |
| **spi_plume_hdf5** | .false. | if 'spi_plume_file' is in HDF5format? |  |  |  | x |  |  |  |  |
| **spi_abl_mag_reduction** | .false. | Whether to use the magnetic reduction effect described in Eq.(27) of Nucl. Fusion 60 066027 |  |  |  | x |  |  |  |  |
| **n_adas** | 1 | Number of species to be traced by ADAS |  |  |  | x |  |  |  | x |
| **spi_tor_rot** | .false. | Flag to turn on a rigid body toroidal plasma rotation for SPI |  |  |  | x |  |  | x | x |
| **spi_num_vol** | .true. | Flag to turn on numerical integration of the gas source volumes from SPI |  |  |  | x |  |  | x | x |
| **adas_dir** | ' ' | The directory of ADAS data file to be read |  |  |  | x |  |  |  | x |
| **output_prad_phi** | .false. | Output Prad(phi) into a file using integrals_3D |  |  |  | x |  |  |  | x |
| **amix** | 0. | Mix Poisson solution with previous one with a given factor | x | x | x | x | x | x | x | x |
| **equil_accuracy** | 1.d-6 | Tolerance of the convergence for the fix-boundary equilibrium | x | x | x | x | x | x | x | x |
| **axis_srch_radius** | 99. | Magnetic axis will be searched inside a circle with this radius | x | x | x | x | x | x | x | x |
| **delta_psi_GS** | 10000. | Expected psi_bnd - psi_axis for the final equilibrium |  |  | x | x | x | x | x | x |
| **newton_GS_fixbnd** | .false. | Newton instead of Picard iterations for fixed-boundary equilibria? |  |  | x | x | x | x | x | x |
| **newton_GS_freebnd** | .true. | Newton instead of Picard iterations for free-boundary equilibria? |  |  | x | x | x | x | x | x |
| **freeboundary_equil** | .false. | use a free or fixed boundary equilibrium? ([[jorek-starwall\|JOREK-STARWALL]]) | x | x | x | x | x | x | x | x |
| **freeboundary** | .false. | use free or fixed boundary conditions in time-evolution? ([[jorek-starwall\|JOREK-STARWALL]]) | x | x | x | x | x | x | x | x |
| **resistive_wall** | .false. | use a resistive or ideal wall? ([[jorek-starwall\|JOREK-STARWALL]]) | x | x | x | x | x | x | x | x |
| **freeb_equil_iterate_area** | .false. | iterate to a target area during freeboundary equilibrium limiter cases [[jorek-starwall-faqs\|jorek_starwall]] | x | x | x | x | x | x | x | x |
| **amix_freeb** | 0.85 | choose amix for freeboundary equilibrium | x | x | x | x | x | x | x | x |
| **equil_accuracy_freeb** | 1.d-6 | Tolerance of the convergence for the freeboundary equilibrium | x | x | x | x | x | x | x | x |
| **freeb_change_indices** | .true. | Exchange grid node indices to parallelize boundary integral | x | x | x | x | x | x | x | x |
| **n_R** | 0 | Number of grid points in R-direction (for rectangular grid) (see also [[grids#tutorials\|here]]) | x | x | x | x | x | x | x | x |
| **n_Z** | 0 | Number of grid points in Z-direction (for rectangular grid) | x | x | x | x | x | x | x | x |
| **R_begin** | -0.1 | Left boundary of grid in R-direction (for rectangular grid) | x | x | x | x | x | x | x | x |
| **R_end** | 0.1 | Right boundary of grid in R-direction (for rectangular grid) | x | x | x | x | x | x | x | x |
| **Z_begin** | -0.1 | Lower boundary of grid in Z-direction (for rectangular grid) | x | x | x | x | x | x | x | x |
| **Z_end** | 0.1 | Upper boundary of grid in Z-direction (for rectangular grid) | x | x | x | x | x | x | x | x |
| **rect_grid_vac_psi** | 0. | Use a vacuum psi-bnd condition for squared-grid, ie. (rect_grid_vac_psi * R**2) |  |  |  | x |  |  |  |  |
| **force_horizontal_Xline** | .false. | Force the grid line through Xpoint to be horizontal (instead of perp. to line between Xpoint and axis) | x | x | x | x | x | x | x | x |
| **n_radial** | 11 | Number of radial grid points (for polar grid) (see also [[grids\|here]]) | x | x | x | x | x | x | x | x |
| **n_pol** | 16 | Number of poloidal grid points (for polar grid) | x | x | x | x | x | x | x | x |
| **R_geo** | 10. | Center of the grid (for polar grid) | x | x | x | x | x | x | x | x |
| **Z_geo** | 0. | Center of the grid (for polar grid) | x | x | x | x | x | x | x | x |
| **psi_axis_init** | -0.1 | Initial guess for Psi at the magnetic axis (for polar grid) | x | x | x | x | x | x | x | x |
| **XR_r** | 999. | Psi_N position of radial grid accumulation (two positions) (for polar grid) (also used for R-position in square-grid) | x | x | x | x | x | x | x | x |
| **SIG_r** | 999. | Width of grid accumulation (two positions) (for polar grid) (also used for R-width in square-grid) | x | x | x | x | x | x | x | x |
| **XR_tht** | 999. | Position of poloidal grid accumulation (0...1, two positions) (for polar grid) | x | x | x | x | x | x | x | x |
| **SIG_tht** | 999. | Width of grid accumulation (two positions) (for polar grid) | x | x | x | x | x | x | x | x |
| **XR_z** | 999. | Z-position of square grid accumulation (two positions) (for square grid) |  |  |  | x |  |  |  |  |
| **SIG_z** | 999. | Z-Width of grid accumulation (two positions) (for square grid) |  |  |  | x |  |  |  |  |
| **bgf_r** | 0.7 |  |  |  |  | x |  |  |  |  |
| **bgf_z** | 0.7 | Background for meshac distribution for R-Z accumulation |  |  |  | x |  |  |  |  |
| **bgf_rpolar** | 0.6 |  |  |  | x | x | x | x | x | x |
| **bgf_tht** | 0.6 | Background for meshac distribution for R-theta accumulation |  |  | x | x | x | x | x | x |
| **n_flux** | 11 | Number of radial grid points (for flux-aligned grid) (see also [[grids#tutorials\|here]]) | x | x | x | x | x | x | x | x |
| **n_tht** | 16 | Number of poloidal grid points (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **xr1** | 9999. | Grid accumulation parameter (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **xr2** | 99999. | Grid accumulation parameter (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **sig1** | 9999. | Grid accumulation parameter (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **sig2** | 99999. | Grid accumulation parameter (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **m_pol_bc** | 1 | Number of poloidal modes for Psi boundary condition in stellarator | x |  |  |  |  |  |  |  |
| **i_plane_rtree** | 1 | The poloidal plane in a stellarator on which the RTree is to be built (RZ_minmax refers to this plane) | x | x |  |  |  |  |  |  |
| **n_open** | 5 | Number of 'radial' grid points in the open flux region - between the two separatrices if double-null | x | x | x | x | x | x | x | x |
| **n_outer** | 0 | Number of 'radial' grid points in the open flux region on the outer side (LFS) if double-null | x | x | x | x | x | x | x | x |
| **n_inner** | 0 | Number of 'radial' grid points in the open flux region on the inner side (HFS) if double-null | x | x | x | x | x | x | x | x |
| **n_private** | 5 | Number of 'radial' grid points in the private flux region at the bottom | x | x | x | x | x | x | x | x |
| **n_leg** | 5 | Number of 'poloidal' grid points along the divertor legs at the bottom | x | x | x | x | x | x | x | x |
| **n_leg_out** | 0 | Number of 'poloidal' grid points along the divertor legs at the bottom on the LFS |  |  |  | x | x | x | x | x |
| **n_up_priv** | 0 | Number of 'radial' grid points in the private flux region at the top (upper Xpoint or double-null) | x | x | x | x | x | x | x | x |
| **n_up_leg** | 0 | Number of 'poloidal' grid points along the divertor legs at the top (upper Xpoint or double-null) | x | x | x | x | x | x | x | x |
| **n_up_leg_out** | 0 | Number of 'poloidal' grid points along the divertor legs on the top on the LFS (upper Xpoint or double-null) |  |  |  | x | x | x | x | x |
| **n_ext** | 0 | Number of 'radial' grid points from the outermost flux surface to wall) | x | x | x | x | x | x | x | x |
| **n_tht_equidistant** | .false. | switch on to get an equidistant poloidal distribution of elements in the core of the grid (psi<0.5) |  |  |  | x | x | x | x | x |
| **xr_closed** | 1.0, 9999., 9999. | Location for grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_closed** | 0.1, 9999., 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_open** | 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_outer** | 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_inner** | 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_private** | 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_up_priv** | 0.1 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_theta** | 0.03 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_theta_up** | 999. | Width with grid accumulation (for flux-aligned grid; only valid for double-null) |  |  | x | x | x | x | x | x |
| **SIG_leg_0** | 0.05 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_leg_1** | 0.2 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_up_leg_0** | 0.05 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **SIG_up_leg_1** | 0.2 | Width with grid accumulation (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **dPSI_open** | 0.11 | Delta Psi grid extends into the open flux region (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **dPSI_outer** | 0.11 | Delta Psi grid extends into the open flux region (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **dPSI_inner** | 0.11 | Delta Psi grid extends into the open flux region (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **dPSI_private** | 0.03 | Delta Psi grid extends into the private flux region (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **dPSI_up_priv** | 0.03 | Delta Psi grid extends into the private flux region (for flux-aligned grid) | x | x | x | x | x | x | x | x |
| **D_perp** | 1.d-5, 0., 0., 99., 99. | Coefficients for perpendicular particle diffusion profile | x | x | x | x | x | x | x | x |
| **D_par** | 0. | Parallel particle diffusion (usually not useful) | x | x | x | x | x | x | x | x |
| **V_pinch_gauss** | 0. | Amplitude of Gaussian inward pinch velocity profile for background fluid (rho only). |  |  |  | x |  |  |  |  |
| **V_pinch_psin** | 0. | Centre of V_pinch Gaussian in normalised poloidal flux (psin). |  |  |  | x |  |  |  |  |
| **V_pinch_sig** | 1. | Width (sigma) of V_pinch Gaussian in psin units. |  |  |  | x |  |  |  |  |
| **D_perp_imp** | 1.d-5, 0., 0., 99., 99. | Coefficients for perpendicular imp particle diffusion profile |  |  |  | x |  |  |  | x |
| **D_par_imp** | 0. | Parallel impurity particle diffusion (usually not useful) |  |  |  | x |  |  |  | x |
| **ZK_perp** | 1.d-5, 0., 0., 99., 99. | Coefficients for perpendicular heat diffusion profile | x | x | x | x | x | x | x | x |
| **ZK_par** | 1. | Parallel heat diffusion value in the plasma center | x | x | x | x | x | x | x | x |
| **ZK_par_max** | 1.d20 | Do not use larger parallel heat diffusion values for numerical reasons | x | x | x | x | x | x | x | x |
| **T_min_ZKpar** | -1.d12 | Do not use smaller parallel heat diffusion values below this MHD temperature (Ti+Te); JOREK units |  |  |  | x |  |  |  |  |
| **Ti_min_ZKpar** | -1.d12 | Do not use smaller parallel heat diffusion values below Ti; JOREK units |  |  |  | x |  |  |  |  |
| **Te_min_ZKpar** | -1.d12 | Do not use smaller parallel heat diffusion values below Te; JOREK units |  |  |  | x |  |  |  |  |
| **ZK_i_perp** | 1.d-5, 0., 0., 99., 99. | Coefficients for perpendicular ion heat diffusion profile |  | x |  | x |  | x | x | x |
| **ZK_e_perp** | 1.d-5, 0., 0., 99., 99. | Coefficients for perpendicular electron heat diffusion profile |  | x |  | x |  | x | x | x |
| **ZK_i_par** | 1. | Ion parallel heat diffusion coefficient in the plasma center |  | x |  | x |  | x | x | x |
| **ZK_e_par** | 1. | Electron parallel heat diffusion coefficient in the plasma center |  | x |  | x |  | x | x | x |
| **D_neutral_x** | 1.d-5 | Neutral particle diffusivity in R-direction |  |  |  | x |  |  | x | x |
| **D_neutral_y** | 1.d-5 | Neutral particle diffusivity in Z-direction |  |  |  | x |  |  | x | x |
| **D_neutral_p** | 1.d-5 | Neutral particle diffusivity in phi-direction |  |  |  | x |  |  | x | x |
| **ZKpar_T_dependent** | .true. | Use a temperature dependent parallel heat diffusivity | x | x | x | x | x | x | x | x |
| **d_perp_file** | 'none' | ASCII file with perpendicular particle diffusion profile | x | x | x | x |  |  |  |  |
| **zk_perp_file** | 'none' | ASCII file with perpendicular heat diffusion profile | x | x | x | x |  |  |  |  |
| **zk_e_perp_file** | 'none' | ASCII file with perpendicular electron heat diffusion profile |  |  |  | x |  |  |  |  |
| **zk_i_perp_file** | 'none' | ASCII file wtih perpendicular ion heat diffusion profile |  |  |  | x |  |  |  |  |
| **v_pinch_file** | 'none' | ASCII file with inward pinch velocity profile (psin, V_pinch columns) |  |  |  | x |  |  |  |  |
| **rho_0** | 1. | Central normalized density (usually 1) | x | x | x | x | x | x | x | x |
| **rho_1** | 1. | SOL normalized density | x | x | x | x | x | x | x | x |
| **rho_coef** | 0.;  rho_coef(1) =  0. | Density profile coefficients | x | x | x | x | x | x | x | x |
| **rho_file** | 'none' | ASCII file the density profile is read from. | x | x | x | x | x | x | x | x |
| **T_0** | 1.d-6 | Central normalized temperature | x | x | x | x | x | x | x | x |
| **T_1** | 1.d-8 | SOL normalized temperature | x | x | x | x | x | x | x | x |
| **T_coef** | 0.;  T_coef(1)   = -1. | Temperature profile coefficients | x | x | x | x | x | x | x | x |
| **Ti_0** | 5.d-7 | Central ion normalized temperature | x | x |  | x |  | x | x | x |
| **Ti_1** | 5.d-9 | SOL ion normalized temperature | x | x |  | x |  | x | x | x |
| **Ti_coef** | 0.;  Ti_coef(1)  = -1. | Ion temperature profile coefficients | x | x |  | x |  | x | x | x |
| **Te_0** | 5.d-7 | Central ion normalized temperature | x | x |  | x |  | x | x | x |
| **Te_1** | 5.d-9 | SOL ion normalized temperature | x | x |  | x |  | x | x | x |
| **Te_coef** | 0.;  Te_coef(1)  = -1. | Ion temperature profile coefficients | x | x |  | x |  | x | x | x |
| **T_file** | 'none' | ASCII file the temperature profile is read from. | x | x | x | x | x | x | x | x |
| **Ti_file** | 'none' | ASCII file the ion temperature profile is read from. |  |  |  | x |  | x | x | x |
| **Te_file** | 'none' | ASCII file the electron temperature profile is read from. |  |  |  | x |  | x | x | x |
| **rhon_0** | 0. | Central value for the initial normalized neutral density |  |  |  | x |  |  |  |  |
| **rhon_1** | 0. | SOL value for the initial normalized neutral density |  |  |  | x |  |  |  |  |
| **rhon_coef** | 0. 0.01 0.01 | Coefficients for the intitial neutral density profile |  |  |  | x |  |  |  |  |
| **Fprofile_file** | 'none' | ASCII file the Fprofile is read from. |  |  |  |  | x | x | x | x |
| **phi_0** | 0. | Central background potential; (usually 1) | x | x |  |  |  |  |  |  |
| **phi_1** | 0. | Edge background potential | x | x |  |  |  |  |  |  |
| **phi_coef** | 0.;  phi_coef(1) =  0.; phi_coef(4) = 1. | potential profile coefficients | x | x |  |  |  |  |  |  |
| **phi_file** | 'none' | ASCII file the potential profile is read from. | x | x |  |  |  |  |  |  |
| **nu_phi_source** | 0. | Friction coefficient of the n=0 background potential profile source term (>~ visco) |  | x |  |  |  |  |  |  |
| **FF_0** | 1. | FF' value in the plasma center | x | x | x | x | x | x | x | x |
| **FF_1** | 0. | FF' value in the SOL | x | x | x | x | x | x | x | x |
| **FF_coef** | 0.;  FF_coef(1)  = -1. | Coefficients for FF' profile | x | x | x | x | x | x | x | x |
| **ffprime_file** | 'none' | ASCII file the FF' profile is read from. | x | x | x | x | x | x | x | x |
| **NEO** | .false. | If .true. neoclassical effects are considered, (see [[neo\|here]]) | x | x | x | x | x | x | x | x |
| **neo_file** | 'none' | ASCII file the aki and amu profiles is read from. | x | x | x | x | x | x | x | x |
| **aki_neo_const** | 0. | if ( (NEO) .and. (neo_file=='none') ), this constant value is used for aki_neo | x | x | x | x | x | x | x | x |
| **amu_neo_const** | 0. | if ( (NEO) .and. (neo_file=='none') ), this constant value is used for amu_neo | x | x | x | x | x | x | x | x |
| **output_bnd_elements** | .false. | If .true., writes bnd nodes and bnd elements in files 'boundary_nodes.dat' and 'boundary_elements.dat' | x | x | x | x | x | x | x | x |
| **RMP_on** | .false. | Activates RMPs on boundary if .true. (the old version without STARWALL) | x | x | x | x |  |  |  |  |
| **RMP_psi_cos_file** | 'none' | ASCII file the profiles of psi_RMP_cos and derivatives are read from | x | x | x | x |  |  |  |  |
| **RMP_psi_sin_file** | 'none' | ASCII file the profiles of psi_RMP_sin and derivatives are read from | x | x | x | x |  |  |  |  |
| **RMP_growth_rate** | 0.011 |  | x | x | x | x |  |  |  |  |
| **RMP_ramp_up_time** | 1000 | parameters for time dependence of psi_RMP: Sigmoid f(t)= 1/ (1 + exp(-RMP_growth_rate*(t-RMP_ramp_up_time/2) )) | x | x | x | x |  |  |  |  |
| **RMP_har_cos** | 2 |  | x | x | x | x |  |  |  |  |
| **RMP_har_sin** | 3 |  | x | x | x | x |  |  |  |  |
| **Number_RMP_harmonics** | 1 | Number_RMP_harmonics < N_RMP_max. If only one harmonic,  Number_RMP_harmonics=1, by default it's =1 in models/preset_parameters.f90 |  |  |  | x |  |  |  |  |
| **RMP_har_cos_spectrum** | RMP_har_cos | If only one harmonic,by default RMP_har_cos_spectrum(1)=RMP_har_cos; |  |  |  | x |  |  |  |  |
| **RMP_har_sin_spectrum** | RMP_har_sin | If only one harmonic,by default RMP_har_sin_spectrum(1)=RMP_har_sin |  |  |  | x |  |  |  |  |
| **V_0** | 0. | analytical parallel rotation profile -- central value | x | x | x | x | x | x | x | x |
| **V_1** | 0. | analytical parallel rotation profile -- SOL value | x | x | x | x | x | x | x | x |
| **V_coef** | 0., 0., 0., 0.1, 1.0 | analytical parallel rotation profile -- coefficients | x | x | x | x | x | x | x | x |
| **R_Z_psi_bnd_file** | 'none' | ASCII file for R_boundary,Z_boundary, psi_boundary, with n_boundary size. | x | x | x | x | x | x | x | x |
| **wall_file** | 'wall.txt' | ASCII file for external wall geometry, if n_ext is greater than zero. | x | x | x | x | x | x | x | x |
| **rot_file** | 'none' | ASCII file the parallel rotation profile is read from (see normalized_velocity_profile) | x | x | x | x |  |  |  |  |
| **normalized_velocity_profile** | .true. | if true, reads the normalized velocity profile as flux function, else Omega_tor is read as flux function. | x | x | x | x |  |  |  |  |
| **domm_file** | 'none' | Namelist file containing the coefficients for Dommaschk potentials | x | x |  |  |  |  |  |  |
| **iter_precon** | 10 | whenever the number of gmres iterations exceeds iter_precon, the preconditioning matrix is updated | x | x | x | x | x | x | x | x |
| **max_steps_noUpdate** | 10000000 | whenever the steps without preconditioning matrix update exceeds max_steps_noUpdate, the preconditioning matrix is updated | x | x | x | x | x | x | x | x |
| **gmres_m** | 20 | gmres restart parameter (dimension) | x | x | x | x | x | x | x | x |
| **gmres_4** | 1.d3 | see gmres manual (error ratio between preconditioned and non-preconditioned error) | x | x | x | x | x | x | x | x |
| **gmres_tol** | 1.d-8 | the tolerance for the gmres iterations to be seen as converged | x | x | x | x | x | x | x | x |
| **tgnum** | 0. | Coefficients for Taylor Galerkin stabilization for each equation separately | x | x | x |  |  |  |  |  |
| **tgnum_psi** | 0. | Same as previous line, but avoiding equation indexing for model families |  |  |  | x |  |  |  |  |
| **tgnum_u** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_zj** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_w** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_rho** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_T** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_Ti** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_Te** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_vpar** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_rhon** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_rhoimp** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_nre** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_AR** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_AZ** | 0. |  |  |  |  | x |  |  |  |  |
| **tgnum_A3** | 0. |  |  |  |  | x |  |  |  |  |
| **keep_current_prof** | .true. | Artificial current source to approximately keep the initial current profile, i.e., $\eta(j-j0)$? | x | x | x | x | x | x | x | x |
| **init_current_prof** | .false. | Initialize the current source from the current profile present | x | x |  |  |  |  |  |  |
| **D_prof_neg** | 1.d-5 | Particle diffusion coefficient in regions with negative background species density | x | x | x | x |  |  |  |  |
| **D_prof_neg_thresh** | 0. | D_prof_neg becomes effective if r0-rimp0 < D_prof_neg_thresh | x | x | x | x |  |  |  |  |
| **D_prof_imp_neg_thresh** | -1.d3 | D_prof_neg becomes effective if rimp0 < D_prof_imp_neg_thresh |  |  |  | x |  |  |  |  |
| **D_prof_tot_neg_thresh** | 0. | D_prof_neg becomes effective if r0 < D_prof_tot_neg_thresh |  |  |  | x |  |  |  |  |
| **ZK_prof_neg** | 1.d-5 | Perp. heat diffusion coefficient in regions with negative temperature | x | x | x | x |  |  |  |  |
| **ZK_par_neg** | 1.d-3 | Parallel diffusion coefficient in regions with negative temperature |  |  |  | x |  |  |  |  |
| **ZK_prof_neg_thresh** | 0. | ZK_prof_neg becomes effective if T < ZK_prof_neg_thresh | x | x | x | x |  |  |  |  |
| **ZK_par_neg_thresh** | 0. | ZK_par_neg becomes effective if T < ZK_par_neg_thresh |  |  |  | x |  |  |  |  |
| **ZK_e_prof_neg** | 1.d-5 | Perp. heat diffusion coefficient in regions with negative temperature |  |  |  | x |  |  |  |  |
| **ZK_e_par_neg** | 1.d-3 | Parallel diffusion coefficient in regions with negative temperature |  |  |  | x |  |  |  |  |
| **ZK_e_prof_neg_thresh** | 0. | ZK_e_prof_neg becomes effective if T < ZK_e_prof_neg_thresh |  |  |  | x |  |  |  |  |
| **ZK_e_par_neg_thresh** | 0. | ZK_e_par_neg becomes effective if T < ZK_e_par_neg_thresh |  |  |  | x |  |  |  |  |
| **ZK_i_prof_neg** | 1.d-5 | Perp. heat diffusion coefficient in regions with negative temperature |  |  |  | x |  |  |  |  |
| **ZK_i_par_neg** | 1.d-3 | Parallel diffusion coefficient in regions with negative temperature |  |  |  | x |  |  |  |  |
| **ZK_i_prof_neg_thresh** | 0. | ZK_i_prof_neg becomes effective if T < ZK_i_prof_neg_thresh |  |  |  | x |  |  |  |  |
| **ZK_i_par_neg_thresh** | 0. | ZK_i_par_neg becomes effective if T < ZK_i_par_neg_thresh |  |  |  | x |  |  |  |  |
| **D_imp_extra_neg_thresh** | -1.d3 | D_imp_extra_neg becomes effective if rho_imp < D_imp_extra_neg_thresh |  |  |  | x |  |  |  |  |
| **T_min** | 1.0d-20 | minimum temperature (limits on the temperature dependence of resistivity etc.) value in jorek units: 2.01d-5*central_density*Tmin_ev (preset central_density = 1, 20 eV) | x | x | x | x | x | x | x | x |
| **rho_min** | 1.0d-20 | minimum density |  |  | x | x | x | x | x | x |
| **ne_SI_min** | 1.d18 | minimum e density (in SI unit) below which we cut-off the radiation loss |  |  |  | x |  |  |  |  |
| **Te_eV_min** | 5. | minimum temperature (in eV) below which we cut-off the radiation loss |  |  |  | x |  |  |  | x |
| **rn0_min** | 1.d-8 | minimum impurity density (in JU) for radiation loss cut-off |  |  |  | x |  |  |  |  |
| **T_min_neg** | -1.d12 | minimum temperature,used for correcting negative values,in jorek units: 2.01d-5*central_density*Tmin_ev (preset central_density = 1, 20 eV) |  |  | x | x | x | x | x | x |
| **rho_min_neg** | -1.d12 | minimum density, used for correcting negative values |  |  | x | x | x | x | x | x |
| **implicit_heat_source** | 0. | Choose = 1.d0 to fully switch on the implicit heat source for numerical stabilization |  |  |  | x |  |  |  |  |
| **n_tor_fft_thresh** | 2 | If n_tor >= n_tor_fft_thresh, element_matrix_fft will be used | x | x | x | x | x | x | x | x |
| **corr_neg_temp_coef** | 0.5, 0.5 | Parameters used in models/corr_neg.f90 | x | x | x | x | x | x | x | x |
| **corr_neg_dens_coef** | 0.5, 0.5 | Parameters used in models/corr_neg.f90 | x | x | x | x | x | x | x | x |
| **thermalization** | .true. | If true turns on the ion-electron thermalization term |  |  |  | x |  | x | x | x |
| **zjz_0** | 0.1173 |  | x | x | x | x | x | x | x | x |
| **zjz_1** | 0.0 |  | x | x | x | x | x | x | x | x |
| **zj_coef** | 0.;  zj_coef(1)  = -1. |  | x | x | x | x | x | x | x | x |
{: .params-table}

## vacuum

| parameter | default | description | 180 | 183 | 199 | 600 | 710 | 711 | 712 | 750 |
|---|---|---|---|---|---|---|---|---|---|---|
| **CARIDDI_mode** | .false. | CARIDDI or STARWALL;True if CARIDDI input file |  |  | x | x | x | x | x | x |
| **vacuum_min** | .false. | Mode to minimalize memory consumption |  |  | x | x | x | x | x | x |
| **wall_resistivity_fact** | 1. | Scaling factor for the wall and coil resistivities specified in STARWALL | x | x | x | x |  |  |  |  |
| **wall_resistivity** | 0. | Resistivity of the external wall | x | x | x | x |  |  |  |  |
| **n_pf_coils** | 0 | number of poloidal field coils | x | x | x | x | x | x | x | x |
| **starwall_equil_coils** | .false. | specify wheter the equilibrium PF coils will be given by STARWALL or not | x | x | x | x | x | x | x | x |
| **find_pf_coil_currents** | .false. | search for optimal pf_coil currents to build a free-bnd equil? [[jorek-starwall-faqs\|fbnd_eq_FAQs]] | x | x | x | x | x | x | x | x |
| **psi_offset_freeb** | 0. | Allows to shift the value of psi by a global constant for freeb_equil (improves convergence) | x | x | x | x | x | x | x | x |
| **current_ref** | 1.d22 | Target total plasma current Ip for the feedback (FB) [[jorek-starwall-faqs\|fbnd_eq_FAQs]] | x | x | x | x | x | x | x | x |
| **FB_Ip_position** | 0.2 | Amplification factor for Ip feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **FB_Ip_integral** | 0.01 | Amplification factor for Ip feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **Z_axis_ref** | 1.d22 | Target magnetic axis vertical position (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **R_axis_ref** | -99. | Optional target magnetic axis radial position (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) |  |  | x | x | x | x | x | x |
| **FB_Zaxis_position** | 1. | Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **FB_Zaxis_derivative** | 0. | Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **FB_Zaxis_integral** | 0. | Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **cte_current_FB_fact** | -1d99 | Constant factor that scales FF'& T profiles before freebnd GS iterations (switches off current FB) |  |  | x | x | x | x | x | x |
| **start_VFB** | 10 | Iteration for starting vertical feedback (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **n_feedback_current** | 2 | Feedback will be performed each n_... iterations (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **n_feedback_vertical** | 1 | Feedback will be performed each n_... iterations (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **n_iter_freeb** | 900 | Number of iterations for freeboundary equilibirum (see [[jorek-starwall-faqs\|fbnd_eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **diag_coils** |  | see [[jorek-starwall-faqs\|jorek_starwall_FAQs]] | x | x | x | x | x | x | x | x |
| **rmp_coils** |  | see [[jorek-starwall-faqs\|jorek_starwall_FAQs]] | x | x | x | x | x | x | x | x |
| **voltage_coils** |  | not ready yet (see [[jorek-starwall-faqs\|jorek_starwall_FAQs]]) | x | x | x | x | x | x | x | x |
| **pf_coils** | 0 | see [[jorek-starwall-faqs\|jorek_starwall_FAQs]] | x | x | x | x | x | x | x | x |
| **vert_FB_amp** |  | Tune direction and magnitude of vert feedback for each poloidal field coil ([[jorek-starwall-faqs\|eq_FAQs]]) | x | x | x | x | x | x | x | x |
| **rad_FB_amp** |  | Tune direction and magnitude of vert feedback for each poloidal field coil ([[jorek-starwall-faqs\|eq_FAQs]]) |  |  | x | x | x | x | x | x |
| **vert_pos_file** |  |  | x | x | x | x | x |  |  |  |
| **start_VFB_ts** | 0. | start time of active VFB during simulation ([JOREK units]) | x | x | x | x | x |  |  |  |
| **vert_FB_amp_ts** | 0. | Amplitude and sign of vert feedback for each coil ([[jorek-starwall-faqs\|eq_FAQs]]);amplification factor (of PF coil) | x | x | x | x | x |  |  |  |
| **I_coils_max** | 1.d99 | Current limit of each coil ([Ampere]);Maximum absolute value for coils | x | x | x | x | x |  |  |  |
| **vert_FB_gain** | 0. | Gain parameters for vertical feedback controller;Proportional, derivative, integral gain of VFB controller | x | x | x | x | x |  |  |  |
| **vert_FB_tact** | 1.d-9 | Time interval between two controller actions ([JOREK units]);Tact of VFB controller | x | x | x | x | x |  |  |  |
{: .params-table}
