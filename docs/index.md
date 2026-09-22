---
title: "JOREK Wiki"
nav_order: 1
nav_exclude: true
layout: default
render_with_liquid: false
---

- **[JOREK Website](https://www.jorek.eu)** → [Send an e-mail](mailto:mhoelzl@ipp.mpg.de) if something is missing \| [Add to the gallery](https://www.jorek.eu)
- **[Information for New JOREK Users](jorek_access.md)** \| **[Code License](code_license.md)**
- **[Publication Rules](publication_rules.md)** \| **[JOREK Logo](jorek_logo.md)** \| **[jorek-preprint@jorek.eu](mailto:jorek-preprint@jorek.eu)**
- **[Frequently Asked Questions (FAQ)](faq.md)**
- **[Contributing to the Docs](contributing.md)**

### [Compiling and Running](compiling/cat_compiling.md)

- [Obtain the Source Code](compiling/get_jorek.md) \| [Compilation and Libraries](compiling/compilation_and_libraries.md) \| [Run](compiling/running.md)
- **[HPC & PC Configurations](compiling/compilation_and_libraries/hpc_pc_configurations.md)**
  - [Leonardo-CPU](compiling/compilation_and_libraries/hpc_pc_configurations/leonardo-cpu.md) \| [LUMI-CPU](compiling/compilation_and_libraries/hpc_pc_configurations/lumi-cpu.md) \| [Pitagora-CPU](compiling/compilation_and_libraries/hpc_pc_configurations/pitagora-cpu.md) \| [Pitagora-GPU](compiling/compilation_and_libraries/hpc_pc_configurations/pitagora-gpu.md) \| [IPP Garching](compiling/compilation_and_libraries/hpc_pc_configurations/run_at_ipp_garching.md) \| [ITER Cluster](compiling/compilation_and_libraries/hpc_pc_configurations/iter_cluster.md) \| [EUROfusion Gateway](compiling/compilation_and_libraries/hpc_pc_configurations/eurofusion_gateway.md) \| [TGCC-CEA](compiling/compilation_and_libraries/hpc_pc_configurations/tgcc-cea.md) \| [macOS](compiling/compilation_and_libraries/hpc_pc_configurations/macos.md)
- **[List of Input Parameters](compiling/input.md)**
- **[Diagnostics and post-processing](diagnostics/cat_diagnostics.md)**
  - [Diagnostics and Scripts](diagnostics/diagnostics.md) \| [JOREK-IMAS](diagnostics/jorek-imas.md)

### Code Development

- **[Development Workflow](code_development/development_workflow.md)**
- **[Regression Tests](code_development/nrt.md)**

### Physics Models

- **[Notation Conventions](physics/notation.md)**
- **[Normalization](physics/normalization.md)**
- **[Vector Identities](physics/vector-identities.md)**
- **[Coordinate Systems](physics/coordinates.md)**
- **[Base Fluid Models](physics/base_fluid_models/base_fluid_models.md)**
  - [Tokamak Reduced MHD](physics/base_fluid_models/reduced_mhd.md)
  - [Tokamak Full MHD](physics/base_fluid_models/full_mhd.md)
  - [Stellarator Reduced MHD](physics/base_fluid_models/jorek3d.md)
- **[Kinetic Particle Module](physics/particles.md)**
- **[Model Extensions](physics/model_extensions/model_extensions.md)**
  - [Free Boundary Extension (STARWALL, CARIDDI)](physics/model_extensions/freebound.md)
  - Impurities: [fluid model](physics/model_extensions/impurities_fluid.md) \| [marker model](physics/model_extensions/impurities_marker.md) \| [kinetic model](physics/model_extensions/impurities_kinetic.md)
    - [ADAS Atomic Data](physics/model_extensions/adas.md)
  - Neutrals: [fluid model](physics/model_extensions/neutrals_fluid.md) \| [kinetic model](physics/model_extensions/neutrals_kinetic.md)
  - REs: [fluid model](physics/model_extensions/runaway_fluid.md) \| [kinetic model](physics/model_extensions/runaway_kinetic.md)

### Numerics and Tools

- [Spatial Discretization](numerics/spatial-discretization.md) \| [Grids](numerics/grids.md)
- [Time Integration](numerics/time-integration.md)
- [Element Matrix FFT](numerics/element_matrix_fft.md)
- [Solver and Preconditioner](numerics/solver.md) \| [Sparse Matrix Format](numerics/sparse-matrix.md) 
- [Random-Number Generators](numerics/rngs.md)
- [HDF5 Tools](numerics/hdf5tools.md)

### Howto...

- Get started with **[JOREK](howto/running_jorek_for_the_first_time.md)**, **[JOREK-STARWALL](howto/running_jorek-starwall_for_the_first_time.md)**, **[JOREK-CARIDDI](howto/running_jorek-cariddi_for_the_first_time.md)** and **[JOREK Diagnostics](howto/introduction_to_jorek_diagnostics.md)**
- Set up a **[JOREK simulation grid](howto/wallgrid_tutorial.md)**
- Convert EFIT equilibrium data into JOREK input with **[eqdsk2jorek](howto/eqdsk2jorek.md)**
- Get **[D_perp and ZK_perp for stationary profiles](howto/physics_options/plasma_functions_and_coefficients/diffusion_coef.md)**
- Check **[energy conservation](diagnostics/energy_conservation.md)**
- Run with **[diamagnetic drift](howto/physics_options/physical_terms/diamag.md)** and **[neoclassical effects](howto/physics_options/physical_terms/neo.md)** and **[include diamagnetic drift in the viscosity term](howto/physics_options/physical_terms/wdia.md)**
- Run with **[Taylor-Galerkin Stabilization](howto/tgnum.md)** and **[correct negative densities / temperatures](howto/corr_neg.md)** (workaround)
- Run with **[RMPs](howto/physics_options/controllers_and_actuators/rmp.md)** (old boundary conditions, without STARWALL)
- Calculate **B** outside and inside the JOREK grid (**[jorek2_fields_xyz](howto/jorek2_fields_xyz.md)**)
- Calculate the **[total wall forces](howto/jorek2_wall_forces.md)** (needs JOREK-STARWALL)
- Run including **[Ohmic heating](howto/physics_options/physical_terms/ohmic_heating.md)**
- Run with **[Sheath heat-flux BC](howto/physics_options/boundary_conditions/sheath_heatflux_bc.md)**
- Run an **[MGI simulation](howto/mgi_tutorial.md)** and an **[SPI simulation](howto/spi_tutorial.md)**
- Set up **[Spitzer resistivity](howto/physics_options/plasma_functions_and_coefficients/spitzer_resistivity.md)** and anisotropic heat diffusion
- Create an **[X-point plasma from a limiter plasma](howto/x-point_from_limiter.md)**
- Run with **[mode groups in preconditioner](howto/mode_groups.md)**
- **[Choose boundary conditions](howto/physics_options/boundary_conditions/choose_boundary_conditions.md)**
- **[How to use shock capturing features](howto/shock_capturing.md)**
- **[Run stellarator simulations](howto/stellarator_setup.md)**
- **[Plot equation terms in VTK](howto/plot_rhs_terms.md)**
- **[Use phase space projections](howto/particles_phase_space.md)**
- **[Assess particle wall loads with the particle tracker](howto/particles_wall_load.md)**
- **[Assess fluid loads on 3D walls with field line tracing](howto/fluid_wall_load.md)**
- **[Generate Poincaré plots with the particle tracker](howto/particles_poincare.md)**
- **[Runaway electron physics in the particle tracker](howto/particles_runaways.md)**
- **[Use the controller module](howto/physics_options/controllers_and_actuators/using_controller_module.md)**
- Run with **[Non-linear time-stepping (Newton iterations)](howto/inexact_newton_solver.md)**
- **[Reconstruct how namelist input parameters were changed](howto/nml2h5.md)**
- **[Show time in Paraview plots](howto/showing_time_in_paraview.md)**
- **[Variable MultiScale (VMS) Stabilization in full MHD model 750](howto/vms.md)**
- Compile and run the **[Particle Fast Camera](howto/particle_fast_camera.md)**
- **[Run the stellarator model with the divertor region](howto/stellarator_with_divertor.md)**
- **[Run with kinetic neutrals and impurities](howto/ncs_ics_tutorial.md)**
- Run particle simulations with **[neutral neutral collisions](howto/neutral_neutral_collisions.md)**
- Activate the **[inward pinch term in the density equation](howto/physics_options/physical_terms/inward_pinch_term.md)**
- **[Computing parallel electric field $E_{\\vert\\vert}$ in JOREK](physics/base_fluid_models/RMHD/expressions_and_derivations/efield_in_jorek.md)**
- **[Helicity conservation](physics/base_fluid_models/RMHD/expressions_and_derivations/helicity_conservation.md)**

### Machines, Coordinates, Geometry, Synthetic Diagnostics, Reference Scenarios

- **[ITER](machines/iter.md)**
- **[JET](machines/jet.md)**
- **[ASDEX Upgrade](machines/asdex_upgrade.md)**
