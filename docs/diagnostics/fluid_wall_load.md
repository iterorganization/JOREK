---
title: "Assess Fluid Loads on 3D Walls"
nav_order: 3
parent: "Diagnostics and post-processing"
layout: default
render_with_liquid: false
---

# How to Assess Fluid Wall Loads with Field-Line Tracing

The Fortran code that calculates fluid loads on 3D walls is located in

```text
particles/examples/fluid_loads_on_3D_wall.f90
```

## What Does This Code Do?

- It calculates fluid fluxes (heat and current loads) on 3D thin walls discretized with triangles.
- It exports them in VTK format (`3D_wall_fluid_loads.vtk`).

## How to Compile and Run

1. Compile the program, similarly to `jorek2vtk`:

   ```bash
   make -j 8 fluid_loads_on_3D_wall
   ```

2. Copy the fluid restart file from which you want to calculate the heat fluxes to `jorek_restart.h5`, for example:

   ```bash
   cp jorek07000.h5 jorek_restart.h5
   ```

3. Execute the program using the usual JOREK input file:

   ```bash
   export OMP_NUM_THREADS=56 # The code only uses OpenMP parallelization (not MPI)
   ./fluid_loads_on_3D_wall < input_file_jorek
   ```

## Necessary Inputs

You must provide the following two HDF5 files in your simulation folder:

- **`wall_to_load.h5`**: The wall on which field-line tracing will be performed and the fluxes calculated.
- **`wall.h5`**: The wall with which field lines intersect. Normally, this includes the structures around the loaded wall as well as the loaded structure itself.

The format of these HDF5 files is explained in the [particle wall-load documentation](/JOREK/howto/particles_wall_load.html).

## Important Things to Know

- For field-line tracing to be effective, **the wall triangles must be within the JOREK domain**. Otherwise, those triangles are considered shadowed (not connected to the plasma), and there are no loads there. Choose your wall carefully.
- Be careful with the ordering of the wall triangles: the code assumes that their normals point towards the plasma.
- `l_par_min` is hardcoded. You may want to adapt it to your tokamak; you need to recompile after changing the code. You may also want to check the other hardcoded parameters near the beginning of the code.
- The code assumes that the fluxes are parallel to the magnetic field lines. Important perpendicular contributions are therefore not taken into account.
- The calculation of the parallel heat flux depends on the `qpar_tot` expression of `jorek2_postproc`, which you can inspect in `diagnostics/new_diag/mod_expression.f90`. It is the sum of the kinetic and thermal heat fluxes of the plasma ions and electrons, as shown in the [MHD sheath heat-flux documentation](/JOREK/howto/physics_options/boundary_conditions/sheath_heatflux_bc.html).

  **Check that the grid resolution is adequate for the chosen parallel and perpendicular heat-diffusion coefficients.** Otherwise, sharp boundary gradients—for example, those resulting from a parallel-conduction coefficient that is too small—may not be resolved, and the calculated parallel heat flux may be completely wrong.

  To verify this, check global energy conservation using the [energy-conservation diagnostics](./plot_live_data.sh.md#check-energy-conservation). When using sheath boundary conditions, it is usually more convenient to replace $q_\parallel$ with the employed $q_\parallel$ formula; this typically works much better. You may also want to add your own fluxes, such as particle or runaway-electron fluxes, to the code.

## Useful Resources

- [Scripts for time loops and calculating the temperature rise with a 1D equation](https://github.com/jorekart/scripts_javier/tree/main/tasks/loads_3d_walls)
